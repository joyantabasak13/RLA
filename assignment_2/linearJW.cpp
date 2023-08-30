#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iostream>
#include <map>
#include <pthread.h>
#include <set>
#include <string>
#include <sys/time.h>
#include <random>
#include <tuple>
#include <time.h>
#include <unistd.h>
#include <utility>
#include <vector>
#include <mutex>
#include <queue>
#include <memory>

// #include <bits/stdc++.h>
// #include <boost/algorithm/string.hpp>
// #include <boost/foreach.hpp>

using namespace std;

// Jaro distance variables
vector<string> recordVector;
int stringMaxLen;
const int alphabetSize = 36;
int* jwMatS1;
int* jwMatS2;
int* s1CharCount;
int* s2CharCount;
int* hash_s1;
int* hash_s2;
int* strLength;
int* recordIntArray;


// Best version till May 5th
// keep two tables of size alphabets * max possible string length; one for each string
// keep two array of alphabet size, each keeps count of characters in corresponding string
double jaroDistanceLinear(int ind1, int ind2){
	int len1 = strLength[ind1];
	int len2 = strLength[ind2];

	if (len1 == 0 || len2 == 0)
		return 0.0;

	int range = floor(max(len1, len2) / 2) - 1;

	int match = 0;
	int ptr1, ptr2;
	int count1, count2;
	int point1 = ind1*stringMaxLen;
	int point2 = ind2*stringMaxLen;

	for(int i=0; i<len1; i++) {
		jwMatS1[recordIntArray[point1+i] * stringMaxLen + s1CharCount[recordIntArray[point1+i]]] = i;
		s1CharCount[recordIntArray[point1+i]]++;
	}

	for(int i=0; i<len2; i++) {
		jwMatS2[recordIntArray[point2+i] * stringMaxLen + s2CharCount[recordIntArray[point2+i]]] = i;
		s2CharCount[recordIntArray[point2+i]]++;
	}

	for(int i=0; i<len1; i++) {
		int charInt = recordIntArray[point1+i];
		count1 = s1CharCount[charInt];
		count2 = s2CharCount[charInt];
		if((count1 != 0) && (count2 != 0)) {
			// Do linear pass to find matchs for charInt character 
			ptr1 = 0;
			ptr2 = 0;
			int offset = charInt * stringMaxLen;
			while((ptr1<count1) && (ptr2<count2))
			{
				if (abs(jwMatS1[offset + ptr1] - jwMatS2[offset + ptr2]) <= range)
				{
					match++;
					hash_s1[jwMatS1[offset + ptr1]] = 1;
					hash_s2[jwMatS2[offset + ptr2]] = 1;
					ptr1++;
					ptr2++;
				} else {
					if (jwMatS2[offset + ptr2] < jwMatS1[offset + ptr1])
					{
						ptr2++;
					} else
					{
						ptr1++;
					}
				}
			}
		}
		s1CharCount[charInt] = 0;
	}

	// memset gives better runtime for smaller alphabet sizes
	// memset(s2CharCount, 0, alphabetSize * sizeof(int));
	for(int i=0; i<len2; i++) {
		s2CharCount[recordIntArray[point2+i]] = 0;
	}

	// Linear pass for transposition calculation

	double t = 0.0;
	int k = 0;
	for (int i = 0; i < len1; i++)
	{
		if (hash_s1[i] == 1)
		{
			int j;
			for (j = k; j < len2; j++)
			{
				if (hash_s2[j] == 1)
				{
					k = j + 1;
					hash_s2[j] = 0;
					break;
				}
			}
			if (recordIntArray[point1+i] != recordIntArray[point2+j])
			{
				t++;
			}
		}
		hash_s1[i] = 0;
	}

	t = t / 2.0;

	return (((double)match) / ((double)len1) + ((double)match) / ((double)len2) + ((double)match - t) / ((double)match)) / 3.0;
}

// geeks for geeks
double jaroDistance_1(string &s1, string &s2)
{
	// If the strings are equal
	if (s1 == s2)
		return 1.0;

	// Length of two strings
	int len1 = s1.length();
	int len2 = s2.length();

	// Maximum distance upto which matching
	// is allowed
	int max_dist = floor(max(len1, len2) / 2) - 1;

	// Count of matches
	int match = 0;

	// Hash for matches
	std::vector<int> hash_s1(len1, 0);
	std::vector<int> hash_s2(len2, 0);

	// Traverse through the first string
	for (int i = 0; i < len1; i++)
	{

		// Check if there is any matches
		for (int j = max(0, i - max_dist);
			 j < min(len2, i + max_dist + 1); j++)

			// If there is a match
			if (s1[i] == s2[j] && hash_s2[j] == 0)
			{
				hash_s1[i] = 1;
				hash_s2[j] = 1;
				match++;
				break;
			}
	}

	// If there is no match
	if (match == 0)
		return 0.0;

	// Number of transpositions
	double t = 0;

	int point = 0;

	// Count number of occurrences
	// where two characters match but
	// there is a third matched character
	// in between the indices
	for (int i = 0; i < len1; i++)
		if (hash_s1[i])
		{

			// Find the next matched character
			// in second string
			while (hash_s2[point] == 0)
				point++;

			if (s1[i] != s2[point++])
				t++;
		}

	t /= 2;

	// cout << "Quad 1 Algo" << endl;
	// cout << "Match: " << match << endl;
	// cout << "Transposition: " << t << endl;
	// cout << endl;

	// Return the Jaro Similarity
	return (((double)match) / ((double)len1) + ((double)match) / ((double)len2) + ((double)match - t) / ((double)match)) / 3.0;
}

// code shared by Prof. Sahni
double jaroDistance_2(std::string &s1, std::string &s2)
{
	double m = 0;
	int low, high, range;
	int k = 0, numTrans = 0;

	// Exit early if either are empty
	if (s1.length() == 0 || s2.length() == 0)
	{
		return 0;
	}

	// Exit early if they're an exact match.
	if (s1 == s2)
	{
		return 1;
	}

	range = (std::max(s1.length(), s2.length()) / 2) - 1;
	int s1Matches[s1.length()];
	int s2Matches[s2.length()];
	memset(s1Matches, 0, s1.length() * sizeof(int));
	memset(s2Matches, 0, s2.length() * sizeof(int));

	for (int i = 0; i < s1.length(); i++)
	{

		// Low Value;
		if (i >= range)
		{
			low = i - range;
		}
		else
		{
			low = 0;
		}

		// High Value;
		if (i + range <= (s2.length() - 1))
		{
			high = i + range;
		}
		else
		{
			high = s2.length() - 1;
		}

		for (int j = low; j <= high; j++)
		{
			if (s1Matches[i] != 1 && s2Matches[j] != 1 && s1[i] == s2[j])
			{
				m += 1;
				s1Matches[i] = 1;
				s2Matches[j] = 1;
				break;
			}
		}
	}

	// Exit early if no matches were found
	if (m == 0)
	{
		return 0;
	}

	// Count the transpositions.
	for (int i = 0; i < s1.length(); i++)
	{
		if (s1Matches[i] == 1)
		{
			int j;
			for (j = k; j < s2.length(); j++)
			{
				if (s2Matches[j] == 1)
				{
					k = j + 1;
					break;
				}
			}

			if (s1[i] != s2[j])
			{
				numTrans += 1;
			}
		}
	}

	// cout << "Quad 2 Algo" << endl;
	// cout << "Match: " << m << endl; 
	// cout << "Transposition: " << (numTrans / 2) << endl;
	// cout<< endl;

	double weight = (m / s1.length() + m / s2.length() + (m - (numTrans / 2)) / m) / 3;
	// double l = 0;
	// double p = 0.1;
	// if (weight > 0.7)
	// {
	// 	while (s1[l] == s2[l] && l < 4)
	// 	{
	// 		l += 1;
	// 	}

	// 	weight += l * p * (1 - weight);
	// }
	return weight;
}

void getData(string& file_path) {
    string line;
    ifstream records(file_path);
	int largest = 0;
	int smallest = INT32_MAX;
	long long int totalChar = 0;

    int ind = 0;
    while (getline (records, line)) {
		vector<string> result;
        // boost::split(result, line, boost::is_any_of(","));
		result.push_back(line);
		recordVector.push_back(result[0]);
		if(result[0].size() > largest) {
			largest = result[0].size();
		} else if (result[0].size() < smallest) {
			smallest = result[0].size();
		}
		totalChar += result[0].size();
    }
    records.close();
	cout<< "Largest: " << largest << " Smallest: " << smallest << " AVG: " << ((double)totalChar)/((double)recordVector.size()) << endl;
	stringMaxLen = largest + 1;
}

void preprocessData() {
	for(int i=0; i< recordVector.size(); i++){
		for(int j=0; j<recordVector[i].size(); j++) {
			if(recordVector[i][j]>=97) {
				recordIntArray[i * stringMaxLen + j] = recordVector[i][j] - 97;
			} else if ((recordVector[i][j] >=48) && (recordVector[i][j] < 58)) {
				recordIntArray[i * stringMaxLen + j] = recordVector[i][j] - 48 + 26;
			}
		}
		strLength[i] = recordVector[i].size();
	}
}

int main(int argc, char **argv)
{
	string filePath = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/";
	// string filePath = "/home/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/";
	string fileName = argv[1];
    filePath = filePath + argv[1];

	getData(filePath);

	double jaroDist_2, jaroDistLinear;

	clock_t currTS_q1	= clock();
	for(int i=0; i<recordVector.size()-1; i++){
		for(int j=i+1; j<recordVector.size(); j++) {
			jaroDist_2 = jaroDistance_2(recordVector[i], recordVector[j]);
			// cout<< jaroDist_2 << endl;
		}
	}
	double quadtotalCompTime	= (double)(clock() - currTS_q1) / CLOCKS_PER_SEC;
    cout<< "QuadAlg Total Time " << quadtotalCompTime << endl;


	clock_t currTS_p1	= clock();

	jwMatS1 = new int[alphabetSize * stringMaxLen];
	jwMatS2 = new int[alphabetSize * stringMaxLen];
	s1CharCount = new int[alphabetSize];
	s2CharCount = new int[alphabetSize];
	hash_s1 = new int[stringMaxLen];
	hash_s2 = new int[stringMaxLen];
	memset(s1CharCount, 0, alphabetSize * sizeof(int));
	memset(s2CharCount, 0, alphabetSize * sizeof(int));
	memset(hash_s1, 0, stringMaxLen * sizeof(int));
	memset(hash_s2, 0, stringMaxLen * sizeof(int));
	recordIntArray = new int[recordVector.size()*stringMaxLen];
	strLength = new int[recordVector.size()];

	preprocessData();

	double preprocessingTime	= (double)(clock() - currTS_p1) / CLOCKS_PER_SEC;
	clock_t currTS_p2	= clock();
	for(int i=0; i<recordVector.size()-1; i++){
		for(int j=i+1; j<recordVector.size(); j++) {
			if(recordVector[i].size() < recordVector[j].size()) {
				jaroDistLinear = jaroDistanceLinear(i,j);
			} else {
				jaroDistLinear = jaroDistanceLinear(j,i);
			}
			// cout<< jaroDistLinear << endl;
		}
	}

	double compTime	= (double)(clock() - currTS_p2) / CLOCKS_PER_SEC;
	double totalCompTime	= (double)(clock() - currTS_p1) / CLOCKS_PER_SEC;
    cout<< "Processing time "<< preprocessingTime << " CompTime " << compTime << " Total Time " << totalCompTime << endl;

	delete[] jwMatS1;
	delete[] jwMatS2;
	delete[] s1CharCount;
	delete[] s2CharCount;
	delete[] hash_s1;
	delete[] hash_s2;
	delete[] strLength;
	delete[] recordIntArray;

	return 0;
}
