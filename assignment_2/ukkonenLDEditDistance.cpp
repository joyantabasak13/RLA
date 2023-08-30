// Implement ukkonen's aproximate string mapping algorithm
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
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
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

using namespace std;

// Jaro distance variables
const int stringMaxLen = 20;
const int alphabetSize = 26;
auto jwMatS1 = std::make_unique<int[]>(alphabetSize * stringMaxLen);
auto jwMatS2 = std::make_unique<int[]>(alphabetSize * stringMaxLen);
int s1CharCount[alphabetSize];
int s2CharCount[alphabetSize];


// Edit distance variables
int threshold = 99;
int matArr[50][50] = {0};

// let us assume all the characters are alphabets
double jaroDistImprovedLinear(string& s1, string& s2){
	int len1 = s1.length();
	int len2 = s2.length();

	if (len1 == 0 || len2 == 0)
		return 0.0;

	int range = floor(max(len1, len2) / 2) - 1;

	int match = 0;
	int ptr1, ptr2;
	int count1, count2;
	int hash_s1[len1];
	int hash_s2[len2];

	memset(hash_s1, 0, len1 * sizeof(int));
	memset(hash_s2, 0, len2 * sizeof(int));
	memset(s1CharCount, 0, alphabetSize * sizeof(int));
	memset(s2CharCount, 0, alphabetSize * sizeof(int));

	for(int i=0; i<len1; i++) {
		jwMatS1[(s1[i]-97) * stringMaxLen + s1CharCount[s1[i]-97]] = i;
		s1CharCount[s1[i]-97]++;
	}

	for(int i=0; i<len2; i++) {
		jwMatS2[(s2[i]-97) * stringMaxLen + s2CharCount[s2[i]-97]] = i;
		s2CharCount[s2[i]-97]++;
	}

	for(int i=0; i<len1; i++) {
		count1 = s1CharCount[s1[i]-97];
		count2 = s2CharCount[s1[i]-97];
		if((count1 != 0) && (count2 != 0)) {
			ptr1 = 0;
			ptr2 = 0;
			while((ptr1<count1) && (ptr2<count2))
			{
				if (abs(jwMatS1[(s1[i]-97) * stringMaxLen + ptr1] - jwMatS2[(s1[i]-97) * stringMaxLen + ptr2]) <= range)
				{
					match++;
					cout<< "Matched Ch: " << s1[i] << " Ind1: " << jwMatS1[(s1[i]-97) * stringMaxLen + ptr1] << " Ind: " << jwMatS2[(s1[i]-97) * stringMaxLen + ptr2] << endl;
					hash_s1[jwMatS1[(s1[i]-97) * stringMaxLen + ptr1]] = 1;
					hash_s2[jwMatS2[(s1[i]-97) * stringMaxLen + ptr2]] = 1;
					ptr1++;
					ptr2++;
				} else {
					cout<< "Not Matched Ch: " << s1[i] << " Ind1: " << jwMatS1[(s1[i]-97) * stringMaxLen + ptr1] << " Ind: " << jwMatS2[(s1[i]-97) * stringMaxLen + ptr2] << endl;
					if (jwMatS2[(s1[i]-97) * stringMaxLen + ptr2] < jwMatS1[(s1[i]-97) * stringMaxLen + ptr1])
					{
						ptr2++;
					} else
					{
						ptr1++;
					}
				}
			}
			s1CharCount[s1[i]-97] = 0;
			s2CharCount[s1[i]-97] = 0;
		}
	}

	int t = 0;
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
					break;
				}
			}
			if (s1[i] != s2[j])
			{
				t++;
			}
		}
	}

	t = t / 2;

	cout << "Improved linear jaro" << endl;
	cout << "Range: " << range << endl;
	cout << "Match: " << match << endl;
	cout << "Total Length: " << len1 + len2 << endl;
	cout << "Transposition: " << t << endl;

	return (((double)match) / ((double)len1) + ((double)match) / ((double)len2) + ((double)match - t) / ((double)match)) / 3.0;


}



// helps edit distance calculation in calculateBasicED()
int calculateBasicED2(string &str1, string &str2, int threshRem)
{
	int row, col, i, j;
	row = str1.length() + 1;
	col = str2.length() + 1;

	for (i = 0; i < row; i++)
	{
		for (j = 0; j < col; j++)
		{
			if (i == 0)
				matArr[i][j] = j;
			else if (j == 0)
				matArr[i][j] = i;
			else
			{
				if ((int)str1[i - 1] == (int)str2[j - 1])
					matArr[i][j] = min(min(matArr[i - 1][j - 1], matArr[i - 1][j] + 1), matArr[i][j - 1] + 1);
				else
					matArr[i][j] = min(min(matArr[i - 1][j - 1] + 1, matArr[i - 1][j] + 1), matArr[i][j - 1] + 1);
				if (i > 1 && j > 1)
				{
					if ((str1[i - 1] == str2[j - 2]) && (str1[i - 2] == str2[j - 1]))
					{
						matArr[i][j] = min(matArr[i][j], matArr[i - 2][j - 2] + 1);
					}
				}
			}

			if ((row - col) == (i - j) && (matArr[i][j] > threshRem))
			{
				return threshold + 1;
			}
		}
	}

	return (matArr[row - 1][col - 1]);
}

// calculates edit distance between two string
// takes two strings and a threshold value as input
// returns global variable threshold + 1 if distance exceeds theshold param
// returns edit distance
// core mechanism is a DP algo
int calculateBasicED(string &str1, string &str2, int threshRem)
{
	int dist = threshRem;

	if (abs((int)(str1.length() - str2.length())) > dist)
		return threshold + 1;
	else if ((2 * dist + 1) >= max(str1.length(), str2.length()))
		return calculateBasicED2(str1, str2, dist);
	else
	{

		string s1, s2;
		int row, col, diagonal;
		int i, j;

		if (str1.length() > str2.length())
		{
			s1 = str2;
			s2 = str1;
		}
		else
		{
			s1 = str1;
			s2 = str2;
		}

		row = s1.length() + 1;
		col = 2 * dist + 1;
		diagonal = dist + s2.length() - s1.length();

		for (i = 0; i < dist + 1; i++)
		{
			for (j = dist - i; j < col; j++)
			{
				if (i == 0)
					matArr[i][j] = j - dist;
				else if (j == (dist - i))
					matArr[i][j] = matArr[i - 1][j + 1] + 1;
				else if (j != (col - 1))
				{
					if ((int)s1[i - 1] == (int)s2[j - (dist - i) - 1])
						matArr[i][j] = min(min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
					else
					{
						matArr[i][j] = min(min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
					}
				}
				else
				{
					if ((int)s1[i - 1] == (int)s2[j - (dist - i) - 1])
						matArr[i][j] = min(matArr[i - 1][j], matArr[i][j - 1] + 1);
					else
						matArr[i][j] = min(matArr[i - 1][j] + 1, matArr[i][j - 1] + 1);
				}

				if ((j == diagonal) && matArr[i][j] > dist)
					return threshold + 1;
			}
		}

		for (i = dist + 1; i < s2.length() - dist + 1; i++)
		{
			for (j = 0; j < col; j++)
			{
				if (j == 0)
				{
					if ((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j] = min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1);
					else
						matArr[i][j] = min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1);
				}
				else if (j != (col - 1))
				{
					if ((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j] = min(min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
					else
						matArr[i][j] = min(min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
				}
				else
				{
					if ((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j] = min(matArr[i - 1][j], matArr[i][j - 1] + 1);
					else
						matArr[i][j] = min(matArr[i - 1][j] + 1, matArr[i][j - 1] + 1);
				}
				if ((j == diagonal) && (matArr[i][j] > dist))
					return threshold + 1;
			}
		}

		for (i = s2.length() - dist + 1; i < row; i++)
		{
			for (j = 0; j < col - i + s2.length() - dist; j++)
			{
				if (j == 0)
				{
					if ((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j] = min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1);
					else
						matArr[i][j] = min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1);
				}
				else
				{
					if ((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j] = min(min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
					else
						matArr[i][j] = min(min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
				}
				if ((j == diagonal) && (matArr[i][j] > dist))
					return threshold + 1;
			}
		}
		return matArr[row - 1][diagonal];
	}
}

ptrdiff_t algo_8(const std::string &A, ptrdiff_t m, const std::string &B, ptrdiff_t n, ptrdiff_t inf, ptrdiff_t k, ptrdiff_t p, size_t num_p, ptrdiff_t *f)
{
	ptrdiff_t t = f[num_p * (k + m) + (p)] + 1; // f(k,p-1) + 1
	ptrdiff_t t1 = -inf;
	if (k + m >= 1 && p >= 0)
	{
		t1 = f[num_p * (k + m - 1) + (p)];
	}

	ptrdiff_t t2 = -inf;
	if (k < n && p < inf + 2)
	{
		t2 = f[num_p * (k + m + 1) + (p)] + 1; // f(k+1,p-1) + 1
	}

	// check whether a(t)a(t+1) = b(k+t+1)b(k+t), but *only* if those
	// ranges are valid. IOW, if:
	//
	//   - t > 0
	//   - t + 1 <= m
	//   - k + t > 0
	//   - k + t + 1 <= n
	ptrdiff_t t3 = -inf;
	if (t > 0 &&
		t + 1 <= m &&
		k + t > 0 &&
		k + t + 1 <= n)
	{
		if (A[t - 1] == B[k + t] && A[t] == B[k + t - 1])
		{
			t3 = t + 1;
		}
	}
	// t := max(t,t1,t2,t3)
	if (t < t1)
	{
		t = t1;
	}
	if (t < t2)
	{
		t = t2;
	}
	if (t < t3)
	{
		t = t3;
	}
	// while a(t+1) = b(t+1+k) do t := t + 1
	while (t + 1 <= m &&
		   t + k + 1 <= n &&
		   A[t] == B[t + k])
		t++;

	if (t > m || t + k > n)
	{
		t = inf;
	}
	return t;
}

ptrdiff_t ukkonen(const std::string &A, const std::string &B, std::size_t D)
{
	ptrdiff_t m = A.length();
	ptrdiff_t n = B.length();
	ptrdiff_t inf = max(m, n); // |A,B| <= inf

	// Allocating too much space, here. This implementation won't satisfy the
	// space bounds.

	// To index into `f' in terms of (i,j), -m <= i <= n,
	// -1 <=j <= inf, do f[i+m][j+1]
	ptrdiff_t f[(m + n + 1) * (inf + 2)];
	for (ptrdiff_t i = 0; i < m + n + 1; ++i)
	{
		for (ptrdiff_t j = 0; j < inf + 2; ++j)
		{
			f[i * (inf + 2) + j] = -inf - 1;
		}
	}

	// Initialize f: f(k,|k|-1) = |k|-1, if k < 0...
	for (ptrdiff_t k = -1; k >= -m; --k)
	{
		f[(k + m) * (inf + 2) - k] = -k - 1;
	}
	for (ptrdiff_t k = 0; k <= n; ++k)
	{
		f[(k + m) * (inf + 2) + k] = -1;
	}

	ptrdiff_t p = -1;
	ptrdiff_t r = p - min(m, n);
	while (f[n * (inf + 2) + p + 1] != m)
	{
		p = p + 1;
		r = r + 1;
		if (r <= 0)
		{
			for (ptrdiff_t k = -p; k <= p; ++k)
			{
				// f(k,p)
				f[(inf + 2) * (k + m) + p + 1] = algo_8(A, m, B, n, inf, k, p, inf + 2, (ptrdiff_t *)f);
			}
		}
		else
		{
			for (ptrdiff_t k = max(-m, -p); k <= -r; ++k)
			{
				// f(k,p)
				f[(inf + 2) * (k + m) * p + 1] = algo_8(A, m, B, n, inf, k, p, inf + 2, (ptrdiff_t *)f);
			}
			for (ptrdiff_t k = r; k <= max(n, p); ++k)
			{
				// f(k,p)
				f[(inf + 2) * (k + m) + p + 1] = algo_8(A, m, B, n, inf, k, p, inf + 2, (ptrdiff_t *)f);
			}
		}
		if (p > D)
		{
			return threshold + 1;
		}
	}
	return p;
}

void sortCharIDArray(pair<int, int> charIDList[], int length)
{
	int alphabets = 36;
	int numRecords = length;
	pair<int, int> tempArr[numRecords];
	int countArr[alphabets];
	memset(countArr, 0, alphabets * sizeof(int));

	for (int j = 0; j < numRecords; ++j)
	{
		countArr[charIDList[j].first]++;
	}
	// Do prefix sum
	for (int k = 1; k < alphabets; ++k)
		countArr[k] += countArr[k - 1];

	for (int j = numRecords - 1; j >= 0; --j)
		tempArr[--countArr[charIDList[j].first]] = charIDList[j];

	for (int j = 0; j < numRecords; ++j)
		charIDList[j] = tempArr[j];
}

double jaro_distance(string &s1, string &s2)
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

	cout << "Range: " << max_dist << endl;
	cout << "Match: " << match << endl;
	cout << "Total Length: " << len1 + len2 << endl;
	cout << "Transposition: " << t << endl;

	// Return the Jaro Similarity
	return (((double)match) / ((double)len1) + ((double)match) / ((double)len2) + ((double)match - t) / ((double)match)) / 3.0;
}

double jaro_dist_2(std::string &s1, std::string &s2)
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

	cout << "Range: " << range << endl;
	cout << "Match: " << m << endl;
	cout << "Total Length: " << s1.length() + s2.length() << endl; 
	cout << "Transposition: " << (numTrans / 2) << endl;

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

double jaro_distance_linear(string &s1, string &s2)
{
	cout << endl;
	cout<< s1 << " " << s2<< endl;
	int len1 = s1.length();
	int len2 = s2.length();

	if (len1 == 0 || len2 == 0)
		return 0.0;

	int range = floor(max(len1, len2) / 2) - 1;

	int match = 0;
	int misMatch = 0;
	pair<int, int> s1charIndPair[len1];
	pair<int, int> s2charIndPair[len2];
	int hash_s1[len1];
	int hash_s2[len2];
	memset(hash_s1, 0, len1 * sizeof(int));
	memset(hash_s2, 0, len2 * sizeof(int));

	pair<int, int> p;

	for (int i = 0; i < len1; i++)
	{
		if ((s1[i] >= 97) && (s1[i] <= 122))
		{
			p.first = s1[i] - 97;
		}
		else if ((s1[i] >= 48) && (s1[i] <= 57))
		{
			p.first = s1[i] - 48 + 26;
		}
		p.second = i;
		s1charIndPair[i] = p;
	}

	for (int i = 0; i < len2; i++)
	{
		if ((s2[i] >= 97) && (s2[i] <= 122))
		{
			p.first = s2[i] - 97;
		}
		else if ((s2[i] >= 48) && (s2[i] <= 57))
		{
			p.first = s2[i] - 48 + 26;
		}
		p.second = i;
		s2charIndPair[i] = p;
		// cout<< p.first << endl;
	}
	sortCharIDArray(s1charIndPair, len1);
	sortCharIDArray(s2charIndPair, len2);
	int j = 0;
	int i = 0;
	while ((i < len1) && (j < len2))
	{
		if (s1charIndPair[i].first == s2charIndPair[j].first)
		{
			if (abs(s1charIndPair[i].second - s2charIndPair[j].second) <= range)
			{
				match++;
				hash_s1[s1charIndPair[i].second] = 1;
				hash_s2[s2charIndPair[j].second] = 1;
				j++;
				i++;
			}
			else
			{
				misMatch++;
				if (s1charIndPair[i].second < s2charIndPair[j].second)
				{
					i++;
				}
				else
				{
					j++;
				}
			}
		}
		else
		{
			misMatch++;
			if (s1charIndPair[i].first < s2charIndPair[j].first)
			{
				i++;
			}
			else
			{
				j++;
			}
		}
	}

	misMatch += len1 + len2 - i - j;
	if (match == 0)
		return 0.0;

	int t = 0;
	int point = 0;
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
					break;
				}
			}
			if (s1[i] != s2[j])
			{
				t++;
			}
		}
	}
	cout<< "transpositions: " << t << endl; 
	t = t / 2;

	cout<< "Non thresholded linear jaro" << endl;
	cout << "Range: " << range << endl;
	cout << "Match: " << match << endl;
	cout << "Mismatch: " << misMatch << endl;
	cout << "Total Length: " << len1 + len2 << endl;
	cout << "Transposition: " << t << endl;

	return (((double)match) / ((double)len1) + ((double)match) / ((double)len2) + ((double)match - t) / ((double)match)) / 3.0;
}

double jaro_distance_thresholded(string &s1, string &s2, double threshold)
{
	cout << endl;
	cout<< s1 << " " << s2 << endl;
	int len1 = s1.length();
	int len2 = s2.length();
	double len_mul = (double)(len1 * len2);
	int misMatchAllowed = ceil((len_mul * (4.0 - 6.0 * threshold) + (double)pow(len1, 2) + (double)pow(len2, 2)) / (len1 + len2));
	if (len1 == 0 || len2 == 0)
		return 0.0;

	int range_ind = floor(max(len1, len2) / 2) - 1;

	int match = 0;
	int misMatch = 0;
	vector<pair<int, int>> s1charIndPair;
	vector<pair<int, int>> s2charIndPair;
	s1charIndPair.resize(len1);
	s2charIndPair.resize(len2);
	int hash_s1[len1];
	int hash_s2[len2];
	memset(hash_s1, 0, len1 * sizeof(int));
	memset(hash_s2, 0, len2 * sizeof(int));
	pair<int, int> p;

	for (int i = 0; i < len1; i++)
	{
		if ((s1[i] >= 97) && (s1[i] <= 122))
		{
			p.first = s1[i] - 97;
		}
		else if ((s1[i] >= 48) && (s1[i] <= 57))
		{
			p.first = s1[i] - 48 + 26;
		}
		p.second = i;
		s1charIndPair[i] = p;
	}

	for (int i = 0; i < len2; i++)
	{
		if ((s2[i] >= 97) && (s2[i] <= 122))
		{
			p.first = s2[i] - 97;
		}
		else if ((s2[i] >= 48) && (s2[i] <= 57))
		{
			p.first = s2[i] - 48 + 26;
		}
		p.second = i;
		s2charIndPair[i] = p;
		// cout<< p.first << endl;
	}

	stable_sort(s1charIndPair.begin(), s1charIndPair.end(), [](const auto& a, const auto& b) {return a.first < b.first; });
	stable_sort(s2charIndPair.begin(), s2charIndPair.end(), [](const auto& a, const auto& b) {return a.first < b.first; });

	int j = 0;
	int i = 0;
	while ((i < len1) && (j < len2))
	{
		if (s1charIndPair[i].first == s2charIndPair[j].first)
		{
			if (abs(s1charIndPair[i].second - s2charIndPair[j].second) <= range_ind)
			{
				match++;
				hash_s1[s1charIndPair[i].second] = 1;
				hash_s2[s2charIndPair[j].second] = 1;
				j++;
				i++;
			}
			else
			{
				misMatch++;
				if (s1charIndPair[i].second < s2charIndPair[j].second)
				{
					i++;
				}
				else
				{
					j++;
				}
			}
		}
		else
		{
			misMatch++;
			if (s1charIndPair[i].first < s2charIndPair[j].first)
			{
				i++;
			}
			else
			{
				j++;
			}
		}
		if (misMatch > misMatchAllowed)
			return 0.0;
	}

	misMatch += len1 + len2 - i - j;
	if (misMatch > misMatchAllowed)
		return 0.0;

	int t = 0;
	int point = 0;
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
					break;
				}
			}
			if (s1[i] != s2[j])
			{
				t++;
			}
		}
	}
	t = t / 2;

	// cout<< "Thresholded Jaro" << endl;
	// cout<< "Range: " << range_ind << endl;
	// cout<< "Match: " << match << endl;
	// cout<< "Mismatch: " << misMatch << endl;
	// cout<< "Total Length: " << len1 + len2 << endl;
	// cout<< "Transposition: " << t << endl;

	return (((double)match) / ((double)len1) + ((double)match) / ((double)len2) + ((double)match - t) / ((double)match)) / 3.0;
}

double jaro_distance_linear2(string &s1, string &s2, double th)
{
	int len1 = s1.length();
	int len2 = s2.length();

	if (len1 == 0 || len2 == 0)
		return 0.0;

	int alphabets = 36;
	int range = floor(max(len1, len2) / 2) - 1;
	double len_mul = (double)(len1 * len2);
	int misMatchAllowed = ceil((len_mul * (4.0 - 6.0 * th) + (double)pow(len1, 2) + (double)pow(len2, 2)) / (len1 + len2));
	
	int match = 0;
	int misMatch = 0;
	pair<int, int> s1charIndPair[len1];
	pair<int, int> s2charIndPair[len2];
	int hash_s1[len1];
	int hash_s2[len2];
	vector<int> c1[alphabets];
	vector<int> c2[alphabets];
	memset(hash_s1, 0, len1 * sizeof(int));
	memset(hash_s2, 0, len2 * sizeof(int));

	int cIndex;
	int charInStrIndex;
	for (int i = 0; i < len1; i++)
	{
		if ((s1[i] >= 97) && (s1[i] <= 122))
		{
			cIndex = s1[i] - 97;
		}
		else if ((s1[i] >= 48) && (s1[i] <= 57))
		{
			cIndex = s1[i] - 48 + 26;
		}
		charInStrIndex = i;
		c1[cIndex].push_back(charInStrIndex);
	}

	for (int i = 0; i < len2; i++)
	{
		if ((s2[i] >= 97) && (s2[i] <= 122))
		{
			cIndex = s2[i] - 97;
		}
		else if ((s2[i] >= 48) && (s2[i] <= 57))
		{
			cIndex = s2[i] - 48 + 26;
		}
		charInStrIndex = i;
		c2[cIndex].push_back(charInStrIndex);
	}

	int j = 0;
	int i = 0;
	for(int k=0;k<alphabets;k++){
		while((i<c1[k].size()) && (j<c2[k].size())){
			if(abs(c1[k][i]-c2[k][j]) <= range) {
				match++;
				hash_s1[c1[k][i]] = 1;
				hash_s2[c2[k][j]] = 1;
				i++;
				j++;
			} else {
				misMatch++;
				if(c1[k][i] < c2[k][j]){
					i++;
				} else {
					j++;
				}
			}
		}
		misMatch += c1[k].size() + c2[k].size() - i -j;

		if(misMatch>misMatchAllowed) return 0.0;

		i=0;
		j=0;
	}

	int t = 0;
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
					break;
				}
			}
			if (s1[i] != s2[j])
			{
				t++;
			}
		}
	}
	t = t / 2;

	cout<< "Non thresholded linear jaro" << endl;
	cout << "Range: " << range << endl;
	cout << "Match: " << match << endl;
	cout << "Mismatch: " << misMatch << endl;
	cout << "Total Length: " << len1 + len2 << endl;
	cout << "Transposition: " << t << endl;

	return (((double)match) / ((double)len1) + ((double)match) / ((double)len2) + ((double)match - t) / ((double)match)) / 3.0;
}



double calculateJaroWinklerDist(string &str1, string &str2)
{
	int min_string_len = (str1.size() - str2.size());

	if (min_string_len < 0)
	{
		min_string_len = -min_string_len;
	}
	if (min_string_len > 1)
	{
		return 0.20;
	}

	double jaro_dist = jaro_distance(str1, str2);

	if (jaro_dist > 0.0)
	{
		int prefix = 0;

		for (int i = 0; i < min(str1.length(), str2.length()); i++)
		{
			// If the characters match
			if (str1[i] == str2[i])
				prefix++;

			// Else break
			else
				break;
		}

		// Maximum of 4 characters are allowed in prefix
		prefix = min(4, prefix);

		// Calculate jaro winkler Similarity
		jaro_dist += 0.1 * prefix * (1 - jaro_dist);
	}
	return jaro_dist;
}

double jaroDistanceLinear(string& s1, string& s2){
	// cout<< s1 << endl;
	int len1 = s1.length();
	int len2 = s2.length();

	if (len1 == 0 || len2 == 0)
		return 0.0;

	int range = floor(max(len1, len2) / 2) - 1;

	int match = 0;
	int ptr1, ptr2;
	int count1, count2;

	// memset(hash_s1, 0, len1 * sizeof(int));
	// memset(hash_s2, 0, len2 * sizeof(int));
	// memset(s1CharCount, 0, alphabetSize * sizeof(int));
	// memset(s2CharCount, 0, alphabetSize * sizeof(int));
	// cout<< "error" << endl;
	for(int i=0; i<len1; i++) {
		jwMatS1[(s1[i]-97) * stringMaxLen + s1CharCount[s1[i]-97]] = i;
		s1CharCount[s1[i]-97]+=1;
	}
	// cout<< "more error" << endl;
	for(int i=0; i<len2; i++) {
		jwMatS2[(s2[i]-97) * stringMaxLen + s2CharCount[s2[i]-97]] = i;
		s2CharCount[s2[i]-97]++;
	}
	// cout<< "no error" << endl;
	// cout<< "Error!" << endl; 
	for(int i=0; i<len1; i++) {
		if(((s1[i]-97) >= alphabetSize) || (s1[i]-97 < 0)) {
			cout<< "Error!" << endl; 
		}
		count1 = s1CharCount[s1[i]-97];
		count2 = s2CharCount[s1[i]-97];
		if((count1 != 0) && (count2 != 0)) {
			ptr1 = 0;
			ptr2 = 0;
			while((ptr1<count1) && (ptr2<count2))
			{
				if (abs(jwMatS1[(s1[i]-97) * stringMaxLen + ptr1] - jwMatS2[(s1[i]-97) * stringMaxLen + ptr2]) <= range)
				{
					match++;
					// cout<< "Matched Ch: " << s1[i] << " Ind1: " << jwMatS1[(s1[i]-97) * stringMaxLen + ptr1] << " Ind: " << jwMatS2[(s1[i]-97) * stringMaxLen + ptr2] << endl;
					hash_s1[jwMatS1[(s1[i]-97) * stringMaxLen + ptr1]] = 1;
					hash_s2[jwMatS2[(s1[i]-97) * stringMaxLen + ptr2]] = 1;
					ptr1++;
					ptr2++;
				} else {
					// cout<< "Not Matched Ch: " << s1[i] << " Ind1: " << jwMatS1[(s1[i]-97) * stringMaxLen + ptr1] << " Ind: " << jwMatS2[(s1[i]-97) * stringMaxLen + ptr2] << endl;
					if (jwMatS2[(s1[i]-97) * stringMaxLen + ptr2] < jwMatS1[(s1[i]-97) * stringMaxLen + ptr1])
					{
						ptr2++;
					} else
					{
						ptr1++;
					}
				}
			}
		}
		s1CharCount[s1[i]-97] = 0;
	}
	for(int i=0; i<len2; i++) {
		s2CharCount[s2[i]-97] = 0;
	}
	int t = 0;
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
			if (s1[i] != s2[j])
			{
				t++;
			}
		}
		hash_s1[i] = 0;
	}

	t = t / 2;

	// cout << "Improved linear jaro" << endl;
	// cout << "Match: " << match << endl;
	// cout << "Transposition: " << t << endl;
	// cout << endl;
	// for(int i=0; i<alphabetSize; i++) {
	// 	if(s1CharCount[i]!=0) {
	// 		cout<< "MASSIVE PROBLEM" << endl;
	// 		exit(0);
	// 	}
	// }

	return (((double)match) / ((double)len1) + ((double)match) / ((double)len2) + ((double)match - t) / ((double)match)) / 3.0;


}


int main(int argc, char **argv)
{
	string A = "marhta";
	string B = "martha";
	// string A = "joyanta";
	// string B = "joyanta";

	double jaro_dist, jaro_dist_linear, jaro_dist_arbitrater, jaro_dist_thresholded;
	double jaro_dist_linear_2, jaro_dist_improvedLinear;

	jaro_dist_improvedLinear = jaroDistImprovedLinear(A,B);
	cout<< endl;
	// jaro_dist = jaro_distance(A, B);
	// jaro_dist_linear = jaro_distance_linear(A, B);
	jaro_dist_arbitrater = jaro_dist_2(A, B);
	// jaro_dist_thresholded = jaro_distance_thresholded(A, B, 0.5);
	// jaro_dist_linear_2 = jaro_distance_linear2(A,B, 0.5);

	// cout << "Jaro Similarity: " << jaro_dist << endl;
	// cout << "linear Jaro Similarity: " << jaro_dist_linear << endl;
	cout << "Arbitrater Jaro Similarity: " << jaro_dist_arbitrater << endl;
	// cout << "linear thresholded Jaro Similarity: " << jaro_dist_thresholded << endl;
	// cout << "Jaro distance Linear 2: " <<  jaro_dist_linear_2 << endl;
	cout << "Jaro Dist Improved Linear: " << jaro_dist_improvedLinear << endl;
	return 0;


	A = "leland";
	B = "lindaa";
	// jaro_dist = jaro_distance(A, B);
	// jaro_dist_linear = jaro_distance_linear(A, B);
	// jaro_dist_arbitrater = jaro_dist_2(A, B);
	jaro_dist_thresholded = jaro_distance_thresholded(A, B, 0.8);
	jaro_dist_linear_2 = jaro_distance_linear2(A,B, 0.8);
	

	// cout << "Jaro Similarity: " << jaro_dist << endl;
	// cout << "linear Jaro Similarity: " << jaro_dist_linear << endl;
	// cout << "Arbitrater Jaro Similarity: " << jaro_dist_arbitrater << endl;
	cout << "linear thresholded Jaro Similarity: " << jaro_dist_thresholded << endl;
	cout << "Jaro distance Linear 2: " <<  jaro_dist_linear_2 << endl;
	return 0;
	A = "king";
	B = "evals";
	jaro_dist = jaro_distance(A, B);
	jaro_dist_linear = jaro_distance_linear(A, B);
	jaro_dist_arbitrater = jaro_dist_2(A, B);
	jaro_dist_thresholded = jaro_distance_thresholded(A, B, 0.8);
	jaro_dist_linear_2 = jaro_distance_linear2(A,B, 0.8);

	cout << "Jaro Similarity: " << jaro_dist << endl;
	cout << "linear Jaro Similarity: " << jaro_dist_linear << endl;
	cout << "Arbitrater Jaro Similarity: " << jaro_dist_arbitrater << endl;
	cout << "linear thresholded Jaro Similarity: " << jaro_dist_thresholded << endl;
	cout << "Jaro distance Linear 2: " <<  jaro_dist_linear_2 << endl;

	return 0;
}
