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
#include <tuple>
#include <time.h>
#include <unistd.h>
#include <utility>
#include <vector>
#include <mutex> 

#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

using namespace std;

int threshold = 99;
int totalRecords;
int lenMax;
int attributes;
int blockFieldIndex = 1;

vector<string> vec1D;
vector<vector<string> > vec2D;
vector<vector<vector<string> > > sameEntities;
vector<pair<int, string> > combinedData;
vector<vector<vector<int> > > distances;

// helps edit distance calculation in calculateBasicED()
int calculateBasicED2(string& str1, string& str2, int threshRem)
{
	int row, col, i, j;
	vector<vector<int> > matArr;

	row		= str1.length() + 1;
	col 	= str2.length() + 1;

	matArr.resize(row);
	for(i = 0; i < row; ++i)
		matArr[i].resize(col, 0);

	for(i = 0; i < row; i++)
	{
		for(j = 0; j < col; j++)
		{
			if (i == 0)
				matArr[i][j] = j;
			else if (j == 0)
				matArr[i][j] = i;
			else
			{
				if((int)str1[i-1] == (int)str2[j-1])
					matArr[i][j]	= min(min(matArr[i - 1][j - 1], matArr[i - 1][j] + 1), matArr[i][j - 1] + 1);
				else
					matArr[i][j] 	= min(min(matArr[i - 1][j - 1] + 1, matArr[i - 1][j] + 1), matArr[i][j - 1] + 1);
			}

			if((row - col) == (i - j) && (matArr[i][j] > threshRem))
			{
				return threshold + 1;
			}
		}
	}

	return (matArr[row-1][col-1]);
}

// calculates edit distance between two string
// takes two strings and a threshold value as input
// returns global variable threshold + 1 if distance exceeds theshold param
// returns edit distance
// core mechanism is a DP algo 
int calculateBasicED(string& str1, string& str2, int threshRem)
{
	int dist	= threshRem;

	if(abs((int)(str1.length() - str2.length())) > dist)
		return threshold + 1;
	else if(str1.compare(str2) == 0)
		return 0;
	else if(dist == 0)
		return threshold + 1;
	else if((2 * dist + 1) >= max(str1.length(), str2.length()))
		return calculateBasicED2(str1, str2, dist);
	else
	{
		string s1, s2;
		int row, col, diagonal;
		int i, j;
		vector<vector<int> > matArr;

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

		row	 		= s1.length() + 1;
		col 		= 2 * dist + 1;
		diagonal 	= dist + s2.length() - s1.length();

		matArr.resize(row);
		for(i = 0; i < row; ++i)
			matArr[i].resize(col, 0);

		//if(procID == 1 && checkTemp == 3164)
			//	cout << str1 << " -- " << str2 << " rt " << dist << endl;

		for(i = 0; i < dist + 1; i++)
		{
			for(j = dist - i; j < col; j++)
			{
				if (i == 0)
					matArr[i][j]	= j - dist;
				else if(j == (dist - i))
					matArr[i][j] 	= matArr[i - 1][j + 1] + 1;
				else if(j != (col - 1))
				{
					if((int)s1[i - 1] == (int)s2[j - (dist - i) - 1])
						matArr[i][j]	= min(min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
					else
						matArr[i][j] 	= min(min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
				}
				else
				{
					if((int)s1[i - 1] == (int)s2[j - (dist - i) - 1])
						matArr[i][j]	= min(matArr[i - 1][j], matArr[i][j - 1] + 1);
					else
						matArr[i][j] 	= min(matArr[i - 1][j] + 1, matArr[i][j - 1] + 1);
				}

				if((j == diagonal) && matArr[i][j] > dist)
					return threshold + 1;
			}
		}

		for(i = dist + 1; i < s2.length() - dist + 1; i++)
		{
			for(j = 0; j < col; j++)
			{
				if(j == 0)
				{
					if((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j]	= min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1);
					else
						matArr[i][j] 	= min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1);
				}
				else if(j != (col - 1))
				{
					if((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j] 	= min(min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
					else
						matArr[i][j] 	= min(min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
				}
				else
				{
					if((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j] 	= min(matArr[i - 1][j], matArr[i][j - 1] + 1);
					else
						matArr[i][j] 	= min(matArr[i - 1][j] + 1, matArr[i][j - 1] + 1);
				}
				if((j == diagonal) && (matArr[i][j] > dist))
					return threshold + 1;
			}
		}

		for(i = s2.length() - dist + 1; i < row; i++)
		{
			for(j = 0; j < col - i + s2.length() - dist; j++)
			{
				if(j == 0)
				{
					if((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j]	= min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1);
					else
						matArr[i][j] 	= min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1);
				}
				else
				{
					if((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j] 	= min(min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
					else
						matArr[i][j] 	= min(min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
				}
				if((j == diagonal) && (matArr[i][j] > dist))
					return threshold + 1;
			}
		}

		//if(procID == 1 && checkTemp == 3164)
			//cout << str1 << " -- " << str2 << " hukjhk " << matArr[row - 1][diagonal] << endl;

		return matArr[row - 1][diagonal];

	}
}

void radixSort(vector<pair<int, string> > &strDataArr){
	int numRecords = strDataArr.size();
	vector<pair<int, string>> tempArr(numRecords);
	
	for (int i = lenMax - 1; i >= 0; --i) {
		vector<int> countArr(256, 0);
		
		for (int j = 0; j < numRecords; ++j) {
			countArr[(strDataArr[j].second)[i]]++;
		}
		
		for (int k = 1; k < 256; ++k)
			countArr[k]	+= countArr[k - 1];
		
		for (int j = numRecords - 1; j >= 0; --j)
			tempArr[--countArr[(strDataArr[j].second)[i]]]	= strDataArr[j];
		
		for (int j = 0; j < numRecords; ++j)
			strDataArr[j]	= tempArr[j];
	}
}

void printSortedRecords() {
    for (int i =0; i< combinedData.size(); i++) {
        cout<< combinedData[i].first << " " << combinedData[i].second << endl;   
    }
}

bool isLinkageOk(vector<string> &a, vector<string> &b, int blockingThreshold, int singleNonBlockingAttrThreshold, int totalNonBlockingAttrThreshold)
{
	int blockField_dist = calculateBasicED(a[1], b[1], blockingThreshold);
    if (blockField_dist <= threshold) {
		int singleAttributeDist = 0;
		int cumulativeNonBlockingDist = 0;
		for (int i = 2; i < a.size(); i++)
		{
			singleAttributeDist = calculateBasicED(a[i], b[i], singleNonBlockingAttrThreshold);
			cumulativeNonBlockingDist += singleAttributeDist;
			if (cumulativeNonBlockingDist > totalNonBlockingAttrThreshold){
				return false;
			}
		}
		return true;     
    }
	return false;
}

// I/O Functions

void getFormattedDataFromCSV(string& file_path) {
    string line;
    ifstream records(file_path);

    int ind = 0;
    while (getline (records, line)) {
		// cout<< line << endl;
        vector<string> result;
        boost::split(result, line, boost::is_any_of(","));
        vector<string> vec;
        vec.push_back(result[0]);
		string blockStr = result[blockFieldIndex];
		//Remove digits from name
		auto last = std::remove_if(blockStr.begin(), blockStr.end(), [](auto ch) {
        								return ::isdigit(ch) || ::ispunct(ch) || ::iswpunct(ch);
    								});
		blockStr.erase(last, blockStr.end()); //Remove junk left by remove_if() at the end of iterator
        boost::to_lower(blockStr);
		vec.push_back(blockStr);
		for(int i = 1; i<result.size(); i++){
			if (i!=blockFieldIndex) {
				vec.push_back(result[i]);
			}
		}
        vec2D.push_back(vec);
    }
    records.close();
    vec1D.resize(vec2D[0].size()*vec2D.size());

    for (size_t i = 0; i < vec2D.size(); i++)
    {
        for(size_t j = 0; j< vec2D[0].size(); j++) {
            vec1D[i*vec2D[0].size()+j] = vec2D[i][j];
        }
    }
    //vec2D = temp2D;
    attributes = vec2D[0].size();
    cout<< "Attributes: "<<attributes << endl;
}

void getCombinedData() {
	string strSample(50, '0');
	combinedData.resize(totalRecords);
	int max = 0;
	for (int i = 0; i< vec2D.size(); i++) {
		pair<int, string> p;
		p.first = i;
		p.second = vec2D[i][0];
		combinedData[i]=p;
		if (max<p.second.size()) {
			max = p.second.size();
		}
	}
	lenMax = max;
	// Padding to make all characters same size
    for (int i = 0; i < totalRecords; ++i) {
		int lenDiff		= lenMax - combinedData[i].second.length();
		if(lenDiff > 0)
			combinedData[i].second	+= strSample.substr(0, lenDiff);
	}
}

void getSameEntities() {
	vector<vector<string> > cluster;
	string recID = vec2D[combinedData[0].first][0];
	for (size_t i = 0; i < totalRecords; i++)
    {
		if (i == 0) {
			cluster.push_back(vec2D[combinedData[i].first]);
			continue;
		}
		string cur_recID = vec2D[combinedData[i].first][0];
        if (!recID.compare(cur_recID)) {
			// cout<< "recID: " << recID << " curID " << cur_recID << endl;
			cluster.push_back(vec2D[combinedData[i].first]);
		} else {
			recID = cur_recID;
			sameEntities.push_back(cluster);
			cluster.clear();
			cluster.push_back(vec2D[combinedData[i].first]);
		}
    }
	sameEntities.push_back(cluster);
	cout << "total true clusters: " << sameEntities.size() << endl;
}

void getDistances() {
	for(int i = 0; i< sameEntities.size(); i++) {
		bool printThis = false;
		if (sameEntities[i].size()>1)
		{
			cout<< "Cluster " << i << " has " << sameEntities[i].size() << " records" << endl;
			printThis = true;
		} else {
			printThis = false;
		}
		
		
		vector<vector<int>> intraClusterDists;
		for (int j = 0; j < sameEntities[i].size()-1; j++)
		{
			for (int k = j+1 ; k < sameEntities[i].size(); k++)
			{
				vector<int> pairwiseAttrDists;
				for (int attr = 1; attr < attributes; attr++)
				{
					int attrDist = calculateBasicED(sameEntities[i][j][attr], sameEntities[i][k][attr], 10);
					pairwiseAttrDists.push_back(attrDist);
				}
				// if (printThis)
				// {
				// 	for ( int z = 0; z < attributes; z++)
				// 	{
				// 		cout<< " " << pairwiseAttrDists[z];
				// 	}
				// 	cout<< endl;
				// 	cout<< "Attributes: " << endl;
				// 	// for ( int z = 0; z < attributes; z++)
				// 	// {
				// 	// 	cout << sameEntities[i][k][z] << "	 " << sameEntities[i][j][z]  << endl;
				// 	// }
				// 	// cout<< "Suspicious attribute " << sameEntities[i][j][attributes-1] << endl;
				// 	// cout<< "Suspicious attribute " << sameEntities[i][k][attributes-1] << endl;
				// 	cout<< endl;
					
				// }
				
				intraClusterDists.push_back(pairwiseAttrDists);
			}
		}
		distances.push_back(intraClusterDists);
	}
}



int main(int argc, char** argv) {

	// IO PATHS
    string filePath = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/";
    string fileName = argv[1];
    filePath = filePath + argv[1];
	getFormattedDataFromCSV(filePath);
	totalRecords = vec2D.size();
	getCombinedData();

	// Outputs
    string out_file_path = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/data/";
	string stat_file_name = "stat_"+ fileName + "_NC_Analysis";

	// Sort the Combined Data
	clock_t currTS_p0	= clock();
    radixSort(combinedData);
    double sorting_p0_t	= (double)(clock() - currTS_p0) / CLOCKS_PER_SEC;
    cout<< "Sorting time "<< sorting_p0_t << endl;
	// printSortedRecords();
	// Get same Records
	clock_t currTS_p1	= clock();
    getSameEntities();
	double exactClustering_p1_t	= (double)(clock() - currTS_p1) / CLOCKS_PER_SEC;
    cout<< "Same Entity Extraction time: "<< exactClustering_p1_t << endl;
	
	clock_t currTS_p2	= clock();
	getDistances();
	double distanceCalculation_p2_t	= (double)(clock() - currTS_p2) / CLOCKS_PER_SEC;
    cout<< "Intra Cluster Distance Calculation time: "<< distanceCalculation_p2_t << endl;
	string s = "505183";
	string t = "5051830";
	if (!s.compare(t)) {
		cout << "yeah same" << endl;
	}
    return 0;
}