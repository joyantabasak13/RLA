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
#include <queue>

#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

using namespace std;

int threshold = 99;
int totalRecords;
int lenMax;
int attributes;
int cumulativeDistanceThreshold = 7;
vector<int> attrDistThreshold{3,3,3,3};
int base = 36;
int kmer = 3;
int prelen = 1;
long long int blockIDRange;
double fullCheckClockTime = 0.0;
double sortingClockTime = 0.0;

vector<vector<string> > vec2D;
vector<vector<vector<string> > > sameEntities;
vector<pair<int, string> > combinedData;
vector<vector<pair<int, int>>> allEntityBlockRecPairs;
vector<pair<int, int>> allRecordBlockRecPairs;

class UnionFind {
  public:
    int numSets;
	vector<int> parentInd;

    UnionFind() { // Constructor with parameters
		numSets = 0;
    }

	void setVariable(int numRecords) {
		this->numSets = numRecords;
      	this->parentInd.resize(numRecords);
		for (int i = 0; i < numRecords; i++)
		{
			this->parentInd[i] = -1;
		}
	}

	int find(int recID) {
		int root = recID;
    	while (parentInd[root] >= 0){
			root = parentInd[root];
		}
		while (recID != root) { 
			int newParent = parentInd[recID]; 
			parentInd[recID] = root; 
			recID = newParent; 
		}
		return root;
	}

	void weightedUnion(int recId1, int recId2) {
		int i = find(recId1); 
		int j = find(recId2); 
		if (i == j) {
			return;
		}
		// make smaller root point to larger one
		if (parentInd[i] > parentInd[j]) { 
			parentInd[j] += parentInd[i]; 
			parentInd[i] = j; 
		}
		else { 
			parentInd[i] += parentInd[j];
			parentInd[j] = i;
		}
		numSets--;
	}

	bool isConnected(int recId1, int recId2) {
		if (find(recId1) == find(recId2)) {
			return true;
		} else {
			return false;
		}
	}

	int getSetCount() {
		return this->numSets;
	}

	void setParent(int recId, int root) {
		this->parentInd[recId] = root;
	}

	int getRootVal(int recId) {
		int root = find(recId);
		return parentInd[root];
	}

	void setRootVal(int recId, int val) {
		int root = find(recId);
		parentInd[root] = val;
	}
};


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
		for(int i = 0; i<result.size(); i++){
			string attrStr = result[i];
			auto last = std::remove_if(attrStr.begin(), attrStr.end(), [](auto ch) {
        								return ::ispunct(ch) || ::iswpunct(ch);
    								});
			attrStr.erase(last, attrStr.end()); //Remove junk left by remove_if() at the end of iterator
			boost::to_lower(attrStr);
			vec.push_back(attrStr);
		}
        vec2D.push_back(vec);
    }
    records.close();
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
			// cout<< "Match: recID: " << recID << " curID " << cur_recID << endl;
			cluster.push_back(vec2D[combinedData[i].first]);
		} else {
			// cout<< "MisMatch: recID: " << recID << " curID " << cur_recID << endl;
			recID = cur_recID;
			sameEntities.push_back(cluster);
			cluster.clear();
			cluster.push_back(vec2D[combinedData[i].first]);
		}
    }
	sameEntities.push_back(cluster);
	cout << "total true clusters: " << sameEntities.size() << endl;
}

string getBlockingString(vector<string>& record) {
	// string blockingStr = record[blockFieldIndex];
	string blockingStr;
	for(int j=0;j<prelen;j++) {	
		for(int i=1; i<attributes-1; i++) {
			if(record[i].size()<=j) {
				// blockingStr+= "a";
			} else {
				blockingStr+= record[i][j];
			}
		}
	}
	if (blockingStr.size() == 0) {
		// cout<< "Got empty blockSTR" << endl;
		for (int i = 0; i < prelen*(attributes-2); i++)
		// for (int i = 0; i < kmer; i++)
		{
			blockingStr += "a";
		}
	}
	rotate(blockingStr.begin(), blockingStr.begin() + 1, blockingStr.end());
	return blockingStr;
}

void doStandardKmerBlocking(string& blockingStr, pair<int, int> recID, vector<pair<int, pair<int, int>>>& blockRecPairList){
	int alphabets = 26;
	long long int blockID;
	pair<int,pair<int, int>> blockRecPair;
	int strSize = blockingStr.size();
	for (int j = 0; j < strSize ; ++j) {
		blockID = 0;
		for (int k = 0; k < kmer; ++k)
		{
			int blockStrCharCode = (int)blockingStr[(j+k)%strSize];
			if ((blockStrCharCode >= 97) && (blockStrCharCode <= 122)) {
				blockID += ((int)blockingStr[(j+k)%strSize] - 97) * pow(base,k);
			} else if ((blockStrCharCode >= 48) && (blockStrCharCode <= 57)){
				blockID += ((int)blockingStr[(j+k)%strSize] - 48 + alphabets) * pow(base,k);
			}
		}
		blockRecPair.first = blockID;
		blockRecPair.second = recID;
		blockRecPairList.push_back(blockRecPair);
	}

	// cout<< "BlockingString: " << blockingStr << " recID: " << recID.second << " BlockID: " << blockID << endl;
}

void doSuperKmerBlocking(string& blockingStr, pair<int, int> recID, vector<pair<int, pair<int, int>>>& blockRecPairList){
	long long int perAlphaBlocks = pow(base,kmer);
	int alphabets = 26;
	long long int blockID;
	pair<int,pair<int, int>> blockRecPair;
	int strSize = blockingStr.size();
	for (int j = 0; j < strSize ; ++j) {
		blockID = 0;
		for (int k = 0; k < kmer; ++k)
		{
			int blockStrCharCode = (int)blockingStr[(j+k)%strSize];
			if ((blockStrCharCode >= 97) && (blockStrCharCode <= 122)) {
				blockID += ((int)blockingStr[(j+k)%strSize] - 97) * pow(base,k);
			} else if ((blockStrCharCode >= 48) && (blockStrCharCode <= 57)){
				blockID += ((int)blockingStr[(j+k)%strSize] - 48 + alphabets) * pow(base,k);
			}
		}
		if ((blockingStr[0] >= 97) && (blockingStr[0] <= 122)) {
			blockID = ((int)blockingStr[0] - 97)*perAlphaBlocks + blockID;
		} else if ((blockingStr[0] >= 48) && (blockingStr[0] <= 57)){
			blockID = ((int)blockingStr[0] - 48 + alphabets)*perAlphaBlocks + blockID;
		} else {
			blockID = 0;
		}

		blockRecPair.first = blockID;
		blockRecPair.second = recID;
		blockRecPairList.push_back(blockRecPair);
	}

	// cout<< "BlockingString: " << blockingStr << " recID: " << recID.second << " BlockID: " << blockID << endl;
}


void doExactBlocking(string& blockingStr, pair<int, int> recID, vector<pair<int, pair<int, int>>>& blockRecPairList){
	int alphabets = 26;
	long long int blockID = 0;
	pair<int,pair<int, int>> blockRecPair;

	for(int i = 0; i< prelen*(attributes-2); i++) {
		if ((blockingStr[i] >= 97) && (blockingStr[i] <= 122)) {
			blockID += ((int)blockingStr[i] - 97)* pow(base,i);
		} else if ((blockingStr[i] >= 48) && (blockingStr[i] <= 57)){
			blockID += ((int)blockingStr[i] - 48 + alphabets)* pow(base,i);
		}
	}
	
	if (blockID>blockIDRange) {
		cout<< "Invalid BlockID: " << blockID << " BlockIDRange: " << blockIDRange << endl;
		// exit();
	}
	// cout<< "BlockingString: " << blockingStr << " recID: " << recID.second << " BlockID: " << blockID << endl;
	blockRecPair.first = blockID;
	blockRecPair.second = recID;
	blockRecPairList.push_back(blockRecPair);
}

void getBlockRecPairs(pair<pair<int, int>, vector<string>>& record, vector<pair<int, pair<int,int>>>& blockRecPairList) {
	blockIDRange = pow(base,kmer);
	string blockStr = getBlockingString(record.second);
	// doExactBlocking(blockStr, record.first, blockRecPairList);
	doStandardKmerBlocking(blockStr, record.first, blockRecPairList);
	// doSuperKmerBlocking(blockStr, record.first, blockRecPairList);
}


void getAllBlockRecPairs() {
	long long int recID = 0;

	for (int i = 0; i < sameEntities.size(); i++)
	{
		vector<pair<int, pair<int,int>>> thisEntityBlockRecPairs;

		for (int j = 0; j < sameEntities[i].size(); j++)
		{
			recID++;
			pair<pair<int, int>, vector<string>> record;
			// record.first = recID;
			record.first.first = j;
			record.first.second = recID;
			record.second = sameEntities[i][j];
			getBlockRecPairs(record, thisEntityBlockRecPairs);
		}
		// cout<< "Blocking Done" << endl;
		vector<pair<int, int>> blockRecPairVec;
		for(int j=0; j<thisEntityBlockRecPairs.size(); j++) {
			pair<int, int> globalRecIDPair;
			pair<int, int> localBlockIDPair;
			globalRecIDPair.first =  thisEntityBlockRecPairs[j].first;
			localBlockIDPair.first = thisEntityBlockRecPairs[j].first;
			globalRecIDPair.second = thisEntityBlockRecPairs[j].second.second;
			localBlockIDPair.second = thisEntityBlockRecPairs[j].second.first;
			if(thisEntityBlockRecPairs[j].first > blockIDRange) {
				cout<< "BlockID: " << thisEntityBlockRecPairs[j].first << " Local RecID: " << thisEntityBlockRecPairs[j].second.first << " Global recID: " << thisEntityBlockRecPairs[j].second.second << endl;
			}
			blockRecPairVec.push_back(localBlockIDPair);
			allRecordBlockRecPairs.push_back(globalRecIDPair);
		}
		allEntityBlockRecPairs.push_back(blockRecPairVec);
	}
	cout<< "Total BlockRecPairs: "<< allRecordBlockRecPairs.size() << endl;
}

void bubbleSortBlockingIDArray(vector<pair<int, int>>& blockingIDList) {
	int numRecords = blockingIDList.size();
	pair<int, int> tempPair;
    for (int i = 0; i < numRecords - 1; i++){
        for (int j = 0; j < numRecords - i - 1; j++){
            if (blockingIDList[j].first > blockingIDList[j + 1].first){
                tempPair.first = blockingIDList[j+1].first;
				tempPair.second = blockingIDList[j+1].second;
				blockingIDList[j+1] = blockingIDList[j];
				blockingIDList[j] = tempPair;
				// cout<< "J:" << blockingIDList[j].first << " J+1: "<< blockingIDList[j+1].first << endl;
			}
		}
	}
}

void sortBlockingIDArray(vector<pair<int, int>>& blockingIDList) {
	int numRecords = blockingIDList.size();
	vector<pair<int, int>> tempArr(numRecords);
	vector<int> countArr(blockIDRange, 0);
	// cout<< "Started Sorting" << endl;
	for (int j = 0; j < numRecords; ++j) {
		if(blockingIDList[j].first>=blockIDRange) {
			cout<< "Invalid blockID: " << blockingIDList[j].first << endl;
		}
		countArr[blockingIDList[j].first]++;
	}
	// cout<< "Not access any problematic index" << endl;
	// Do prefix sum
	for (int k = 1; k < blockIDRange; ++k)
		countArr[k]	+= countArr[k - 1];

	for (int j = numRecords - 1; j >= 0; --j)
		tempArr[--countArr[blockingIDList[j].first]] = blockingIDList[j];
	
	for (int j = 0; j < numRecords; ++j)
		blockingIDList[j] = tempArr[j];
}

void removeRedundentBlockingID(vector<pair<int, int>>& blockingIDList) {
	int numRecords = blockingIDList.size();
	vector<pair<int, int>> tempArr;
	// totalUniqueBlocks = 1;
	// int copy_count = 1;
	tempArr.push_back(blockingIDList[0]);
	for (int i = 1; i<numRecords; i++) {
		// if (blockingIDList[i].first != blockingIDList[i-1].first) {
		// 	totalUniqueBlocks++;
		// }
		if ( ! ((blockingIDList[i].first == blockingIDList[i-1].first) && (blockingIDList[i].second == blockingIDList[i-1].second))) {
			tempArr.push_back(blockingIDList[i]);
			// copy_count++;
		}
	}
	// totalBlockedKmers = copy_count;
	blockingIDList = tempArr;
	// cout << "Total Length: "<< numRecords << " total copies: "<< copy_count << " Total Unique Blocks: "<< totalUniqueBlocks << endl;
}

bool isLinkageOk(vector<string> &a, vector<string> &b)
{
	// This condition is for the dataset under investigation only
	int cumulativeDist = 0;
	for (int i = 1; i < attributes-1; i++)
	{	
		int dist = calculateBasicED(a[i], b[i], attrDistThreshold[i-1]);
		if (dist > attrDistThreshold[i-1]){
			return false;
		} else {
			cumulativeDist += dist;
			if (cumulativeDist > cumulativeDistanceThreshold){
				return false;
			}
		}
	}
	return true;
}

void getEdgesFromBlockedRecords(vector<pair<int, int>>& blockRecPairs, UnionFind& uf, int clusterID) {
	for (int i = 0; i < blockRecPairs.size() - 1; i++)
	{
		
		int blockID_i = blockRecPairs[i].first;
		int blockID_i_plus_1 = blockRecPairs[i+1].first;
		if(blockID_i != blockID_i_plus_1) {
			continue;
		}
		int recID_i = blockRecPairs[i].second;
		int recID_i_plus_1 = blockRecPairs[i+1].second;
		
		if (!uf.isConnected(recID_i, recID_i_plus_1)) {
			if(isLinkageOk(sameEntities[clusterID][recID_i], sameEntities[clusterID][recID_i_plus_1])) {
				uf.weightedUnion(recID_i, recID_i_plus_1);
			}
		}
	}
}

void doSortedComp(int clusterID, UnionFind& uf){
	vector<pair<int,string> > headlessCopies;
	headlessCopies.resize(2*sameEntities[clusterID].size());
	int blockingStrMaxSize = prelen*(attributes-2);
	string strSample(blockingStrMaxSize, '0');

	for (int i=0; i<sameEntities[clusterID].size(); i++ ) {
		string blockStr = getBlockingString(sameEntities[clusterID][i]);
		int lenDiff	= blockingStrMaxSize - blockStr.length();
		if(lenDiff > 0) {
			blockStr += strSample.substr(0, lenDiff);
		}
		headlessCopies[i].first = i;
		headlessCopies[i].second = blockStr;
		headlessCopies[sameEntities[clusterID].size()+i].first = i;
		headlessCopies[sameEntities[clusterID].size()+i].second = blockStr.substr(1, blockStr.length()-1) + '0';
	}

	lenMax = blockingStrMaxSize;
	radixSort(headlessCopies);

	for (int i = 1; i < headlessCopies.size(); i++) {
		if (headlessCopies[i-1].second.compare(headlessCopies[i].second) == 0) {
			int recID_i = headlessCopies[i-1].first;
            int recID_j = headlessCopies[i].first;
            if (!uf.isConnected(recID_i, recID_j)) {
				uf.weightedUnion(recID_i, recID_j);
			}
		}
	}
}

bool checkIfConnected(vector<pair<int, int>>& blockRecPairs, int clusterID) {
	UnionFind uf;
	uf.setVariable(sameEntities[clusterID].size());
	clock_t currTS_3	= clock();
	bubbleSortBlockingIDArray(blockRecPairs);
	sortingClockTime += (double)(clock() - currTS_3) / CLOCKS_PER_SEC;
	// cout<< "Pairs Sorted: " << blockRecPairs.size() << endl;
	removeRedundentBlockingID(blockRecPairs);
	// cout<< "After Redundent pairs Removal remaining Pairs : " << blockRecPairs.size() << endl;
	getEdgesFromBlockedRecords(blockRecPairs, uf, clusterID);
	// doSortedComp(clusterID, uf);
	// cout<< "Calculating if Cluster is connected or not" << endl;
	int numCluster = uf.getSetCount();
	if(numCluster == 1) {
		// cout<< "Yes cluster is connected" << endl;
		return true;
	} else {
		// cout<< "No cluster is not connected" << endl;
		return false;
	}
}

void calculateClusterCoverage() {
	int coveredClusters = 0;
	int coveredRecords = 0;
	for(int i=0; i< allEntityBlockRecPairs.size(); i++) {
		clock_t currTS_1	= clock();
		bool isConnected = checkIfConnected(allEntityBlockRecPairs[i], i);
		fullCheckClockTime += (double)(clock() - currTS_1) / CLOCKS_PER_SEC;
		// if(i%100000 == 0) {
		// 	cout<< "Clusteres processed: "<< i << endl;
		// 	cout<< "Full checkClockTime: " << fullCheckClockTime << endl;
		// 	cout<< "Sorting Clock Time: " << sortingClockTime << endl;
		// }
		if (isConnected) {
			coveredClusters++;
			coveredRecords += sameEntities[i].size();
		}
	}
	cout<< "Total Clusters: " << allEntityBlockRecPairs.size() << endl;
	cout<< "Covered Clusters: " << coveredClusters << endl;
	cout<< "Covered Record: " << coveredRecords << " %Total: " << ((double)coveredRecords)/((double)vec2D.size());
	cout<< "Coverage: " << (double)(((double)coveredClusters) / ((double)allEntityBlockRecPairs.size())) << endl;
}

void calculateReductionRatio() {
	long long int comparisions = 0;
	sortBlockingIDArray(allRecordBlockRecPairs);
	removeRedundentBlockingID(allRecordBlockRecPairs);
	int start = 0;
	for (int i = 0; i < allRecordBlockRecPairs.size(); i++)
	{
		if(allRecordBlockRecPairs[start].first != allRecordBlockRecPairs[i].first) {
			long long int blockSize = i-start;
			if (blockSize>0) {
				long long int compsNeededforThisBlock = (blockSize*(blockSize-1))/2;
				comparisions += compsNeededforThisBlock;
			}
			start = i;
		}
	}
	long long int blockSize = allRecordBlockRecPairs.size()-start;
	if (blockSize>0) {
		long long int compsNeededforThisBlock = (blockSize*(blockSize-1))/2;
		comparisions += compsNeededforThisBlock;
	}
	double totalComparisions = (double)((double)comparisions) / ((double)1000000000.0);
	cout<< "Total Comps needed: " << totalComparisions << " Billion" << endl;
}

int main(int argc, char** argv) {
	// IO PATHS
    string filePath = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/";
    // string filePath = "/home/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/";

	string fileName = argv[1];
    filePath = filePath + argv[1];
	getFormattedDataFromCSV(filePath);
	totalRecords = vec2D.size();
	getCombinedData();

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
	getAllBlockRecPairs();
	calculateClusterCoverage();
	// calculateReductionRatio();

	double totalTime = (double)(clock() - currTS_p1) / CLOCKS_PER_SEC;
	cout<< "Total time required: " << totalTime << endl;

    return 0;
}