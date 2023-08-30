// Parallel Record Linkage
// Distributed Radix Sort
// MultiCharacter Radix Sort
// Broadcast all records
// Dynamic Load Balancing
// MPI Status Ignore is okay for now, but have to dealt with later

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <ctime> 
#include <fstream>
#include <iostream>
#include <map>
#include <mpi.h>
#include <set>
#include <stdio.h>
#include <string>
#include <sys/time.h>
#include <random>
#include <time.h>
#include <tuple>
#include <unistd.h>
#include <utility>
#include <vector>


using namespace std;
using namespace boost;

#define ROOT_ID				0
#define UNIONFIND_ARR		1
#define KMER				3
#define WORK_COMPLETE		4
#define WORK_ASSINGED		5
#define TEMP_SORT_DATA		6
#define MAX_REC_SIZE		7
#define SHUFFLED_IND_INT	8
#define SHUFFLED_IND_ARR	9

// auto rng = std::default_random_engine {};

int attributes;
int assignedBlocksCount =0;
int base = 36;
int blockingAttrIndex = 1;
int blockIDRange = pow(base,KMER+1);
int characterBlockSize = 3;
int compDividerMultiplier = 100;
int coreRank;
int distanceThreshold = 1;
int l;
int lenMax;
int numCores;
int m;
int perCoreStaticAllocatedChunk = 10;
int threshold = 99;
int totalRecords;
int totalUniqueRecords;
int totalUniqueBlocks;
int world;

long int datasizeAsChar;

long long int totalCompRequired;
long long int totalEdges = 0;

char* superString;
int* assignedRecordIndices;
int* countArr;
int* myIndicesToSort;
int* uniqueRecordIDs;

int matArr[50][50] = {0};

map<int, vector<int> > approxConnectedComponents;

vector<pair<int, int> > blockingIDList;
vector<pair<int, int> > boundaryArr;
vector<pair<int, string> > combinedData;
vector<pair<int, string> > combinedDataToBeSorted;
vector<pair<int,string> > headlessCopies;
vector<pair<int, vector<int>>> preProcessedCombinedData;

vector<string> uniqueCombinedData;
vector<vector<int> > edgeArr;
vector<vector<int> > exactmatches;
vector<vector<int> > finalConnectedComponents;
vector<vector<pair<int, int> > > assignedBlocklists;
vector<vector<pair<int, int> > > assignedRecordDomain;
vector<vector<pair<int, int> > > assignedRecordRange;
vector<vector<string>> uniqueRecords;
vector<vector<string> > vec2D;


class UnionFind {
  public:
    int numSets;
	int* parentInd;

    UnionFind() { // Constructor with parameters
		numSets = 0;
    }

	void setVariable(int numRecords) {
		this->numSets = numRecords;
      	this->parentInd = new int[numRecords];
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

	void setParentArr(int parentArr[]){
		this->parentInd = parentArr;
	}
};

UnionFind uf;

// helps edit distance calculation in calculateBasicED()
int calculateBasicED2(string& str1, string& str2, int threshRem)
{
	int row, col, i, j;

	row		= str1.length() + 1;
	col 	= str2.length() + 1;


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

		return matArr[row - 1][diagonal];

	}
}


void radixSort_basic(vector<pair<int, string> > &strDataArr){
	int numRecords = strDataArr.size();
	vector<pair<int, string> > tempArr(numRecords);
	
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


void radixSort(vector<pair<int, string> > &strDataArr){
	int numRecords = strDataArr.size();
	vector<pair<int, string> > tempArr(numRecords);
	int i = lenMax - 1;
	int range = pow(base,characterBlockSize);
	countArr = new int[range];
	memset(countArr, 0, range * sizeof(int));
	while(i>=0) {
		characterBlockSize = (i-characterBlockSize+1)>=0 ? characterBlockSize: i+1;
		range = pow(base,characterBlockSize);
		memset(countArr, 0, range * sizeof(int));
		for (int j = 0; j < numRecords; ++j) {
			int countInd = 0;
			for(int l=0; l<characterBlockSize; l++){
				int chCode = (strDataArr[j].second)[i-l];
				if ((chCode >= 97) && (chCode <= 122)) {
					countInd += (chCode - 97) * pow(base,l);
				} else if ((chCode >= 48) && (chCode <= 57)){
					countInd += (chCode - 48 + 26) * pow(base,l);
				}
			}
			countArr[countInd]++;
		}
		for (int k = 1; k < range; ++k)
			countArr[k]	+= countArr[k - 1];
		for (int j = numRecords - 1; j >= 0; --j){
			int countInd = 0;
			for(int l=0; l<characterBlockSize; l++){
				int chCode = (strDataArr[j].second)[i-l];
				if ((chCode >= 97) && (chCode <= 122)) {
					countInd += (chCode - 97) * pow(base,l);
				} else if ((chCode >= 48) && (chCode <= 57)){
					countInd += (chCode - 48 + 26) * pow(base,l);
				}
			}
			tempArr[--countArr[countInd]]	= strDataArr[j];
		}
		for (int j = 0; j < numRecords; ++j)
			strDataArr[j]	= tempArr[j];

		i = i-characterBlockSize;
	}

	// for(int i= 0; i<50; i++){
	// 	cout<< strDataArr[i].second << endl;

	// }
}


void radixSort_onPreProcessedData(vector<pair<int, vector<int>>> &strDataArr){
	int numRecords = strDataArr.size();
	int *pointerArr = new int(numRecords);
	int *tempArr = new int(numRecords);
	for(int i=0;i<numRecords;i++){
		pointerArr[i] = i;
	}
	int range = pow(base,characterBlockSize);
	cout<< "Range: " << range << endl;
	countArr = new int[range];
	memset(countArr, 0, range * sizeof(int));
	for (int i = ceil(lenMax/characterBlockSize) - 2; i >= 0; --i) {
		for (int j = 0; j < numRecords; ++j) {
			countArr[strDataArr[pointerArr[j]].second[i]]++;
		}

		for (int k = 1; k < range; ++k)
			countArr[k]	+= countArr[k - 1];

		for (int j = numRecords - 1; j >= 0; --j)
			tempArr[--countArr[strDataArr[pointerArr[j]].second[i]]]	= pointerArr[j];
		
		for(int j= 0; j<numRecords; j++){
			pointerArr[j] = tempArr[j];
		}
		memset(countArr, 0, range * sizeof(int));
		// memset(tempArr, 0, numRecords * sizeof(int));
	}
	vector<pair<int, string>> tempCombData(numRecords);
	for(int i=0; i<numRecords; i++){
		tempCombData[i] = combinedData[pointerArr[i]];
	}
	combinedData = tempCombData;

	// for(int i= 0; i<10; i++){
	// 	cout<< combinedData[i].second << endl;
	// }
	delete countArr;
	delete tempArr;
}


bool isLinkageOk(vector<string> &a, vector<string> &b, int distance)
{
	int singleAttributeDist = 0;
	int cumulativeNonBlockingDist = 0;
	for (int i = 1; i < a.size(); i++)
	{
		singleAttributeDist = calculateBasicED(a[i], b[i], distance);
		if (singleAttributeDist > distance){
			return false;
		}

		cumulativeNonBlockingDist += singleAttributeDist;
		if (cumulativeNonBlockingDist > distance){
			return false;
		}
	}
	return true;
}


double getWallTime() {
	struct timeval time;
    if (gettimeofday(&time,NULL)){
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}


void getFormattedDataFromCSV(string& file_path) {
    string line;
    ifstream records(file_path);
    int ind = 0;
	lenMax = 0;

    while (getline (records, line)) {
        vector<string> result;
        boost::split(result, line, boost::is_any_of(","));
        vector<string> vec;
		int recCharSize = 0;
		for(int i=0; i<result.size(); i++) {
			auto last = std::remove_if(result[i].begin(), result[i].end(), [](auto ch) {
        								return ::ispunct(ch) || ::iswpunct(ch);
    								});
			result[i].erase(last, result[i].end()); //Remove junk left by remove_if() at the end of iterator
			boost::to_lower(result[i]);
			recCharSize += result[i].size();
			vec.push_back(result[i]);
		}
        vec2D.push_back(vec);
		lenMax = recCharSize>lenMax?recCharSize:lenMax;
    }
    records.close();

	attributes = vec2D[0].size();
	totalRecords = vec2D.size();
}


void getCombinedData() {
	string strSample(50, '0');
	combinedData.resize(totalRecords);
	for (int i = 0; i< vec2D.size(); i++) {
		pair<int, string> p;
		p.first = i;
		p.second = vec2D[i][1] + vec2D[i][2] + vec2D[i][3]+vec2D[i][4];
		combinedData[i]=p;
	}

	// Padding to make all characters same size
    for (int i = 0; i < totalRecords; ++i) {
		int lenDiff		= lenMax - combinedData[i].second.length();
		if(lenDiff > 0)
			combinedData[i].second	+= strSample.substr(0, lenDiff);
	}
}


void preprocessCombinedDataForSorting(vector<pair<int, string> > &strDataArr){
	int numInts = ceil(lenMax/characterBlockSize);
	preProcessedCombinedData.resize(totalRecords);
	for(int i=0; i<totalRecords;i++){
		preProcessedCombinedData[i].second.resize(numInts);
		preProcessedCombinedData[i].first = i;
	}
	for(int i=0; i<totalRecords;i++){
		for(int k=0; k<numInts-1;k++){
			int val = 0;
			for(int l=0; l<characterBlockSize; l++){
				int chCode = (strDataArr[i].second)[k*characterBlockSize+l];
				if ((chCode >= 97) && (chCode <= 122)) {
					val += (chCode - 97) * pow(base,characterBlockSize-l-1);
				} else if ((chCode >= 48) && (chCode <= 57)){
					val += (chCode - 48 + 26) * pow(base,characterBlockSize-l-1);
				}
			}
			preProcessedCombinedData[i].second[k] = val;
			if(val>pow(base,characterBlockSize)){
				cout<< "INVALID VAL: " << val << endl;
			}
		}
		int k = numInts-1;
		int lastVal = 0;
		// for(int l=k*characterBlockSize; l<lenMax; l++){
		// 	int chCode = (strDataArr[i].second)[l];
		// 	if ((chCode >= 97) && (chCode <= 122)) {
		// 		lastVal += (chCode - 97) * pow(base,lenMax-l);
		// 	} else if ((chCode >= 48) && (chCode <= 57)){
		// 		lastVal += (chCode - 48 + 26) * pow(base,lenMax-l);
		// 	}
		// }
		preProcessedCombinedData[i].second[k] = lastVal;
	}
}


void getExactMatches() {
	vector<int> tempVec;
	tempVec.push_back(combinedData[0].first);

	for (int i = 1; i < totalRecords; ++i) {
		if(combinedData[i].second.compare(combinedData[i - 1].second) == 0)
			tempVec.push_back(combinedData[i].first);
		else {
			exactmatches.push_back(tempVec);
			tempVec.clear();
			tempVec.push_back(combinedData[i].first);
		}
	}
	exactmatches.push_back(tempVec);
	totalUniqueRecords = exactmatches.size();
	cout << "total exact clusters: " << totalUniqueRecords << endl;
}


void getUniqueEntries() {
	uniqueRecordIDs = new int[totalUniqueRecords];
	uniqueRecords.resize(totalUniqueRecords);
	for(size_t i=0; i<totalUniqueRecords; i++){
		uniqueRecords[i].resize(attributes);
	}

	for (size_t i = 0; i < totalUniqueRecords; i++) {
        uniqueRecordIDs[i] = exactmatches[i][0];
		for(size_t j=0; j<attributes; j++){
			uniqueRecords[i][j] = vec2D[exactmatches[i][0]][j];
		}
    }
}


void getCombinedDataSuperString(){
	datasizeAsChar = (lenMax+1)*totalRecords;
	int curSize = 10;
	int totSize = 10;
	int numDigits = 1;
	// while(totSize<totalRecords){
	// 	datasizeAsChar+=curSize*numDigits;
	// 	totSize = totSize*10;
	// 	curSize = totSize - curSize;
	// 	numDigits++;
	// }
	// curSize = totalRecords - totSize/10;
	// datasizeAsChar+=curSize*numDigits;

	superString = new char[datasizeAsChar];
	int ind = 0;
	int totalTokens = 0;
	for (size_t i = 0; i < totalRecords; i++) {
		// string index = std::to_string(combinedData[i].first);
		// for(size_t j = 0; j< index.size(); j++){
		// 	superString[ind++] = index[j];
		// }
		superString[ind++] = '|';
        for(size_t j = 0; j< lenMax; j++) {
			superString[ind++] = combinedData[i].second[j];
        }
		superString[ind++] = '|';
		totalTokens++;
    }
	cout<< "Total Tokens: " << totalTokens << endl;
}


void getSuperStringSize(){
	datasizeAsChar = 0;
	for (size_t i = 0; i < totalUniqueRecords; i++) {
        for(size_t j = 0; j< attributes; j++) {
            datasizeAsChar += std::strlen(uniqueRecords[i][j].c_str()) +1;
        }
    }
	cout<< "Unique Records Super String Size: " << datasizeAsChar << endl;
}

// Super string is a single string generated 
// cancatenating all attributes of unique records seperated by a special sign 
void getSuperString(){
	getSuperStringSize();
	superString = new char[datasizeAsChar];
	int ind = 0;
	int totalTokens = 0;
	for (size_t i = 0; i < totalUniqueRecords; i++) {
        for(size_t j = 0; j< attributes; j++) {
			for(size_t k=0; k< std::strlen(uniqueRecords[i][j].c_str()); k++){
				superString[ind++] = uniqueRecords[i][j][k];
			}
			superString[ind++] = '|';
			totalTokens++;
        }
    }
}


void populateCombinedRecordsFromSuperString(){
	string superStringStr(superString);
	combinedData.resize(totalRecords);
	vector<string> attrStrs;
    boost::split(attrStrs, superStringStr, boost::is_any_of("|"));
	int ind = 0;

	// split creates an empty token at the end(after last terminating sign)
	for(int i=0; i<attrStrs.size()-1; i++) {
		combinedData[i].first = i;
		combinedData[i].second = attrStrs[i];
		// if((i%2 == 0)){
		// 	combinedData[ind].first = stoi(attrStrs[i]);
		// } else {
		// 	combinedData[ind].second = attrStrs[i];
		// 	ind++;
		// }
	}
}

void populateUniqueRecordsFromSuperString(){
    string superStringStr(superString);
	uniqueRecords.resize(totalUniqueRecords);
	for(size_t i=0; i<totalUniqueRecords; i++){
		uniqueRecords[i].resize(attributes);
	}
	vector<string> attrStrs;
    boost::split(attrStrs, superStringStr, boost::is_any_of("|"));
	int ind = 0;
	// cout<< "Total tokens: " << attrStrs.size() << endl;
	// split creates an empty token at the end(after last terminating sign)
	for(int i=0; i<attrStrs.size()-1; i++) {
		if((i != 0) && (i%attributes == 0)){
			ind++;
		}
		// cout<< "Index accessing: " << ind << " Token: " << attrStrs[i] << " substring no: " << i <<endl;
		uniqueRecords[ind][i%attributes] = attrStrs[i];
	}
}


// Blocking Functions


void getNormalBlockingIDArray() {
	int blockID;
	string blockingStr;
	for (int i = 0; i < totalUniqueRecords; i++) {
		blockingStr = uniqueRecords[i][blockingAttrIndex];
		string temp_str = uniqueRecords[i][blockingAttrIndex];
		while(blockingStr.size() < KMER) {
			if (blockingStr.size() == 0) {
				blockingStr = "a";
				temp_str = "a";
			}
			blockingStr = blockingStr + temp_str;
		}
		int blkstrSize = blockingStr.size();
		for (int j = 0; j < blkstrSize; ++j) {
			blockID = 0;
			for (int k = 0; k < KMER; ++k)
			{
				int blockStrCharCode = (int)blockingStr[(j+k)%blkstrSize];
				if ((blockStrCharCode >= 97) && (blockStrCharCode <= 122)) {
					blockID += (blockStrCharCode - 97) * pow(base,k);
				} else if ((blockStrCharCode >= 48) && (blockStrCharCode <= 57)){
					blockID += (blockStrCharCode - 48 + 26) * pow(base,k);
				}
			}
			pair<int,int> blkRecPair;
			blkRecPair.first = blockID;
			blkRecPair.second = i;
			blockingIDList.push_back(blkRecPair);

			if(blockID<0 || blockID>=blockIDRange) {
				cout<< "Invalid BlockID: " << blockID << endl;
				cout<< "BlockingString: " << blockingStr << endl;
			}
		}	
	}
}


void getSuperBlockingIDArray() {
	int blockID;
	string blockingStr;
	int perAplhaBlocks = pow(base,KMER);
	for (int i = 0; i < totalUniqueRecords; i++) {
		blockingStr = uniqueRecords[i][blockingAttrIndex];
		string temp_str = uniqueRecords[i][blockingAttrIndex];
		while(blockingStr.size() < KMER) {
			if (blockingStr.size() == 0) {
				blockingStr = "a";
				temp_str = "a";
			}
			blockingStr = blockingStr + temp_str;
		}
		int blkstrSize = blockingStr.size();
		for (int j = 0; j < blkstrSize; ++j) {
			blockID = 0;
			for (int k = 0; k < KMER; ++k)
			{
				int blockStrCharCode = (int)blockingStr[(j+k)%blkstrSize];
				if ((blockStrCharCode >= 97) && (blockStrCharCode <= 122)) {
					blockID += (blockStrCharCode - 97) * pow(base,k);
				} else if ((blockStrCharCode >= 48) && (blockStrCharCode <= 57)){
					blockID += (blockStrCharCode - 48 + 26) * pow(base,k);
				}
			}
			pair<int,int> blkRecPair;
			blockID = (blockingStr[0]-97)*perAplhaBlocks + blockID;
			blkRecPair.first = blockID;
			blkRecPair.second = i;
			blockingIDList.push_back(blkRecPair);

			if(blockID<0 || blockID>=blockIDRange) {
				cout<< "Invalid BlockID: " << blockID << endl;
				cout<< "BlockingString: " << blockingStr << endl;
			}
		}	
	}
}


void sortBlockingIDArray() {
	int numRecords = blockingIDList.size();
	vector<pair<int, int>> tempArr(numRecords);
	vector<int> countArr(blockIDRange, 0);
	for (int j = 0; j < numRecords; ++j) {
		countArr[blockingIDList[j].first]++;
	}
	// Do prefix sum
	for (int k = 1; k < blockIDRange; ++k)
		countArr[k]	+= countArr[k - 1];
	for (int j = numRecords - 1; j >= 0; --j)
		tempArr[--countArr[blockingIDList[j].first]] = blockingIDList[j];
	for (int j = 0; j < numRecords; ++j)
		blockingIDList[j] = tempArr[j];
}


void removeRedundentBlockingID() {
	int numRecords = blockingIDList.size();
	vector<pair<int, int>> tempArr;
	totalUniqueBlocks = 1;
	tempArr.push_back(blockingIDList[0]);
	for (int i = 1; i<numRecords; i++) {
		if (blockingIDList[i].first != blockingIDList[i-1].first) {
			totalUniqueBlocks++;
		}
		if ( ! ((blockingIDList[i].first == blockingIDList[i-1].first) && (blockingIDList[i].second == blockingIDList[i-1].second))) {
			tempArr.push_back(blockingIDList[i]);
		}
	}
	blockingIDList = tempArr;

	// cout<< "Total Unique Blocks: "<< totalUniqueBlocks << endl;
	// cout<< "Total unique block-rec pairs: " << numRecords << endl;
	// cout<< "Removed redundant block-rec pairs: " << numRecords - blockingIDList.size() << endl; 
}


// Load Balancing Functions


void randomlyShuffleBlocks(){
	auto rng = std::default_random_engine {};
	std::shuffle(boundaryArr.begin(), boundaryArr.end(),rng);
}


void findBlockBoundaries() {
	//just keep starting ind and range.
	// boundaryArr.resize(totalUniqueBlocks);
	int numRecords = blockingIDList.size();
	totalCompRequired = 0;
	int startInd = 0;
	int range = 0;
	int curBlockInd = 0;
	for (int i = 1; i<numRecords; i++) {
		if (blockingIDList[i].first != blockingIDList[i-1].first) {
			range = i-startInd;
			pair<int,int> p;
			p.first = startInd;
			p.second = range;
			boundaryArr.push_back(p);
			totalCompRequired = totalCompRequired + ceil((range*(range-1))/2);
			curBlockInd++;
			startInd = i;
		}
	}
	// Enter last Block info
	range = numRecords-startInd;
	totalCompRequired = totalCompRequired + ceil((range*(range-1))/2);
	pair<int,int> p;
	p.first = startInd;
	p.second = range;
	boundaryArr.push_back(p);
	// cout<< "Total Unique blocks found: " << curBlockInd << endl;
	// cout<< "Total comp required: " << totalCompRequired << endl;
}


void findBlockAssignments() {
	int totalChunks = (numCores-1)*compDividerMultiplier;
	assignedRecordDomain.resize(totalChunks);
	assignedRecordRange.resize(totalChunks);

	long long int compThreshold = ceil(totalCompRequired/totalChunks);
	long long int curAssignmentSize = 0;
	int recDomainEndInd = 0;
	int remRecCount = 0;
	int curChunk = 0;
	for (int j = 0; j < boundaryArr.size(); j++)
	{
		int range = boundaryArr[j].second;

		long long int curBlockSize = ceil((range*(range-1))/2);

		// assign all remaining comparisions to the last core
		if(curChunk == totalChunks - 1) {
			compThreshold = LLONG_MAX;
		}

		if ((curAssignmentSize+curBlockSize) < compThreshold) {
			assignedRecordDomain[curChunk].push_back(boundaryArr[j]);
			assignedRecordRange[curChunk].push_back(boundaryArr[j]);
			curAssignmentSize = curAssignmentSize + curBlockSize;
		} else {
			remRecCount = boundaryArr[j].second;
			while(curAssignmentSize + remRecCount < compThreshold){
				recDomainEndInd++;
				curAssignmentSize+=remRecCount-1;
				remRecCount--;
			}
			if(recDomainEndInd != 0){
				pair<int,int> dom;
				dom.first = boundaryArr[j].first;
				dom.second = recDomainEndInd;
				assignedRecordDomain[curChunk].push_back(dom);
				assignedRecordRange[curChunk].push_back(boundaryArr[j]);
				pair<int,int> newBound;
				newBound.first = boundaryArr[j].first+recDomainEndInd;
				newBound.second = boundaryArr[j].second-recDomainEndInd;
				boundaryArr.push_back(newBound);
				recDomainEndInd = 0;
				remRecCount = 0;
			} else {
				j--;
			}
			// last core will not reach this condition as no more blocks will remain
			// if(coreRank==curChunk){
			// 	// cout<< endl;
			// 	cout<< "By Rank: "<< coreRank << " coreRank: "<< curChunk << " was assigned: " << curAssignmentSize << " comparisions where threshold is: " << compThreshold << endl;
			// 	// cout<< endl;
			// }
			curChunk++;
			curAssignmentSize = 0;
		}
	}

	// if(coreRank==numCores-1){
	// 	cout<< endl;
	// 	cout<< "By Rank: "<< coreRank << " coreRank: "<< curChunk << " was assigned: " << curAssignmentSize << " comparisions where threshold is: " << compThreshold << endl;
	// 	cout<< endl;
	// }
}


void printWorkChunkSizes(){
	long long int totalCompEst = 0;
	long long int curCompEst = 0;
	long long int minChunk = LLONG_MAX;
	long long int maxChunk = 0;
	int totalChunks = (numCores-1)*compDividerMultiplier;
	cout<< "Total chunks: " << totalChunks << endl;
	for(int c = 0; c < totalChunks; c++){
		curCompEst = 0;
		for(int i = 0; i < assignedRecordDomain[c].size(); i++) {
			int n = assignedRecordDomain[c][i].second;
			curCompEst+= n * assignedRecordRange[c][i].second - ceil(((n-2)*(n-1))/2) ;
		}
		cout<< "Chunk: " << c << " estimated comp: " << curCompEst << endl;
		totalCompEst += curCompEst;
		minChunk = minChunk>curCompEst?curCompEst:minChunk;
		maxChunk = maxChunk<curCompEst?curCompEst:maxChunk;
	}
	cout<<"Loop ended" << endl;
	cout<< "Total Estimated Comp: " << totalCompEst << endl;
	cout<< "Max Estimated Comp: " << maxChunk << endl;
	cout<< "Min Estimated Comp: " << minChunk << endl;
	cout<< "Load Balancing difference: " << maxChunk-minChunk << endl;
}


void doAssignedRecordComp() {
	long long int comparisionsDone = 0;
	long long int totalComparisionAssigned = 0;
	for(int i = 0; i < assignedRecordDomain[coreRank].size(); i++) {
		for(int j=assignedRecordDomain[coreRank][i].first; j < assignedRecordDomain[coreRank][i].first + assignedRecordDomain[coreRank][i].second; j++) {
			for(int k=j+1; k < assignedRecordRange[coreRank][i].first + assignedRecordRange[coreRank][i].second; k++){
				int recid_j = blockingIDList[j].second;
				int recid_k = blockingIDList[k].second;
				totalComparisionAssigned++;
				if(!uf.isConnected(recid_j, recid_k)){
					comparisionsDone++;
				 	if(isLinkageOk(uniqueRecords[recid_j], uniqueRecords[recid_k], distanceThreshold)) {
						uf.weightedUnion(recid_j, recid_k);
				 		totalEdges++;
				 	}
				}
			}
		}
	}
	// cout<< "Processor coreRank: " << coreRank << " found edges: " << totalEdges << endl;
	// cout<< "Rank: " << coreRank << " Out of comparisions:"<< totalComparisionAssigned << " Actual Comparision done: " << comparisionsDone << endl;
}


void doRequestedRecordComp(int requestedChunk) {
	for(int i = 0; i < assignedRecordDomain[requestedChunk].size(); i++) {
		for(int j=assignedRecordDomain[requestedChunk][i].first; j < assignedRecordDomain[requestedChunk][i].first + assignedRecordDomain[requestedChunk][i].second; j++) {
			for(int k=j+1; k < assignedRecordRange[requestedChunk][i].first + assignedRecordRange[requestedChunk][i].second; k++){
				int recid_j = blockingIDList[j].second;
				int recid_k = blockingIDList[k].second;
				if(!uf.isConnected(recid_j, recid_k)){
				 	if(isLinkageOk(uniqueRecords[recid_j], uniqueRecords[recid_k], distanceThreshold)) {
						uf.weightedUnion(recid_j, recid_k);
				 		totalEdges++;
				 	}
				}
			}
		}
	}
}


void mergeEdges(UnionFind &ufTemp) {
	for(int i=0; i< totalUniqueRecords; i++) {
		int root_i = uf.find(i);
		int root_mt = ufTemp.find(i);
		if( root_i != root_mt){
			uf.weightedUnion(root_i, root_mt);
		}
	}
}


void getUniqueCombinedData() {
	string strSample(50, '0');
	uniqueCombinedData.resize(totalUniqueRecords);
	int max = 0;
	for (int i = 0; i< uniqueRecords.size(); i++) {
		// Current hard coded for dataset under investigation
		// Needs to be generalised
		uniqueCombinedData[i] = uniqueRecords[i][1] + uniqueRecords[i][2] + uniqueRecords[i][3]+ uniqueRecords[i][4];
		if (max<uniqueCombinedData[i].size()) {
			max = uniqueCombinedData[i].size();
		}
	}
	lenMax = max;
	// Padding to make all characters same size
    for (int i = 0; i < totalUniqueRecords; ++i) {
		int lenDiff		= lenMax - uniqueCombinedData[i].length();
		if(lenDiff > 0)
			uniqueCombinedData[i]	+= strSample.substr(0, lenDiff);
	}
}


void doSortedComp() {
	getUniqueCombinedData();
	headlessCopies.resize(2*totalUniqueRecords);
	for(int i=0; i< uniqueRecords.size(); i++) {
		headlessCopies[i].first = i;
		headlessCopies[i].second = uniqueCombinedData[i];
		headlessCopies[totalUniqueRecords+i].first = i;
		headlessCopies[totalUniqueRecords+i].second = uniqueCombinedData[i].substr(1,lenMax) + '0';
	}
	// cout<< "Got Headless Copies" << endl;
	radixSort(headlessCopies);
	int extraEdges = 0;
	for (int i = 1; i < headlessCopies.size(); i++) {
		if (headlessCopies[i-1].second.compare(headlessCopies[i].second) == 0) {
			int recID_i = headlessCopies[i-1].first;
            int recID_j = headlessCopies[i].first;
            if (!uf.isConnected(recID_i, recID_j)) {
                uf.weightedUnion(recID_i, recID_j);
				extraEdges++;
			}
		}
	}
	cout<< "Edges Added by Sorting: "<< extraEdges << endl;
}


void findConnComp()
{
	int root, edgeTotal;

    for (int i = 0; i < totalUniqueRecords; i++)
	{
		root = uf.find(i);
		approxConnectedComponents[root].push_back(i);
	}

    cout<< " Single Linkage Connected Components: " << approxConnectedComponents.size()<<endl;
}


void findFinalConnectedComp() {
    int totalNodes = 0;
    int pairsAccessed = 0;
    for (auto const& p : approxConnectedComponents) {
        pairsAccessed++;
        int numComponents = p.second.size();
		totalNodes+=numComponents;
        bool distmat[numComponents][numComponents];
        vector<vector<string>> dataArr(numComponents); // to make cache-efficient, keep records in a row

		// Copy Data in cluster
        for(int c=0; c<p.second.size(); c++) {
            dataArr[c] = uniqueRecords[p.second[c]];
        };

		// generate a 2D matrix filled with all pair comparision results
        for (int i =0; i<numComponents; i++) {
            distmat[i][i] = true;
            for (int j = i+1; j < numComponents; j++)
            {
                if (isLinkageOk(dataArr[i], dataArr[j], distanceThreshold)) {
                    distmat[i][j] = true;
                    distmat[j][i] = true;
                } else {
                    distmat[i][j] = false;
                    distmat[j][i] = false;
                }
            }
        }

        bool nodesConsidered[numComponents];
        for(int i=0; i<numComponents; i++) {
            nodesConsidered[i] = false;
        }

        for(int i=0; i<numComponents; i++) {
            if(nodesConsidered[i] == false) {
                vector<int> connectedComponent;
                connectedComponent.push_back(p.second[i]);
                nodesConsidered[i] = true;
                for(int j=0; j<numComponents; j++) {
                    if ((distmat[i][j] == true) && (nodesConsidered[j]==false)) {
                        connectedComponent.push_back(p.second[j]);
                        nodesConsidered[j] = true;
                    }
                }
                finalConnectedComponents.push_back(connectedComponent);
            }
        }
    }
    // cout<< "Total Nodes: "<< totalNodes << " unique records: " << totalUniqueRecords << endl;
    cout<< "Complete Linkage Connected Components: "<< finalConnectedComponents.size()<<endl;
}


void writeFinalConnectedComponentToFile(string& result_file_name) {
	ofstream out_file;
    out_file.open(result_file_name);

	for(int i = 0; i< finalConnectedComponents.size(); i++) {
        for(int j = 0; j< finalConnectedComponents[i].size(); j++) {
            for(int k=0; k< exactmatches[finalConnectedComponents[i][j]].size(); k++) {
                out_file << uniqueRecords[finalConnectedComponents[i][j]][0] << ",";
			}
		}
        out_file<< "\n";
	}
	out_file.close();
}


int main(int argc, char** argv) {
	// Start timing from here
	clock_t currTS_p0	= clock();

    MPI_Init(NULL, NULL);

    MPI_Comm_rank(MPI_COMM_WORLD, &coreRank);
    MPI_Comm_size(MPI_COMM_WORLD, &world);
	string filePath = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/firstName_LastName_DS/";
	// string filePath = "/home/job22011/RecordLinkage/data/";
	// string filePath = "/gpfs/scratchfs1/sar02010/job22011/RecordLinkage/data/";
	string fileName = argv[1];
	filePath = filePath + argv[1];

	numCores = world;
	l = numCores - 1;
	m = numCores - 1;


    if (coreRank == ROOT_ID) {
		cout<< "Total Cores: " << world << endl;
		// Read data only on root core;
		// Find out unique records and distribute them to other cores
		getFormattedDataFromCSV(filePath);

        double dataRead_t	= (double)(clock() - currTS_p0) / CLOCKS_PER_SEC;
        cout<< "DataRead time "<< dataRead_t << endl;
		cout<< "Max record size in char: " << lenMax << endl;
		// clock_t currTS_p1	= clock();
		// concatenate all attributes for a record
		getCombinedData();

		// Sort the concatenated records
		cout<< "Sorting With " << characterBlockSize <<" characters" << endl;
		clock_t currTS_p101	= clock();
        preprocessCombinedDataForSorting(combinedData);
		// radixSort_onPreProcessedData(preProcessedCombinedData);
        double sorting_t1	= (double)(clock() - currTS_p101) / CLOCKS_PER_SEC;
        cout<< "Sorting time "<< sorting_t1 << endl;
	}

	MPI_Finalize();
	return 0;

	if(coreRank == ROOT_ID) { 
		clock_t currTS_p2	= clock();
		// Get Unique Records
        getExactMatches();
        getUniqueEntries();
        double exactClustering_t	= (double)(clock() - currTS_p2) / CLOCKS_PER_SEC;
        cout<< "De-duplication Time "<< exactClustering_t << endl;
		clock_t currTS_p3	= clock();
		getSuperString();
		double superStrGen_t	= (double)(clock() - currTS_p3) / CLOCKS_PER_SEC;
        cout<< "Super String Generation Time "<< superStrGen_t << endl;
	}
	clock_t currTS_p4	= clock();
	MPI_Bcast(&totalUniqueRecords, 1, MPI_INT, ROOT_ID, MPI_COMM_WORLD);
	MPI_Bcast(&attributes, 1, MPI_INT, ROOT_ID, MPI_COMM_WORLD);
	MPI_Bcast(&datasizeAsChar, 1, MPI_LONG, ROOT_ID, MPI_COMM_WORLD);
	if(coreRank != ROOT_ID) {
		superString = new char[datasizeAsChar];
	}
	MPI_Bcast(superString, datasizeAsChar, MPI_CHAR, ROOT_ID, MPI_COMM_WORLD);

	if(coreRank == ROOT_ID) {
		double broadcasting_t	= (double)(clock() - currTS_p4) / CLOCKS_PER_SEC;
        cout<< "Master Broadcasting Time "<< broadcasting_t << endl;
	}

	uf.setVariable(totalUniqueRecords);

	if(coreRank!= ROOT_ID){
		clock_t currTS_temp;
		clock_t currTS_p5	= clock();
		populateUniqueRecordsFromSuperString();
		// cout<< "Unique records populated by rank: " << coreRank << endl;
		double superStringDecode_t	= (double)(clock() - currTS_p5) / CLOCKS_PER_SEC;
        cout<< "Rank: " << coreRank << " superStringDecoding Time "<< superStringDecode_t << endl;
		// Get Blocking ID Array
		clock_t currTS_p6	= clock();
		getSuperBlockingIDArray();
		double blockingArray_t	= (double)(clock() - currTS_p6) / CLOCKS_PER_SEC;
		cout<< "Rank: " << coreRank << " Blocking Time "<< blockingArray_t << endl;

		// Get Sorted Blocking ID Array
		clock_t currTS_p7	= clock();
		sortBlockingIDArray();
		// Remove Redundent blocks [same <blockID, recID>]
		removeRedundentBlockingID();
		// Find Block assignments
		findBlockBoundaries();
		findBlockAssignments();

		double loadBalancing_t	= (double)(clock() - currTS_p7) / CLOCKS_PER_SEC;
		cout<< "Rank: " << coreRank << " Load Balancing Time "<< loadBalancing_t << endl;
		int chunkProcessed = perCoreStaticAllocatedChunk;
		// printWorkChunkSizes();
		clock_t currTS_p8	= clock();
		for(int i=0; i<perCoreStaticAllocatedChunk; i++){
			int selectedChunk = (coreRank-1)*perCoreStaticAllocatedChunk + i;
			doRequestedRecordComp(selectedChunk);
		}
		auto timeNow = std::chrono::system_clock::now();
		std::time_t cTimeNow = std::chrono::system_clock::to_time_t(timeNow);
		// cout<< "Rank: " << coreRank << " finished static work at: "<< std::ctime(&cTimeNow)  << endl;
		int selectedChunk;
		while(selectedChunk != -1){
			// cout<< "Rank: " << coreRank << " will send a workDone msg" << endl;
			MPI_Send(&coreRank, 1, MPI_INT, ROOT_ID, WORK_COMPLETE, MPI_COMM_WORLD);
			MPI_Recv(&selectedChunk, 1, MPI_INT, ROOT_ID, WORK_ASSINGED, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			// cout<< "Rank: " << coreRank << " got new assigned work of chunk: " << selectedChunk << endl;
			if(selectedChunk != -1){
				// cout<< "Rank: " << coreRank << " received work req chunk: " << selectedChunk << endl;
				doRequestedRecordComp(selectedChunk);
				chunkProcessed++;
				timeNow = std::chrono::system_clock::now();
				cTimeNow = std::chrono::system_clock::to_time_t(timeNow);
				// cout<< "Rank: " << coreRank << " done processing chunk: " << selectedChunk << " at: " << std::ctime(&cTimeNow) << endl;
				currTS_temp	= clock();
			}
		}
		double comp_t	= (double)(clock() - currTS_p8) / CLOCKS_PER_SEC;
		double idle_t	= (double)(clock() - currTS_temp) / CLOCKS_PER_SEC;
		cout<< "Rank: "<< coreRank << " Processed: " << chunkProcessed << " Record Comparision Time "<< comp_t << " Idle time: " << idle_t<< endl;
	}

	if(coreRank == ROOT_ID){
		int received=-1;
		int assignedChunk = (numCores-1)*perCoreStaticAllocatedChunk;
		MPI_Request requests[numCores-1];
		for(int i = 0; i< numCores-1; i++){
			MPI_Irecv(&received, 1, MPI_INT, i+1, WORK_COMPLETE, MPI_COMM_WORLD, &requests[i]);
		}
		while(assignedChunk < (numCores-1)*compDividerMultiplier){
			int index_count;
            int indices[numCores-1];
			// cout<< endl << "Start waiting..." << endl;
            MPI_Waitsome(numCores-1, requests, &index_count, indices, MPI_STATUSES_IGNORE);
			// cout<< "Root got reply from num of nodes: " << index_count << endl;
            for(int i = 0; i < index_count; i++)
            {
				// cout<< "Root replying to core: " << indices[i]+1 << endl;
				if(assignedChunk < (numCores-1)*compDividerMultiplier){
					MPI_Send(&assignedChunk, 1, MPI_INT, indices[i]+1, WORK_ASSINGED, MPI_COMM_WORLD);
					// cout<< "Assigned chunk: " << assignedChunk << " to rank: " << indices[i]+1 << endl;
					MPI_Irecv(&received, 1, MPI_INT, indices[i]+1, WORK_COMPLETE, MPI_COMM_WORLD, &requests[indices[i]]);
					// cout<< "Posted asynchronus receive for rank: " << indices[i]+1 << " with req handle: " <<  indices[i] << endl;
					assignedChunk++;
				}
            }
		}
		// cout<< "while loop terminated" << endl;
		// cout<< endl << "Start waiting for ALL to finish" << endl;
        MPI_Waitall(numCores-1, requests, MPI_STATUSES_IGNORE);
		int endFlag = -1;
		for(int i=1; i<numCores; i++ ){
			MPI_Send(&endFlag, 1, MPI_INT, i, WORK_ASSINGED, MPI_COMM_WORLD);
		}
	}

	// cout<< "Rank: " << coreRank << " found edges: " << totalEdges << endl;

	// double comp_t	= (double)(clock() - currTS_p8) / CLOCKS_PER_SEC;
	// cout<< "Rank: "<< coreRank << " Record Comparision Time "<< comp_t << endl;


	clock_t currTS_p9	= clock();
	if(coreRank != ROOT_ID) {
		MPI_Send(uf.parentInd, totalUniqueRecords, MPI_INT, 0, UNIONFIND_ARR, MPI_COMM_WORLD);
	}

	if(coreRank == ROOT_ID) {
		if(world>1) {
			for(int i=1; i<numCores; i++){
				int tempArr[totalUniqueRecords];
				MPI_Recv(tempArr, totalUniqueRecords, MPI_INT, i, UNIONFIND_ARR, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				// Merge UnionFind arr
				UnionFind ufTemp;
				ufTemp.setVariable(totalUniqueRecords);
				ufTemp.setParentArr(tempArr);
				mergeEdges(ufTemp);
			}
		}
		double merge_t	= (double)(clock() - currTS_p9) / CLOCKS_PER_SEC;
		cout<< "UnionFind Merging Time "<< merge_t << endl;

		clock_t currTS_p10	= clock();
		doSortedComp();
		double sortedComp_t	= (double)(clock() - currTS_p10) / CLOCKS_PER_SEC;
		cout<< "Superblocking Sorting Time "<< sortedComp_t << endl;

		clock_t currTS_p11	= clock();
		findConnComp();

		findFinalConnectedComp();

		double connCompFind_t	= (double)(clock() - currTS_p11) / CLOCKS_PER_SEC;
		cout<< "Complete Linkage Find Time "<< connCompFind_t << endl;

		clock_t currTS_p12	= clock();
		string outname = "Out_" + std::to_string(world) + "_"+ argv[1] + ".csv";
		writeFinalConnectedComponentToFile(outname);
		double fileWriteTime	= (double)(clock() - currTS_p12) / CLOCKS_PER_SEC;
		cout<< "Data Write Time "<< fileWriteTime << endl;
	}

	// MPI_Barrier(MPI_COMM_WORLD);
	if(coreRank == ROOT_ID){
		double total_t	= (double)(clock() - currTS_p0) / CLOCKS_PER_SEC;
		cout<< "TOTAL TIME: "<< total_t << endl;
	}
	MPI_Finalize();
	return 0;
}