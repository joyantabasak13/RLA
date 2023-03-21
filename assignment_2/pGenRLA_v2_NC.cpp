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
int blockingDistanceThreshold = 2;
int blockingDistanceThresholdForClusteredRecords = 3;
int singleNonBlockingDistanceThreshold = 5;
int cumulativeNonBlockingDistanceThreshold = 10;
int clusterSizeThreshold = 1;
int totalRecords;
int lenMax;
int totalUniqueRecords;
int totalUniqueBlocks;
int totalBlockedKmers;
int attributes;
int base = 36;
int kmer = 3;
int blockIDRange = pow(base,kmer+1);
int blockIDRangeForClusteredRecords = pow(base,kmer);
int extraEdges = 0;
int numThreads = 6;
int clusterSizeCutoff = 1;
int blockFieldIndex = 4;
long long int totalCompRequired;
std::mutex mtx;

vector<int> recordsInSmallClusters;
vector<int> recordsInLargeClusters;
vector<vector<int> > exactmatches;
map<int, vector<int> > approxConnectedComponents;
vector<vector<int> > finalConnectedComponents;
map<int, vector<int> > approxConnectedComponentsOnClusteredRecords;
vector<vector<int> > finalConnectedComponentsOnClusteredRecords;
vector<string> vec1D;
vector<vector<string> > vec2D;
vector<vector<int> > clusterExactIndArr;
vector<pair<int,string> > uniqueRecords;
vector<pair<int,string> > headlessCopies;
vector<pair<int, string> > combinedData;
vector<pair<int,int> > blockingIDList;
vector<pair<int,int> > blockingIDListForRecordsInLargeClusters;
vector<pair<int,int> > blockingIDListForRecordsInSmallClusters;
vector<pair<int, int>> boundaryArrForLCBlocks;
vector<pair<int, int>> boundaryArrForSCBlocks;
vector<pair<int, int>> boundaryArr;
vector<vector<pair<int, int>>> assignedBlocklists;
vector<vector<int> > edgeArr;

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

vector<UnionFind> uf;
UnionFind uf_finalClusters;

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
        if ((combinedData[i].second[0] == 'p') && (combinedData[i].second[1]== 'a')){
            cout<< combinedData[i].first << " " << combinedData[i].second << endl;
        }
    }
}

int getKmerCount() {
	long long int totalKmerCount = 0;
	string blockingStr;
	for (int i = 0; i < totalUniqueRecords; i++) {
		totalKmerCount += vec1D[uniqueRecords[i].first*attributes + 1].size() - kmer + 1;
	}
	return totalKmerCount;
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

double getWallTime() {
	struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

// Debug Functions

void printRecordsInSmallClusters(){
	for(int i=0; i<recordsInSmallClusters.size(); i++){
		cout<< vec2D[uniqueRecords[recordsInSmallClusters[i]].first][0] << endl;
	}
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
		
		auto last = std::remove_if(blockStr.begin(), blockStr.end(), [](auto ch) {
        								return ::ispunct(ch) || ::iswpunct(ch);
    								});
		blockStr.erase(last, blockStr.end()); //Remove junk left by remove_if() at the end of iterator
        boost::to_lower(blockStr);
		vec.push_back(blockStr);
		for(int i = 1; i<result.size(); i++){
			if (i!=blockFieldIndex) {
				boost::to_lower(result[i]);
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

void writeApproximateConnectedComponentToFile(string& result_file_name, string& index_file_name) {
	ofstream out_file;
	ofstream index_file;
    out_file.open(result_file_name);
	index_file.open(index_file_name);
	for (auto const& p : approxConnectedComponents) {
        for (int i=0; i<p.second.size(); i++) {
			for (int j=0; j<exactmatches[p.second[i]].size(); j++) {
				out_file<< vec2D[combinedData[exactmatches[p.second[i]][j]].first][0] << ",";
				index_file<< combinedData[exactmatches[p.second[i]][j]].first << ",";
			}
        }
        out_file<< "\n";
		index_file<< "\n";
	}
	out_file.close();
}

void writeFinalConnectedComponentToFile(string& result_file_name) {
	ofstream out_file;
    out_file.open(result_file_name);

	for(int i = 0; i< finalConnectedComponents.size(); i++) {
        for(int j = 0; j< finalConnectedComponents[i].size(); j++) {
            for(int k=0; k< exactmatches[finalConnectedComponents[i][j]].size(); k++) {
                // out_file << vec2D[uniqueRecords[finalConnectedComponents[i][j]].first][0] << ",";
            	out_file<< vec2D[combinedData[exactmatches[finalConnectedComponents[i][j]][k]].first][0] << ",";
			}
		}
        out_file<< "\n";
	}
	out_file.close();
}

void writeFinalConnectedComponentOnClusteredToFile(string& result_file_name) {
	ofstream out_file;
    out_file.open(result_file_name);

	for(int i = 0; i< finalConnectedComponentsOnClusteredRecords.size(); i++) {
        for(int j = 0; j< finalConnectedComponentsOnClusteredRecords[i].size(); j++) {
                //cout<< vec2D[exactmatches[finalConnectedComponents[i][j]][0]][0] << endl;
            out_file << vec2D[uniqueRecords[finalConnectedComponentsOnClusteredRecords[i][j]].first][0] << ",";
		}
        out_file<< "\n";
	}
	out_file.close();
}


// Deduplication and Data preprocessing

void getCombinedData() {
	string strSample(50, '0');
	combinedData.resize(totalRecords);
	int max = 0;
	for (int i = 0; i< vec2D.size(); i++) {
		pair<int, string> p;
		p.first = i;
		for(int j = 1; j<attributes; j++ ) {
			p.second = p.second + vec2D[i][j];
		}
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

// Do exact clustering from lexically sorted vector of int,string pair
void getExactMatches() {
	vector<int> tempVec;

	tempVec.push_back(0);

	for (int i = 1; i < totalRecords; ++i) {
		// if(combinedData[i].second.compare(combinedData[i - 1].second) == 0)
		// 	tempVec.push_back(i);
		// else {
			exactmatches.push_back(tempVec);
			tempVec.clear();
			tempVec.push_back(i);
		// }
	}
	exactmatches.push_back(tempVec);
	totalUniqueRecords = exactmatches.size();
	cout << "total exact clusters: " << totalUniqueRecords << endl;
}

void getUniqueEntries() {
	uniqueRecords.resize(totalUniqueRecords);

	for (size_t i = 0; i < totalUniqueRecords; i++)
    {
        uniqueRecords[i] = combinedData[exactmatches[i][0]];
        // cout<< uniqueRecords[i].first << "\t" << uniqueRecords[i].second << endl;
    }
}

// Blocking Functions

void getBlockingIDArray() {
	int perAlphaBlocks = pow(base,kmer);
	int alphabets = 26;
	int numericals = 10;
	int ind = 0;
	int blockID = 0;
	int indATUnique = 0;
	string blockingStr;
	for (int i = 0; i < totalUniqueRecords; i++) {
		indATUnique = i;
		blockingStr = vec1D[uniqueRecords[i].first*attributes + 1];
		string temp_str = vec1D[uniqueRecords[i].first*attributes + 1];
		while(blockingStr.size() < kmer) {
			if (blockingStr.size() == 0) {
				blockingStr = "a";
				temp_str = "a";
			}
			blockingStr = blockingStr + temp_str;
		}
		for (int j = 0; j < blockingStr.size() - kmer + 1 ; ++j) {
			blockID = 0;
			for (int k = 0; k < kmer; ++k)
			{
				int blockStrCharCode = (int)blockingStr[j+k];
				if ((blockStrCharCode >= 97) && (blockStrCharCode <= 122)) {
					blockID += ((int)blockingStr[j+k] - 97) * pow(base,k);
				} else if ((blockStrCharCode >= 48) && (blockStrCharCode <= 57)){
					blockID += ((int)blockingStr[j+k] - 48 + alphabets) * pow(base,k);
				}
			}
			if ((blockingStr[0] >= 97) && (blockingStr[0] <= 122)) {
					blockID = ((int)blockingStr[0] - 97)*perAlphaBlocks + blockID;
			} else if ((blockingStr[0] >= 48) && (blockingStr[0] <= 57)){
					blockID = ((int)blockingStr[0] - 48 + alphabets)*perAlphaBlocks + blockID;
			} else {
				blockID = 0;
			}
			// blockID = (blockingStr[0]-97)*perAlphaBlocks + blockID;
			pair <int, int> p;
			p.first = blockID;
			p.second = indATUnique;
			blockingIDList.push_back(p);
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
	int copy_count = 1;
	tempArr.push_back(blockingIDList[0]);
	for (int i = 1; i<numRecords; i++) {
		if (blockingIDList[i].first != blockingIDList[i-1].first) {
			totalUniqueBlocks++;
		}
		if ( ! ((blockingIDList[i].first == blockingIDList[i-1].first) && (blockingIDList[i].second == blockingIDList[i-1].second))) {
			tempArr.push_back(blockingIDList[i]);
			copy_count++;
		}
	}
	totalBlockedKmers = copy_count;
	blockingIDList = tempArr;
	cout << "Total Length: "<< numRecords << " total copies: "<< copy_count << " Total Unique Blocks: "<< totalUniqueBlocks << endl;
}

void doSortedComp() {
	headlessCopies.resize(2*totalUniqueRecords);
	int main_tid = numThreads - 1 ;
	for(int i=0; i< uniqueRecords.size(); i++) { 
		headlessCopies[i].first = i;
		headlessCopies[i].second = uniqueRecords[i].second;
		headlessCopies[totalUniqueRecords+i].first = i;
		headlessCopies[totalUniqueRecords+i].second = uniqueRecords[i].second.substr(1,uniqueRecords[i].second.size()-1) + '0';
	}
	cout<< "No SegFault in Copying Headless Copies" << endl;

	radixSort(headlessCopies);

	cout<< "No SegFault in radixSorting Headless Copies" << endl;

	for (int i = 1; i < headlessCopies.size(); i++) {
		if (headlessCopies[i-1].second.compare(headlessCopies[i].second) == 0) {
			int recID_i = headlessCopies[i-1].first;
            int recID_j = headlessCopies[i].first;
            if (!uf[main_tid].isConnected(recID_i, recID_j)) {
                uf[main_tid].weightedUnion(recID_i, recID_j);
				extraEdges++;
			}
		}
	}
	cout<< "No SegFault in edge Adding" << endl;
	cout<< "Edges Added: "<< extraEdges << endl;
}


// Blocking Functions for Clustered Data

void getBlockingIDArrayForRecordsInLargeClusters() {
	int ind = 0;
	int blockID = 0;
	int indATUnique = 0;
	string blockingStr;
	for (int i = 0; i < recordsInLargeClusters.size(); i++) {
		indATUnique = recordsInLargeClusters[i];
		blockingStr = vec1D[uniqueRecords[indATUnique].first*attributes + 1];
		string temp_str = vec1D[uniqueRecords[indATUnique].first*attributes + 1];
		while(blockingStr.size() < kmer) {
			if (blockingStr.size() == 0) {
				blockingStr = "a";
				temp_str = "a";
			}
			blockingStr = blockingStr + temp_str;
		}
		for (int j = 0; j < blockingStr.size() - kmer + 1 ; ++j) {
			blockID = 0;
			for (int k = 0; k < kmer; ++k)
			{
				blockID += ((int)blockingStr[j+k] - 97) * pow(base,k);
			}
			pair <int, int> p;
			p.first = blockID;
			p.second = indATUnique;
			blockingIDListForRecordsInLargeClusters.push_back(p);
		}
	}
}

void getBlockingIDArrayForRecordsInSmallClusters() {
	int ind = 0;
	int blockID = 0;
	int indATUnique = 0;
	string blockingStr;
	for (int i = 0; i < recordsInSmallClusters.size(); i++) {
		indATUnique = recordsInSmallClusters[i];
		blockingStr = vec1D[uniqueRecords[indATUnique].first*attributes + 1];
		string temp_str = vec1D[uniqueRecords[indATUnique].first*attributes + 1];
		while(blockingStr.size() < kmer) {
			if (blockingStr.size() == 0) {
				blockingStr = "a";
				temp_str = "a";
			}
			blockingStr = blockingStr + temp_str;
		}
		for (int j = 0; j < blockingStr.size() - kmer + 1 ; ++j) {
			blockID = 0;
			for (int k = 0; k < kmer; ++k)
			{
				blockID += ((int)blockingStr[j+k] - 97) * pow(base,k);
			}
			pair <int, int> p;
			p.first = blockID;
			p.second = indATUnique;
			blockingIDListForRecordsInSmallClusters.push_back(p);
		}
	}
}

void sortLargeClusterRecordsBlockingIDArray() {
	int numRecords = blockingIDListForRecordsInLargeClusters.size();
	vector<pair<int, int>> tempArr(numRecords);
	vector<int> countArr(blockIDRangeForClusteredRecords, 0);

	for (int j = 0; j < numRecords; ++j) {
		countArr[blockingIDListForRecordsInLargeClusters[j].first]++;
	}
	// Do prefix sum
	for (int k = 1; k < blockIDRangeForClusteredRecords; ++k)
		countArr[k]	+= countArr[k - 1];

	for (int j = numRecords - 1; j >= 0; --j)
		tempArr[--countArr[blockingIDListForRecordsInLargeClusters[j].first]] = blockingIDListForRecordsInLargeClusters[j];
	
	for (int j = 0; j < numRecords; ++j)
		blockingIDListForRecordsInLargeClusters[j] = tempArr[j];
}

void sortSmallClusterRecordsBlockingIDArray() {
	int numRecords = blockingIDListForRecordsInSmallClusters.size();
	vector<pair<int, int>> tempArr(numRecords);
	vector<int> countArr(blockIDRangeForClusteredRecords, 0);

	for (int j = 0; j < numRecords; ++j) {
		countArr[blockingIDListForRecordsInSmallClusters[j].first]++;
	}
	// Do prefix sum
	for (int k = 1; k < blockIDRangeForClusteredRecords; ++k)
		countArr[k]	+= countArr[k - 1];

	for (int j = numRecords - 1; j >= 0; --j)
		tempArr[--countArr[blockingIDListForRecordsInSmallClusters[j].first]] = blockingIDListForRecordsInSmallClusters[j];
	
	for (int j = 0; j < numRecords; ++j)
		blockingIDListForRecordsInSmallClusters[j] = tempArr[j];
}

void removeRedundentLargeClusterRecordsBlockingID() {
	int numRecords = blockingIDListForRecordsInLargeClusters.size();
	vector<pair<int, int>> tempArr;
	totalUniqueBlocks = 1;
	int copy_count = 1;
	tempArr.push_back(blockingIDListForRecordsInLargeClusters[0]);
	for (int i = 1; i<numRecords; i++) {
		if ( ! ((blockingIDListForRecordsInLargeClusters[i].first == blockingIDListForRecordsInLargeClusters[i-1].first) && (blockingIDListForRecordsInLargeClusters[i].second == blockingIDListForRecordsInLargeClusters[i-1].second))) {
			tempArr.push_back(blockingIDListForRecordsInLargeClusters[i]);
			copy_count++;
		}
		if (blockingIDListForRecordsInLargeClusters[i].first != blockingIDListForRecordsInLargeClusters[i-1].first) {
			totalUniqueBlocks++;
		}
	}
	totalBlockedKmers = copy_count;
	blockingIDListForRecordsInLargeClusters = tempArr;
	cout << "For Records in Large Clusters Total Length: "<< numRecords << " total copies: "<< copy_count << " Total Unique Blocks: "<< totalUniqueBlocks << endl;
}

void removeRedundentSmallClusterRecordsBlockingID() {
	int numRecords = blockingIDListForRecordsInSmallClusters.size();
	vector<pair<int, int>> tempArr;
	totalUniqueBlocks = 1;
	int copy_count = 1;
	tempArr.push_back(blockingIDListForRecordsInSmallClusters[0]);
	for (int i = 1; i<numRecords; i++) {
		if ( ! ((blockingIDListForRecordsInSmallClusters[i].first == blockingIDListForRecordsInSmallClusters[i-1].first) && (blockingIDListForRecordsInSmallClusters[i].second == blockingIDListForRecordsInSmallClusters[i-1].second))) {
			tempArr.push_back(blockingIDListForRecordsInSmallClusters[i]);
			copy_count++;
		}
		if (blockingIDListForRecordsInSmallClusters[i].first != blockingIDListForRecordsInSmallClusters[i-1].first) {
			totalUniqueBlocks++;
		}
	}
	totalBlockedKmers = copy_count;
	blockingIDListForRecordsInSmallClusters = tempArr;
	cout << "For Records in Small Clusters Total Length: "<< numRecords << " total copies: "<< copy_count << " Total Unique Blocks: "<< totalUniqueBlocks << endl;
}

void findBlockBoundariesForLargeClusterRecordsBlockList(){
	boundaryArrForLCBlocks.resize(blockIDRangeForClusteredRecords);
	for (int i = 0; i < blockIDRangeForClusteredRecords; i++)
	{
		boundaryArrForLCBlocks[i].first = -1;
		boundaryArrForLCBlocks[i].second = 0;
	}
	
	int numRecords = blockingIDListForRecordsInLargeClusters.size();
	int startInd = 0;
	int range = 0;
	int curBlockId = blockingIDListForRecordsInLargeClusters[0].first;
	for (int i = 1; i<numRecords; i++) {
		if (blockingIDListForRecordsInLargeClusters[i].first != blockingIDListForRecordsInLargeClusters[i-1].first) {
			range = i-startInd;
			boundaryArrForLCBlocks[curBlockId].first = startInd;
			boundaryArrForLCBlocks[curBlockId].second = range;
			curBlockId = blockingIDListForRecordsInLargeClusters[i].first;
			startInd = i;
		}
	}
	// Enter last Block info
	range = numRecords-startInd;
	boundaryArrForLCBlocks[curBlockId].first = startInd;
	boundaryArrForLCBlocks[curBlockId].second = range;
}

void findBlockBoundariesForSmallClusterRecordsBlockList(){
	boundaryArrForSCBlocks.resize(blockIDRangeForClusteredRecords);
	for (int i = 0; i < blockIDRangeForClusteredRecords; i++)
	{
		boundaryArrForSCBlocks[i].first = -1;
		boundaryArrForSCBlocks[i].second = 0;
	}
	
	int numRecords = blockingIDListForRecordsInSmallClusters.size();
	int startInd = 0;
	int range = 0;
	int curBlockId = blockingIDListForRecordsInSmallClusters[0].first;
	for (int i = 1; i<numRecords; i++) {
		if (blockingIDListForRecordsInSmallClusters[i].first != blockingIDListForRecordsInSmallClusters[i-1].first) {
			range = i-startInd;
			boundaryArrForSCBlocks[curBlockId].first = startInd;
			boundaryArrForSCBlocks[curBlockId].second = range;
			curBlockId = blockingIDListForRecordsInSmallClusters[i].first;
			startInd = i;
		}
	}
	// Enter last Block info
	range = numRecords-startInd;
	boundaryArrForSCBlocks[curBlockId].first = startInd;
	boundaryArrForSCBlocks[curBlockId].second = range;
}

// Edge extraction functions for blocked clustered records
void getEdgesFromBlockedRecordsInSmallClusters() {
	//cout<< "Hit small cluster edges extraction" << endl;
	//cout<< "Blocks: "<< boundaryArrForSCBlocks.size() << endl;
	int edgesAdded = 0;
	for (int i = 0; i < boundaryArrForSCBlocks.size(); i++)
	{
		int startInd = boundaryArrForSCBlocks[i].first;
		int endInd = startInd + boundaryArrForSCBlocks[i].second;
		//cout<< "Start: " << startInd << " End: " << endInd << endl;
		if (startInd < 0 || endInd < 1)
		{
			// just skip
			continue;
		}
		
		vector<pair<int, vector<string>>> recordList;
		for (int j = startInd; j < endInd; j++)
		{
			int recID = blockingIDListForRecordsInSmallClusters[j].second;
		//	cout<< "Copying: " << blockingIDListForRecordsInSmallClusters[j].second << endl;
			pair<int, vector<string>> record;
			record.first = recID;
			record.second = vec2D[uniqueRecords[recID].first];
			recordList.push_back(record);
		}
		// cout<< "Records Copied: " << recordList.size() << endl; 
		for (int j = 0; j < recordList.size() - 1; j++)
		{
			pair<int, vector<string>> rec_j = recordList[j];

			for (int k = j+1; k < recordList.size(); k++)
			{
			//	cout<< "Checking for j:" << j << " and k " << k << endl;
				pair<int, vector<string>> rec_k = recordList[k];
				if (!uf_finalClusters.isConnected(rec_j.first, rec_k.first)) {
					if (isLinkageOk(rec_j.second, rec_k.second, blockingDistanceThresholdForClusteredRecords, singleNonBlockingDistanceThreshold , cumulativeNonBlockingDistanceThreshold)) {
						uf_finalClusters.weightedUnion(rec_j.first, rec_k.first);
						edgesAdded++;
					}
				}	
			}
		}
	}
	cout<< "Edges added by comparing records in Small Clusters: " << edgesAdded << endl;
	
}

void getEdgesFromBlockedRecordsInSmallClustersToLargeClusters() {
	int edgesAdded = 0;
	for (int i = 0; i < boundaryArrForSCBlocks.size(); i++)
	{
		int startInd = boundaryArrForSCBlocks[i].first;
		int endInd = startInd + boundaryArrForSCBlocks[i].second;
		if(startInd < 0) {
			continue;
			// just skip this block
		}
		vector<pair<int, vector<string>>> recordList_small;
		for (int j = startInd; j < endInd; j++)
		{
			int recID = blockingIDListForRecordsInSmallClusters[j].second;
			pair<int, vector<string>> record;
			record.first = recID;
			record.second = vec2D[uniqueRecords[recID].first];
			recordList_small.push_back(record);
		}

		int startInd_large = boundaryArrForLCBlocks[i].first;
		int endInd_large = startInd_large + boundaryArrForLCBlocks[i].second;
		if (startInd_large < 0)
		{
			continue;
			// skip this block too. No records to compare to
		}
		
		vector<pair<int, vector<string>>> recordList_large;
		for (int j = startInd_large; j < endInd_large; j++)
		{
			int recID = blockingIDListForRecordsInLargeClusters[j].second;
			pair<int, vector<string>> record;
			record.first = recID;
			record.second = vec2D[uniqueRecords[recID].first];
			recordList_large.push_back(record);
		}
		
		for (int j = 0; j < recordList_small.size(); j++)
		{
			pair<int, vector<string>> rec_j = recordList_small[j];

			for (int k = 0;  k < recordList_large.size(); k++)
			{
				pair<int, vector<string>> rec_k = recordList_large[k];
				if (!uf_finalClusters.isConnected(rec_j.first, rec_k.first)) {
					if (isLinkageOk(rec_j.second, rec_k.second, blockingDistanceThresholdForClusteredRecords, singleNonBlockingDistanceThreshold, cumulativeNonBlockingDistanceThreshold)) {
						uf_finalClusters.weightedUnion(rec_j.first, rec_k.first);
						edgesAdded++;
					}
				}	
			}
		}
	}
	cout<< "Edges added by comparing records in Small to Large Clusters: " << edgesAdded << endl;
}

// Load Balancing Functions

void findBlockBoundaries() {
	//just keep starting ind and range.
	boundaryArr.resize(totalUniqueBlocks);
	int numRecords = blockingIDList.size();
	totalCompRequired = 0;
	int startInd = 0;
	int range = 0;
	int curBlockInd = 0;
	for (int i = 1; i<numRecords; i++) {
		if (blockingIDList[i].first != blockingIDList[i-1].first) {
			range = i-startInd;
			boundaryArr[curBlockInd].first = startInd;
			boundaryArr[curBlockInd].second = range;
			totalCompRequired += ceil((pow(range,2) - range) / 2);
			curBlockInd++;
			startInd = i;
		}
	}
	// Enter last Block info
	range = numRecords-startInd;
	totalCompRequired += ceil((pow(range,2) - range) / 2);
	boundaryArr[curBlockInd].first = startInd;
	boundaryArr[curBlockInd].second = range;
	curBlockInd++;
	cout<< "Total Unique blocks found: " << curBlockInd << endl;
}

void sortByBlockSizes() {
	std::sort(boundaryArr.begin(), boundaryArr.end(), [](auto &left, auto &right) {
    	return left.second < right.second;
	});
}

void findBlockAssignments() {
	assignedBlocklists.resize(numThreads);
	cout << "Total comp required: " << totalCompRequired << endl;
	long long int threshold = (long long int)(totalCompRequired/numThreads);
	long long int curAssignmentSize = 0;
	int lastInd = boundaryArr.size();
	int startInd = -1;
 
	for(int i=0; i<numThreads; i++) {
		curAssignmentSize = 0;
		//cout<< "Thread: "<< i << " was assigned: " << curAssignmentSize << " comparisions where threshold is: " << threshold << endl;
		for (int j = lastInd-1; j > startInd; j--)
		{
			//cout<< "SegFault for j " << j << endl;
			long long int curBlockSize = pow(boundaryArr[j].second,2);
			if ((curAssignmentSize+curBlockSize) < threshold) {
				assignedBlocklists[i].push_back(boundaryArr[j]);
				curAssignmentSize = curAssignmentSize + curBlockSize;
				lastInd = j;
			} else {
				break;
			}
		}
		// cout<< "Thread: "<< i << " was assigned: " << curAssignmentSize << " comparisions where threshold is: " << threshold << endl;
		if (curAssignmentSize < threshold) {
			for(int j = startInd+1; j<lastInd; j++) {
				if (curAssignmentSize < threshold) {
					assignedBlocklists[i].push_back(boundaryArr[j]);
					curAssignmentSize = curAssignmentSize + pow(boundaryArr[j].second,2);
					startInd = j;
				} else {
					break;
				}
			}
		}
		cout<< "Thread: "<< i << " was assigned: " << curAssignmentSize << " comparisions where threshold is: " << threshold << endl;
	}
	cout<< "Start ind: " << startInd << " last ind " << lastInd << endl; 
	if (lastInd - startInd > 1) {
		for (int i = startInd + 1; i<lastInd; i++) {
			assignedBlocklists[i%numThreads].push_back(boundaryArr[i]);
		}
	}
}

// Connected Component Extraction And Revision

//  find clusters as connected components in a graph where edges are connection among records 
//  and vertices are record index
void findConnComp()
{
	int i, root, edgeTotal;
	int main_tid = numThreads - 1 ;
    for (int i = 0; i < totalUniqueRecords; i++)
	{
		root = uf[main_tid].find(i);
		approxConnectedComponents[root].push_back(i);
	}
	
    cout<< "Single Linkage Connected Components: " << approxConnectedComponents.size()<<endl;
}

void findFinalConnectedComp() {
    int totalNodes = 0;
    int pairsAccessed = 0;
	// accessing all <key, value> pair where,
	// key is the root of each approx. cluster
	// value is a list of integers corresponding to unique record indices
    
	for (auto const& p : approxConnectedComponents) {
        pairsAccessed++;
        int numComponents = p.second.size();
		cout<< "Cluster: "<< pairsAccessed << " has " << numComponents << " Components"<< endl;
		
        totalNodes += numComponents;
        bool distmat[numComponents][numComponents];
        vector<vector<string>> dataArr(numComponents); // to make cache-efficient, keep records in a row
		// Copy Data in cluster into a vector
		cout<< "Data copy Started" << endl;
        for(int c=0; c<p.second.size(); c++) {
            dataArr[c] = vec2D[combinedData[exactmatches[p.second[c]][0]].first];
        };
		cout<< "Data copied Successfully" << endl;
		if(pairsAccessed == 126306) {
			for (int m = 0; m < dataArr.size(); m++)
			{
				for (int n = 0; n < dataArr[m].size(); n++)
				{
					cout<< " " << dataArr[m][n] << " ";
				}
				cout<<endl;
			}
			
		}

		// generate a 2D matrix filled with all pair comparision results
        for (int i =0; i<numComponents; i++) {
            distmat[i][i] = true;
            for (int j = i+1; j < numComponents; j++)
            {
                if (isLinkageOk(dataArr[i], dataArr[j], blockingDistanceThreshold, singleNonBlockingDistanceThreshold, cumulativeNonBlockingDistanceThreshold)) {
                    distmat[i][j] = true;
                    distmat[j][i] = true;
                } else {
                    distmat[i][j] = false;
                    distmat[j][i] = false;
                }
            }
        }

		cout<< "Cluster: "<< pairsAccessed << " all pair comparision done"<< endl;

        bool nodesConsidered[numComponents];
        for(int i=0; i<numComponents; i++) {
            nodesConsidered[i] = false;
        }

        for(int i=0; i<numComponents; i++) {
            if(nodesConsidered[i] == false) { // false means its not in a cluster yet
                vector<int> connectedComponentCandidates;
				vector<int> selectedCandidates;
				vector<int> connectedComponent;
                connectedComponent.push_back(p.second[i]);
                nodesConsidered[i] = true;
                for(int j=0; j<numComponents; j++) {
                    if ((distmat[i][j] == true) && (nodesConsidered[j]==false)) {
						connectedComponentCandidates.push_back(j);
					}
				}

				if (connectedComponentCandidates.size() >= 1)
				{
					int triviallySelected = connectedComponentCandidates[connectedComponentCandidates.size()-1];
					selectedCandidates.push_back(triviallySelected);
					connectedComponentCandidates.pop_back();
					connectedComponent.push_back(p.second[triviallySelected]);
                    nodesConsidered[triviallySelected] = true;
					cout<< "Trivial Component Extracted" << endl;
				}
				
				for (int j = 0; j < connectedComponentCandidates.size(); j++)
				{
					bool isSelected = true;
					for (int k = 0; k < selectedCandidates.size(); k++)
					{
						if (distmat[connectedComponentCandidates[j]][selectedCandidates[k]] == false) {
							isSelected = false;
							cout<< "Is Selected is False" << endl;
							break;
						}
					}
					
					if(isSelected) {
						cout<< "Is Selected is true for " << connectedComponentCandidates[j] << endl;
						selectedCandidates.push_back(connectedComponentCandidates[j]);
						connectedComponent.push_back(p.second[connectedComponentCandidates[j]]);
						nodesConsidered[connectedComponentCandidates[j]] = true;
						cout<< "New Connected Component found" << endl;
					}
				}
                finalConnectedComponents.push_back(connectedComponent);
            }
        }
    }
    cout<< "Total Nodes: "<< totalNodes << " unique records: " << totalUniqueRecords << endl;
    cout<< "Complete Linkage Components: "<< finalConnectedComponents.size()<<endl;
}

void findConnCompOnGeneralizedClustering()
{
	int i, root, edgeTotal;
    for (int i = 0; i < totalUniqueRecords; i++)
	{
		root = uf_finalClusters.find(i);
		approxConnectedComponentsOnClusteredRecords[root].push_back(i);
	}
	
    cout<< "Single Linkage Connected Components On Clustered Records: " << approxConnectedComponentsOnClusteredRecords.size()<<endl;
}

void findFinalConnectedCompOnGeneralizedClustering() {
    int totalNodes = 0;
    int pairsAccessed = 0;
    for (auto const& p : approxConnectedComponentsOnClusteredRecords) {
        pairsAccessed++;
        int numComponents = p.second.size();
        totalNodes+=numComponents;
        bool distmat[numComponents][numComponents];
        vector<vector<string>> dataArr(numComponents); // to make cache-efficient, keep records in a row
		// Copy Data in cluster into a vector
        for(int c=0; c<p.second.size(); c++) {
            dataArr[c] = vec2D[uniqueRecords[p.second[c]].first];
        };

		// generate a 2D matrix filled with all pair comparision results
        for (int i =0; i<numComponents; i++) {
            distmat[i][i] = true;
            for (int j = i+1; j < numComponents; j++)
            {
                if (isLinkageOk(dataArr[i], dataArr[j], blockingDistanceThresholdForClusteredRecords, singleNonBlockingDistanceThreshold, cumulativeNonBlockingDistanceThreshold)) {
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
            if(nodesConsidered[i] == false) { // false means its not in a cluster yet
                vector<int> connectedComponentCandidates;
				vector<int> selectedCandidates;
				vector<int> connectedComponent;
                connectedComponent.push_back(p.second[i]);
                nodesConsidered[i] = true;
                for(int j=0; j<numComponents; j++) {
                    if ((distmat[i][j] == true) && (nodesConsidered[j]==false)) {
						connectedComponentCandidates.push_back(j);
					}
				}
				if (connectedComponentCandidates.size() >= 1)
				{
					int triviallySelected = connectedComponentCandidates[connectedComponentCandidates.size()-1];
					connectedComponentCandidates.pop_back();
					selectedCandidates.push_back(triviallySelected);
					connectedComponent.push_back(p.second[triviallySelected]);
                    nodesConsidered[triviallySelected] = true;
				}

				for (int j = 0; j < connectedComponentCandidates.size(); j++)
				{
					bool isSelected = true;
					for (int k = 0; k < selectedCandidates.size(); k++)
					{
						if (distmat[connectedComponentCandidates[j]][selectedCandidates[k]] == false) {
							isSelected = false;
							break;
						}
					}
					if(isSelected) {
						selectedCandidates.push_back(connectedComponentCandidates[j]);
						connectedComponent.push_back(p.second[connectedComponentCandidates[j]]);
						nodesConsidered[connectedComponentCandidates[j]] = true;
					}
					
				}
                finalConnectedComponentsOnClusteredRecords.push_back(connectedComponent);
            }
        }
    }
    cout<< "Total Nodes: "<< totalNodes << " unique records: " << totalUniqueRecords << endl;
    cout<< "Complete Linkage Components On Connected Records: "<< finalConnectedComponentsOnClusteredRecords.size()<<endl;
}


// Reset and rebuild union find sets
void setUnionFindToCompleteClustering(){
	uf_finalClusters.setVariable(totalUniqueRecords);
	for(int i = 0; i< finalConnectedComponents.size(); i++) {
		int root = finalConnectedComponents[i][0];
		for(int j = 1; j< finalConnectedComponents[i].size(); j++) {
			uf_finalClusters.weightedUnion(root,finalConnectedComponents[i][j]);
		}
	}
}

// extract records belonging to clusters
void extractRecordsInClusters() {
	for(int i = 0; i< finalConnectedComponents.size(); i++) {
		if (finalConnectedComponents[i].size() <= clusterSizeThreshold) {
			for (int j = 0; j < finalConnectedComponents[i].size(); j++)
			{
				recordsInSmallClusters.push_back(finalConnectedComponents[i][j]);
			}
			
		} else {
			for (int j = 0; j < finalConnectedComponents[i].size(); j++)
			{
				recordsInLargeClusters.push_back(finalConnectedComponents[i][j]);
			}
		}
	}
	cout<< "Records in small clusters: " << recordsInSmallClusters.size() << endl;
	cout<< "Records in Large Clusters: " << recordsInLargeClusters.size() << endl;
	// printRecordsInSmallClusters();
}

// Threading related Functions

void mergeEdges() {
	int mainTid = numThreads-1;
	for(int i= 0; i<numThreads-1; i++ ) {
		for(int j=0; j< totalUniqueRecords; j++) {
			int root_i = uf[i].find(j);
			int root_mt = uf[mainTid].find(j);
			if( root_i != root_mt){
				uf[mainTid].weightedUnion(root_i, root_mt);
			}
		}
	}
}

void getBlockRecords(pair<int, int> &blockInfo, vector<pair<int, vector<string>>> &recordList) {
	recordList.resize(blockInfo.second);
	for (int i = 0; i < blockInfo.second; i++) {
		recordList[i].first = blockingIDList[blockInfo.first+i].second;
		recordList[i].second = vec2D[uniqueRecords[recordList[i].first].first];
	}
}

void getEdgesFromBlockedRecords(int id, vector<pair<int, vector<string>>> &blockRecords) {
	for (int i = 0; i < blockRecords.size() - 1; i++)
	{
        for (int j = i+1; j < blockRecords.size(); j++)
        {
            int recID_i = blockRecords[i].first;
            int recID_j = blockRecords[j].first;
            if (!uf[id].isConnected(recID_i, recID_j)) {
				if (isLinkageOk(blockRecords[i].second, blockRecords[j].second, blockingDistanceThreshold, singleNonBlockingDistanceThreshold, cumulativeNonBlockingDistanceThreshold)) {
                    uf[id].weightedUnion(recID_i, recID_j);
                }
			}
        }  
	}
}

void doNormalThreadedBlocking(int tID) {
	cout<< "Doing Super Blocking for Thread: " << tID << endl;
	cout<< "Thread: "<< tID << " has " << assignedBlocklists[tID].size() << " blocks" << endl;
	for(int i = 0; i< assignedBlocklists[tID].size(); i++) {
		pair<int, int> block = assignedBlocklists[tID][i];
		vector<pair<int, vector<string>>> blockRecords;
		getBlockRecords(block, blockRecords);
		getEdgesFromBlockedRecords(tID, blockRecords);
	}
	cout << "Thread "<< tID << " Total Edges: "<< totalUniqueRecords - uf[tID].getSetCount() << endl;
}

// main function for threads
void *threadDriver(void* ptr) {
	// int threadID = *static_cast<int*>(ptr);
	int threadID = *static_cast<int*>(ptr);
	cout<< "Thread: " << threadID << " Running" << endl;
	doNormalThreadedBlocking(threadID);
	return 0;
}

int main(int argc, char** argv) {

	// IO PATHS
	// Input
    // string filePath = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/";
    string filePath = "/home/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/";
    string fileName = argv[1];
    filePath = filePath + argv[1];
	getFormattedDataFromCSV(filePath);
	totalRecords = vec2D.size();
	getCombinedData();

	// Outputs
    // string out_file_path = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/data/";
    string out_file_path = "/home/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/data/";
	string fileNameSuffix = "_pGEN_NC_lastName_6_threads_dist_2_5_9";
	string out_name1 = out_file_path + "out_SB_RLA_SingleLinkage_NO_DEDUP_"+ fileName + fileNameSuffix;
	string out_name2 = out_file_path + "out_SB_RLA_CompleteLinkage_NO_DEDUP_"+ fileName + fileNameSuffix;
	string out_name3 = out_file_path + "out_SB_RLA_SingleLinkage_RecInd_NO_DEDUP_"+ fileName + fileNameSuffix;
	string out_name4 = out_file_path + "out_general_RLA_"+ fileName + fileNameSuffix;
	string stat_file_name = "stat_"+ fileName + fileNameSuffix;

	// Sort the Combined Data
	clock_t currTS_p0	= clock();
	double currWallT_p0 = getWallTime();
    radixSort(combinedData);
    double sorting_p0_t	= (double)(clock() - currTS_p0) / CLOCKS_PER_SEC;
    cout<< "Sorting time "<< sorting_p0_t << endl;

	// Get Unique Records
	clock_t currTS_p1	= clock();
	getExactMatches();
    getUniqueEntries();
	double exactClustering_p1_t	= (double)(clock() - currTS_p1) / CLOCKS_PER_SEC;
    cout<< "De-duplication Time "<< exactClustering_p1_t << endl;

	// Get Blocking ID Array
	clock_t currTS_p2	= clock();
	getBlockingIDArray();
	double blockingArray_p2_t	= (double)(clock() - currTS_p2) / CLOCKS_PER_SEC;
    cout<< "Getting Blocking Array Time "<< blockingArray_p2_t << endl;

	// Get Sorted Blocking ID Array
	clock_t currTS_p3	= clock();
	sortBlockingIDArray();
	double sortBlockingArray_p3_t	= (double)(clock() - currTS_p3) / CLOCKS_PER_SEC;
    cout<< "Getting Sorted Blocking Array Time "<< sortBlockingArray_p3_t << endl;

	// Remove Redundent blocks [same <blockID, recID>]
	clock_t currTS_p4	= clock();
	removeRedundentBlockingID();

	// Find Block assignments
	// findBlockBoundaries() and sortByBlockSizes() can be merged. But might not get much speedup.
	findBlockBoundaries();
	sortByBlockSizes();
	findBlockAssignments();
	double loadBalancing_p4_t	= (double)(clock() - currTS_p4) / CLOCKS_PER_SEC;
    cout<< "Get Load Balancing Time "<< loadBalancing_p4_t << endl;

	// Thread Working
	uf.resize(numThreads);
	for (int i=0; i<numThreads; i++) {
		uf[i].setVariable(totalUniqueRecords);
	}

	pthread_t threads[numThreads-1];
	clock_t currTS_p5	= clock();
	double currWallT_p5 = getWallTime();
	for (int i = 0; i < numThreads-1; i++)
	{
		int threadID = static_cast<int>(i);
		int iret = pthread_create(&threads[threadID], NULL, threadDriver, &threadID);
		usleep(100);
	}

	doNormalThreadedBlocking(numThreads-1);

	for (int i = 0; i < numThreads-1; i++)
	{
		int threadID = i;
		pthread_join(threads[threadID], NULL);
	}

	double comparisionsDone_p5_t	= (double)(clock() - currTS_p5) / CLOCKS_PER_SEC;
	double comparisionsDone_p5_Wt = getWallTime();
    cout<< "Get Camparision time done in threads Processor Time "<< comparisionsDone_p5_t << endl;
	cout<< "Get Camparision time done in threads Wall Time "<< (double)(comparisionsDone_p5_Wt - currWallT_p5) << endl;

	// Merge edges
	clock_t currTS_p6	= clock();
	cout<< "No core dump till now" << endl;
	mergeEdges();
	cout<< "No core dump after merging edges" << endl;
	doSortedComp();
	cout<< "No core dump after Sorted Comparisions" << endl;
	double edgeMergingDone_p6_t	= (double)(clock() - currTS_p6) / CLOCKS_PER_SEC;
	cout<< "Get edge Merging Time "<< edgeMergingDone_p6_t << endl;
	cout<< "Total Edges: " << totalUniqueRecords-uf[numThreads-1].getSetCount() << endl;

	// Find Connected components (Single Linkage)
	clock_t currTS_p7	= clock();
    findConnComp();
    double findComp_p7_t	= (double)(clock() - currTS_p7) / CLOCKS_PER_SEC;
    cout<< "Connected Comp Find Time "<< findComp_p7_t << endl;

	writeApproximateConnectedComponentToFile(out_name1, out_name3);

	cout<< "Single Linkage Connected Components are writen to file" << endl;
	// Find Connected components (Complete Linkage)
    clock_t currTS_p8_t	= clock();
    findFinalConnectedComp();
    double findFinalComp_t	= (double)(clock() - currTS_p8_t) / CLOCKS_PER_SEC;
    cout<< "Final Connected Comps Find Time "<< findFinalComp_t << endl;

	// Total Time Required For SuperBlocking
    double total_SB_t	= (double)(clock() - currTS_p0) / CLOCKS_PER_SEC;
	cout<< "Super Blocking: Total processor run time "<< total_SB_t << endl;
	double SB_done_pX_Wt = getWallTime();
	cout<< "SuperBlocking: Get Total Wall Time "<< (double)(SB_done_pX_Wt - currWallT_p0) << endl;

	// Write the Intermediate Soln
	writeFinalConnectedComponentToFile(out_name2);

	return 0;
	
	// Add general case functions
	cout<< "ReInitialize functions"<< endl;
	clock_t currTS_c0	= clock();
	double currWallT_c0 = getWallTime();

	setUnionFindToCompleteClustering();
	extractRecordsInClusters();
	cout<< "Records Extracted" << endl;

	getBlockingIDArrayForRecordsInLargeClusters();
	getBlockingIDArrayForRecordsInSmallClusters();
	cout<< "Got blockingID array" << endl;

	sortLargeClusterRecordsBlockingIDArray();
	sortSmallClusterRecordsBlockingIDArray();
	cout<< "Sorted BlockingID Array" << endl;

	removeRedundentLargeClusterRecordsBlockingID();
	removeRedundentSmallClusterRecordsBlockingID();
	cout<< "Removed redundent blockID, recID pairs" << endl;

	findBlockBoundariesForLargeClusterRecordsBlockList();
	findBlockBoundariesForSmallClusterRecordsBlockList();
	cout<< "Found Block Boundaries" << endl;

	clock_t currTS_c1	= clock();
	double doBlocking_t	= (double)(clock() - currTS_c0) / CLOCKS_PER_SEC;
	double currWallT_c1 = getWallTime();
	cout<< "GenRLA: Blocking Processor Time "<< doBlocking_t << endl;
	cout<< "GenRLA: Blocking Wall Time "<< (double)(currWallT_c1 - currWallT_c0) << endl;

	// DO PGEN-RLA HERE
	// Find Block assignments
	
	// sortByBlockSizes();
	// findBlockAssignments();
	// double loadBalancing_p4_t	= (double)(clock() - currTS_p4) / CLOCKS_PER_SEC;
    // cout<< "Get Load Balancing Time "<< loadBalancing_p4_t << endl;


	getEdgesFromBlockedRecordsInSmallClusters();
	cout<< "Found edges from records in small clusters " << endl;
	getEdgesFromBlockedRecordsInSmallClustersToLargeClusters();
	cout<< "Found the new edges" << endl;

	findConnCompOnGeneralizedClustering();
	cout<< "Found the new approximate clusters" << endl;
	findFinalConnectedCompOnGeneralizedClustering();
	cout<< "Found the final clusters" << endl;

	// Total Time Required
    double total_ct	= (double)(clock() - currTS_c0) / CLOCKS_PER_SEC;
	double total_t = total_ct + total_SB_t;
	cout<< "Total processor run time "<< total_t << endl;
	double allDone_pX_Wt = getWallTime();
	cout<< "Get Total Wall Time "<< (double)(allDone_pX_Wt - currWallT_c0 + SB_done_pX_Wt) << endl;

	
	// writeApproximateConnectedComponentToFile(out_name1);

	writeFinalConnectedComponentOnClusteredToFile(out_name2);

	string stat_file_path = out_file_path + stat_file_name;
    ofstream stat_file;
	stat_file.open(stat_file_path);
	stat_file << "DataSize: "<< vec2D.size() << endl;
	// stat_file << "Number of Possible comparison: " << tot_possible_com << endl;
	// stat_file << "Number of pairs compared: " << total_comp  << endl;
	// stat_file << "Reduction Ratio:" << ((long double)total_comp / (long double) tot_possible_com) << endl;
	stat_file << "Number of Edges: "<< totalUniqueRecords-uf[numThreads-1].getSetCount() << endl;
	stat_file << "Total Single Clusters: " << approxConnectedComponents.size()<< endl;
	stat_file << "Total Complete Clusters " << finalConnectedComponents.size() << endl;
	stat_file << "Total Processor Time taken: " << total_t << " Seconds" << endl;
	stat_file << "Total Wall Time Taken: " << (double)(allDone_pX_Wt - currWallT_c0 + SB_done_pX_Wt) << " Seconds" << endl;
	stat_file.close();

    return 0;
}