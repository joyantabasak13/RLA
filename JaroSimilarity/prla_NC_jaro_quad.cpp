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

double distanceThreshold = 0.80;
vector<double> attrDistThreshold{0.80,0.80,0.80,0.80};
int totalRecords;
int lenMax;
int totalUniqueRecords;
int totalUniqueBlocks;
int totalBlockedKmers;
int attributes;
int base = 36;
int kmer = 3;
int blockIDRange = pow(base,kmer+1);
int extraEdges = 0;
int numThreads = 6;
int numSources = 5;
int prelen = 1;
long long int totalCompRequired;

vector<vector<int> > fullExactMatches;
vector<vector<int> > candidateExactmatches;
map<int, vector<int> > approxConnectedComponents;
map<int, vector<int> > exactMatches;
vector<string> vec1D;
vector<vector<string> > vec2D;
vector<vector<int> > clusterExactIndArr;
vector<int> uniqueRecords;
vector<pair<int,string> > headlessCopies;
vector<pair<int, string> > combinedData;
vector<pair<int, string> > noConflictCombinedData;
vector<pair<int,int> > blockingIDList;
vector<pair<int, int>> boundaryArr;
vector<vector<pair<int, int>>> assignedBlocklists;
vector<vector<int> > edgeArr;
vector<int> conflictRecords;
vector<int> noConflictRecords;

int** hash_s1;
int** hash_s2;
int stringMaxLen = 0;

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
		// cout<< "called " << parentInd.size()  <<endl;
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
		// cout<< "called" << endl;
		if (find(recId1) == find(recId2)) {
			// cout<< "returned true" << endl;
			return true;
		} else {
			// cout<< "returned false" << endl;
			return false;
		}
		// cout<< "returning" << endl;
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

	int getNumset() {
		return this->numSets;
	}

	vector<int> getParentArr() {
		return this->parentInd;
	}

	void copyUF (UnionFind uf_c) {
		this->numSets = uf_c.getNumset();
		vector<int> parentIndArr = uf_c.getParentArr();
		for (int i=0; i<parentIndArr.size(); i++) {
			this->parentInd[i] = parentIndArr[i];
		}
	}
};

vector<UnionFind> uf;

void sortCharIDArray(pair<int, int> charIDList[], int length)
{
	int alphabets = 36;
	int numRecords = length;
	pair<int, int> tempArr[numRecords];
	int countArr[alphabets];
	memset(countArr, 0, alphabets*sizeof(int));

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

double jaroDistance_quad(string &s1, string &s2, int id)
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

	// Traverse through the first string
	for (int i = 0; i < len1; i++)
	{

		// Check if there is any matches
		for (int j = max(0, i - max_dist);
			 j < min(len2, i + max_dist + 1); j++)

			// If there is a match
			if (s1[i] == s2[j] && hash_s2[id][j] == 0)
			{
				hash_s1[id][i] = 1;
				hash_s2[id][j] = 1;
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
		if (hash_s1[id][i])
		{

			// Find the next matched character
			// in second string
			while (hash_s2[id][point] == 0)
				point++;
			hash_s2[id][point] = 0;
			if (s1[i] != s2[point++])
				t++;
			
			hash_s1[id][i] = 0;
		}

	t /= 2.0;

	// Return the Jaro Similarity
	return (((double)match) / ((double)len1) + ((double)match) / ((double)len2) + ((double)match - t) / ((double)match)) / 3.0;
}

bool isLinkageOk(vector<string> &a, vector<string> &b, int id)
{
	double cumulativeDist = 0.0;
	for (int i = 1; i < attributes-1; i++)
	{
		double dist = jaroDistance_quad(a[i], b[i], id);
		if(dist < attrDistThreshold[i]) {
			return false;
		} else {
			cumulativeDist += dist;
		}
	}
	cumulativeDist = cumulativeDist/((double)attrDistThreshold.size());
	if (cumulativeDist < distanceThreshold){
		return false;
	}
	return true;
}

double getWallTime() {
	struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
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
			if(attrStr.size() > stringMaxLen) {
				stringMaxLen = attrStr.size();
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
    attributes = vec2D[0].size();
    cout<< "Attributes: "<<attributes << endl;
}

void writeSingleLinkageFullOutput(string& recID_file, string& recInd_file) {
	ofstream out_file;
	ofstream index_file;
    out_file.open(recID_file);
	index_file.open(recInd_file);
	for (auto const& p : approxConnectedComponents) {
        for (int i=0; i<p.second.size(); i++) {
			out_file<< vec2D[p.second[i]][0] << ",";
			index_file<< p.second[i] << ",";
        }
        out_file<< "\n";
		index_file<< "\n";
	}
	out_file.close();
}

void writeExactClusters(string& fileName) {
	ofstream out_file;
    out_file.open(fileName);
	for (auto const& p : exactMatches) {
		for (int j=0; j<p.second.size(); j++) {
			out_file<< p.second[j] << ",";
		}
		out_file<< "\n";
	}
	out_file.close();
}

// Deduplication and Data preprocessing

void getCombinedData() {
	string strSample(50, 'a');
	combinedData.resize(totalRecords);
	int max = 0;
	for (int i = 0; i< vec2D.size(); i++) {
		pair<int, string> p;
		p.first = i;
		for(int j = 1; j<attributes-1; j++ ) {
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
		int lenDiff		= lenMax - combinedData[i].second.size();
		if(lenDiff > 0)
			combinedData[i].second	+= strSample.substr(0, lenDiff);
	}
}

// Do exact clustering from lexically sorted vector of nonconflicting <int,string> pair
void getExactMatches() {
	vector<int> tempVec;

	tempVec.push_back(combinedData[0].first);

	for (int i = 1; i < combinedData.size(); ++i) {
		if(combinedData[i].second.compare(combinedData[i - 1].second) == 0)
			tempVec.push_back(combinedData[i].first);
		else {
			if(!exactMatches.count(tempVec[0])) {
				exactMatches[tempVec[0]] = tempVec;
			} else {
				cout<< "ERROR: Element already EXISTS !! " << endl;
			}
			tempVec.clear();
			tempVec.push_back(combinedData[i].first);
		}
	}
	if(!exactMatches.count(tempVec[0])) {
		exactMatches[tempVec[0]] = tempVec;
	}
	totalUniqueRecords = exactMatches.size();
	cout << "total exact clusters: " << totalUniqueRecords << endl;
}

void getUniqueEntries() {
	uniqueRecords.resize(totalUniqueRecords);
	int i=0;
	for (auto const& p : exactMatches)
    {
        uniqueRecords[i] = p.first;
		i++;
    }
}

// Blocking Functions

string getBlockingString(int recInd) {
	string blockingStr;
	for(int j=0;j<prelen;j++) {	
		for(int i=1; i<attributes-1; i++) {
			if(vec1D[(recInd*attributes) + i].size()){
				blockingStr+= vec1D[(recInd*attributes) + i][0];
			}
		}
	}

	while(blockingStr.size() < kmer) {
		blockingStr += "a";
	}
	rotate(blockingStr.begin(), blockingStr.begin() + 1, blockingStr.end());
	return blockingStr;
}

void getBlockingIDArray() {
	int perAlphaBlocks = pow(base,kmer);
	int alphabets = 26;
	int blockID = 0;
	int indATUnique = 0;
	string blockingStr;
	for (int i = 0; i < totalUniqueRecords; i++) {
		indATUnique = uniqueRecords[i];
		blockingStr = getBlockingString(uniqueRecords[i]);
		int blkstrSize = blockingStr.size();
		for (int j = 0; j < blkstrSize; ++j) {
			blockID = 0;
			for (int k = 0; k < kmer; ++k)
			{
				int blockStrCharCode = (int)blockingStr[(j+k)%blkstrSize];
				if ((blockStrCharCode >= 97) && (blockStrCharCode <= 122)) {
					blockID += ((int)blockingStr[(j+k)%blkstrSize] - 97) * pow(base,k);
				} else if ((blockStrCharCode >= 48) && (blockStrCharCode <= 57)){
					blockID += ((int)blockingStr[(j+k)%blkstrSize] - 48 + alphabets) * pow(base,k);
				}
			}
			if ((blockingStr[0] >= 97) && (blockingStr[0] <= 122)) {
					blockID = ((int)blockingStr[0] - 97)*perAlphaBlocks + blockID;
			} else if ((blockingStr[0] >= 48) && (blockingStr[0] <= 57)){
					blockID = ((int)blockingStr[0] - 48 + alphabets)*perAlphaBlocks + blockID;
			} else {
				blockID = 0;
			}

			pair <int, int> p;
			p.first = blockID;
			p.second = indATUnique;
			blockingIDList.push_back(p);
		}
	}
	cout<< "Blocking ID Pairs: " << blockingIDList.size() << endl;

	// cout<< "Printing bloking ID list"<< endl;
	// cout<< endl;
	// for(int i=0; i<blockingIDList.size(); i++) {
	// 	cout<< "BlockID: " << blockingIDList[i].first << " RecID: " << blockingIDList[i].second << endl;
	// }
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
	// cout<< endl;
	// cout<< "Printing UNIQUE blocking ID list"<< endl;
	// cout<< endl;
	// for(int i=0; i<blockingIDList.size(); i++) {
	// 	cout<< "BlockID: " << blockingIDList[i].first << " RecID: " << blockingIDList[i].second << endl;
	// }
}

void doSortedComp() {
	headlessCopies.resize(2*totalUniqueRecords);
	int blockingStrMaxSize = prelen*(attributes-2);
	cout<< "Max Blocking string size: " << blockingStrMaxSize << endl;
	string strSample(blockingStrMaxSize, '0');

	int main_tid = numThreads - 1 ;
	for(int i=0; i< uniqueRecords.size(); i++) { 
		string blockStr = getBlockingString(uniqueRecords[i]);
		// Pad blockstring to make uniform size
		int lenDiff	= blockingStrMaxSize - blockStr.length();
		if(lenDiff > 0) {
			blockStr += strSample.substr(0, lenDiff);
		}
		// if (i<10){
		// 	cout<< "i: " << i << " BlockString: " << blockStr << " Headless: " << (blockStr.substr(1, blockStr.length()-1) + '0') << endl;
		// }
		headlessCopies[i].first = uniqueRecords[i];
		headlessCopies[i].second = blockStr;
		headlessCopies[totalUniqueRecords+i].first = uniqueRecords[i];
		headlessCopies[totalUniqueRecords+i].second = blockStr.substr(1, blockStr.length()-1) + '0';
	}
	// cout<< "No SegFault in Copying Headless Copies" << endl;

	lenMax = blockingStrMaxSize;
	radixSort(headlessCopies);

	// cout<< "No SegFault in radixSorting Headless Copies" << endl;

	for (int i = 1; i < headlessCopies.size(); i++) {
		if (headlessCopies[i-1].second.compare(headlessCopies[i].second) == 0) {
			int recID_i = headlessCopies[i-1].first;
            int recID_j = headlessCopies[i].first;
			// cout<< "RecID I: " << recID_i << " RecID j: " << recID_j << endl;
            if (!uf[main_tid].isConnected(recID_i, recID_j)) {
				if (isLinkageOk(vec2D[recID_i], vec2D[recID_j], main_tid)) {
                	uf[main_tid].weightedUnion(recID_i, recID_j);
					extraEdges++;
				}
			}
		}
	}
	// cout<< "No SegFault in edge Adding" << endl;
	cout<< "Sorting Edges Added: "<< extraEdges << endl;
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
	cout<< "Total Comp Required: "<< totalCompRequired << endl;
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
			long long int curBlockSize = (long long int) ceil((pow(boundaryArr[j].second,2)-boundaryArr[j].second)/2) ;
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
    for (int i = 0; i < totalRecords; i++)
	{
		root = uf[main_tid].find(i);
		approxConnectedComponents[root].push_back(i);
	}
	
    cout<< "Single Linkage Connected Components: " << approxConnectedComponents.size()<<endl;
}


// Threading related Functions

void mergeEdges() {
	int mainTid = numThreads-1;
	int edgesFromExactClusters = 0;
	for(auto const& p : exactMatches) {
		for (int i = 0; i < p.second.size(); i++)
		{
			int recID_i = p.first;
			int recID_j = p.second[i];
			uf[mainTid].weightedUnion(recID_i, recID_j);
			edgesFromExactClusters++;
		}
	}

	cout<< "Edges added for exact clustering: " << edgesFromExactClusters << endl;

	for(int i= 0; i<numThreads-1; i++ ) {
		for(int j=0; j< totalRecords; j++) {
			int root_i = uf[i].find(j);
			int root_mt = uf[mainTid].find(j);
			if( root_i != root_mt){
				uf[mainTid].weightedUnion(root_i, root_mt);
			}
		}
	}

	// for(int i=0; i<fullExactMatches.size(); i++) {
	// 	int recID_i = fullExactMatches[i][0];
	// 	for(int j = 0; j<fullExactMatches[i].size(); j++) {
	// 		int recID_j = fullExactMatches[i][j];
	// 		uf[mainTid].weightedUnion(recID_i, recID_j);
	// 	}
	// }
}

void getBlockRecords(pair<int, int> &blockInfo, vector<pair<int, vector<string>>> &recordList) {
	recordList.resize(blockInfo.second);
	for (int i = 0; i < blockInfo.second; i++) {
		recordList[i].first = blockingIDList[blockInfo.first+i].second;
		recordList[i].second = vec2D[recordList[i].first];
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
				if (isLinkageOk(blockRecords[i].second, blockRecords[j].second, id)) {
                    uf[id].weightedUnion(recID_i, recID_j);
                }
			}
        }
	}
}

void doNormalThreadedBlocking(int tID) {
	cout<< "Thread: "<< tID << " was assigned " << assignedBlocklists[tID].size() << " blocks" << endl;
	for(int i = 0; i< assignedBlocklists[tID].size(); i++) {
		pair<int, int> block = assignedBlocklists[tID][i];
		vector<pair<int, vector<string>>> blockRecords;
		getBlockRecords(block, blockRecords);
		getEdgesFromBlockedRecords(tID, blockRecords);
	}
	cout << "Thread "<< tID << " Total Edges: "<< totalRecords - uf[tID].getSetCount() << endl;
}

// main function for threads
void *threadDriver(void* ptr) {
	int threadID = *static_cast<int*>(ptr);
	doNormalThreadedBlocking(threadID);
	return 0;
}

int main(int argc, char** argv) {

	// IO PATHS
	// Input

    // string filePath = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/";
    // string out_file_path = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/data/";

	string filePath = "/home/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/";
    string out_file_path = "/home/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/data/";

	string fileName = argv[1];
    filePath = filePath + argv[1];
	getFormattedDataFromCSV(filePath);
	totalRecords = vec2D.size();
	getCombinedData();

	//Initialize
	hash_s1 = new int*[numThreads];
	hash_s2 = new int*[numThreads];
	for(int i = 0; i < numThreads; ++i) {
		hash_s1[i] = new int[stringMaxLen];
		hash_s2[i] = new int[stringMaxLen];
		memset(hash_s1[i], 0, stringMaxLen * sizeof(int));
		memset(hash_s2[i], 0, stringMaxLen * sizeof(int));
	}

	// Outputs
    string fileNameSuffix = "quad_jaro_";
	string out_name1 = out_file_path + "CompleteLinkage_"+ fileName + fileNameSuffix;
	string out_name2 = out_file_path + "Unique_SingleLinkage_"+ fileName + fileNameSuffix;
	string out_name3 = out_file_path + "Unique_SingleLinkage_RecInd_"+ fileName + fileNameSuffix;
	string out_name4 = out_file_path + "ExactClustering_"+ fileName + fileNameSuffix;
	string out_name5 = out_file_path + "ALL_SingleLinkage_"+ fileName + fileNameSuffix;
	string out_name6 = out_file_path + "ALL_RECID_SingleLinkage_" + fileNameSuffix;

	string stat_file_name = "stat_"+ fileName + fileNameSuffix;

	// Sort the Combined Data
	clock_t currTS_p0	= clock();
	double currWallT_p0 = getWallTime();
    radixSort(combinedData);
    double sorting_p0_t	= (double)(clock() - currTS_p0) / CLOCKS_PER_SEC;
    cout<< "Sorting time "<< sorting_p0_t << endl;

	// Get Unique Records
	clock_t currTS_p1	= clock();
	radixSort(noConflictCombinedData);
	getExactMatches();
	// getFullClustersFromExactMatchs();
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
	for (int i=0; i<=numThreads-1; i++) {
		uf[i].setVariable(totalRecords);
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
	mergeEdges();
	doSortedComp();
	double edgeMergingDone_p6_t	= (double)(clock() - currTS_p6) / CLOCKS_PER_SEC;
	cout<< "Get edge Merging Time "<< edgeMergingDone_p6_t << endl;
	cout<< "Total Edges: " << totalRecords-uf[numThreads-1].getSetCount() << endl;

	// Find Connected components (Single Linkage)
	clock_t currTS_p7	= clock();
    findConnComp();
    double findComp_p7_t	= (double)(clock() - currTS_p7) / CLOCKS_PER_SEC;
    cout<< "Connected Comp Find Time "<< findComp_p7_t << endl;

	// Write The clusters to Files
	// writeCompleteClustersFromDeduplication(out_name1);
	// writeUniqueRecordConnectedComponents(out_name2, out_name3);
	writeExactClusters(out_name4);
	writeSingleLinkageFullOutput(out_name5, out_name6);


    double total_SB_t	= (double)(clock() - currTS_p0) / CLOCKS_PER_SEC;
	cout<< "Super Blocking: Total processor run time "<< total_SB_t << endl;
	double SB_done_pX_Wt = getWallTime();
	cout<< "SuperBlocking: Get Total Wall Time "<< (double)(SB_done_pX_Wt - currWallT_p0) << endl;

	string stat_file_path = out_file_path + stat_file_name;
    ofstream stat_file;
	stat_file.open(stat_file_path);
	stat_file << "DataSize: "<< vec2D.size() << endl;
	stat_file << "Number of Edges: "<< totalRecords-uf[numThreads-1].getSetCount() << endl;
	stat_file << "Total Single Clusters: " << approxConnectedComponents.size()<< endl;
	stat_file << "Total Complete Clusters (Deduplication) " << approxConnectedComponents.size() + conflictRecords.size() + fullExactMatches.size() << endl;
	stat_file << "Total Processor Time taken: " << total_SB_t << " Seconds" << endl;
	stat_file << "Total Wall Time Taken: " << (double)(SB_done_pX_Wt) << " Seconds" << endl;
	stat_file.close();

    return 0;
}