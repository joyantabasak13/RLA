#include <mpi.h>
#include <stdio.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
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

#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

using namespace std;

#define ROOT_ID				0

int threshold = 99;
int blockingDistanceThreshold = 1;
int nonBlockingDistanceThreshold = 1;
int totalNonBlockingAttrThreshold = 1;
int totalRecords;
int lenMax;
int totalUniqueRecords;
int totalUniqueBlocks;
int totalBlockedKmers;
int attributes;
int base = 36;
int kmer = 3;
int blockIDRange = pow(base,kmer);
int blockingAttrIndex = 1;
int extraEdges = 0;
int numThreads = 0; // will be reset to world
int totalEdges = 0;
long long int totalCompRequired;

vector<set<int> > block_list;
vector<vector<int> > exactmatches;
map<int, vector<int> > approxConnectedComponents;
vector<vector<int> > finalConnectedComponents;
vector<string> vec1D;
vector<vector<string> > vec2D;
vector<vector<int> > clusterExactIndArr;
vector<pair<int,string> > uniqueRecords;
vector<pair<int,string> > headlessCopies;
vector<pair<int, string> > combinedData;
vector<pair<int,int> > blockingIDList;
vector<pair<int, int> > boundaryArr;
vector<vector<pair<int, int> > > assignedBlocklists;
vector<vector<int>> assignedRecordIndices;
vector<vector<int> > edgeArr;

int matArr[50][50] = {0};

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

	void setParentArr(vector<int> &parentArr){
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


void radixSort(vector<pair<int, string> > &strDataArr){
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


bool isLinkageOk(vector<string> &a, vector<string> &b, int singleNonBlockingAttrThreshold)
{
	int singleAttributeDist = 0;
	int cumulativeNonBlockingDist = 0;
	for (int i = 1; i < a.size(); i++)
	{
		singleAttributeDist = calculateBasicED(a[i], b[i], singleNonBlockingAttrThreshold);
		if (singleAttributeDist > singleNonBlockingAttrThreshold){
			return true;
		}

		cumulativeNonBlockingDist += singleAttributeDist;
		if (cumulativeNonBlockingDist > totalNonBlockingAttrThreshold){
			return false;
		}
	}    
    
	return false;
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
    while (getline (records, line)) {
        vector<string> result;
        boost::split(result, line, boost::is_any_of(","));
        vector<string> vec;
		for(int i=0; i<result.size(); i++) {
			auto last = std::remove_if(result[i].begin(), result[i].end(), [](auto ch) {
        								return ::ispunct(ch) || ::iswpunct(ch);
    								});
			result[i].erase(last, result[i].end()); //Remove junk left by remove_if() at the end of iterator
			boost::to_lower(result[i]);
			vec.push_back(result[i]);
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


void getCombinedData() {
	string strSample(50, '0');
	combinedData.resize(totalRecords);
	int max = 0;
	for (int i = 0; i< vec2D.size(); i++) {
		pair<int, string> p;
		p.first = i;
		p.second = vec2D[i][1] + vec2D[i][2] + vec2D[i][3];
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


void getExactMatches() {
	vector<int> tempVec;

	tempVec.push_back(0);

	for (int i = 1; i < totalRecords; ++i) {
		if(combinedData[i].second.compare(combinedData[i - 1].second) == 0)
			tempVec.push_back(i);
		else {
			exactmatches.push_back(tempVec);
			tempVec.clear();
			tempVec.push_back(i);
		}
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
    }
}

// Blocking Functions
int getKmerCount() {
	long long int totalKmerCount = 0;
	string blockingStr;
	for (int i = 0; i < totalUniqueRecords; i++) {
		totalKmerCount += vec1D[uniqueRecords[i].first*attributes + 1].size() - kmer + 1;
	}
	return totalKmerCount;
}

void getBlockingIDArray() {
	int totalKMers = getKmerCount();
	blockingIDList.resize(totalKMers);
	int ind = 0;
	int blockID = 0;
	string blockingStr;
	for (int i = 0; i < totalUniqueRecords; i++) {
		blockingStr = vec1D[uniqueRecords[i].first*attributes + blockingAttrIndex];
		string temp_str = vec1D[uniqueRecords[i].first*attributes + blockingAttrIndex];
		while(blockingStr.size() < kmer) {
			if (blockingStr.size() == 0) {
				blockingStr = "a";
				temp_str = "a";
			}
			blockingStr = blockingStr + temp_str;
		}

		int blkstrSize = blockingStr.size();
		for (int j = 0; j < blkstrSize; ++j) {
			blockID = 0;
			for (int k = 0; k < kmer; ++k)
			{
				int blockStrCharCode = (int)blockingStr[(j+k)%blkstrSize];
				if ((blockStrCharCode >= 97) && (blockStrCharCode <= 122)) {
					blockID += (blockStrCharCode - 97) * pow(base,k);
				} else if ((blockStrCharCode >= 48) && (blockStrCharCode <= 57)){
					blockID += (blockStrCharCode - 48 + 26) * pow(base,k);
				}
			}
			blockingIDList[ind].first = blockID;

            //changed here
			// blockingIDList[ind].second = indATUnique;
            blockingIDList[ind].second = uniqueRecords[i].first;

			ind++;
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
			totalCompRequired = totalCompRequired + ceil((range*(range-1))/2);
			curBlockInd++;
			startInd = i;
		}
	}
	// Enter last Block info
	range = numRecords-startInd;
	totalCompRequired = totalCompRequired + ceil((range*(range-1))/2);
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
	numThreads--;

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
            int range = boundaryArr[j].second;
			long long int curBlockSize = ceil((range*(range-1))/2);
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
                    int range = boundaryArr[j].second;
					curAssignmentSize = curAssignmentSize + ceil((range*(range-1))/2);
					startInd = j;
				} else {
					break;
				}
			}
		}
		cout<< "Thread: "<< (i+1) << " was assigned: " << curAssignmentSize << " comparisions where threshold is: " << threshold << endl;
	}
	cout<< "Start ind: " << startInd << " last ind " << lastInd << endl; 
	if (lastInd - startInd > 1) {
		for (int i = startInd + 1; i<lastInd; i++) {
			assignedBlocklists[i%numThreads].push_back(boundaryArr[i]);
		}
	}
	numThreads++;
}


void doAssignedRecordComp(int pid) {
	for(int i = 0; i < assignedRecordIndices.size(); i++) {
		for(int j=0; j < assignedRecordIndices[i].size(); j++) {
			for(int k=j+1; k<assignedRecordIndices[i].size(); k++){
				int recid_j = assignedRecordIndices[i][j];
				int recid_k = assignedRecordIndices[i][k];
				if(!uf.isConnected(recid_j, recid_k)){
				 	if(isLinkageOk(vec2D[recid_j], vec2D[recid_k], blockingDistanceThreshold)) {
				 		uf.weightedUnion(recid_j, recid_k);
				 		totalEdges++;
				 	}
				}
			}
		}
	}
	// cout<< "Processor Rank: " << pid << " found edges: " << totalEdges << endl;
}

int mergeEdges(UnionFind &ufTemp) {
	for(int i=0; i< totalRecords; i++) {
		int root_i = uf.find(i);
		int root_mt = ufTemp.find(i);
		if( root_i != root_mt){
			uf.weightedUnion(root_i, root_mt);
		}
	}
}

void findConnComp()
{
	int root, edgeTotal;

	for(int i=0; i<exactmatches.size(); i++){
		for(int j=1; j<exactmatches[i].size(); j++){
			uf.weightedUnion(combinedData[exactmatches[i][j]].first, combinedData[exactmatches[i][0]].first);
		}
	}
	cout<< "Done merging Exact Matches" << endl;
    for (int i = 0; i < totalRecords; i++)
	{
		root = uf.find(i);
		approxConnectedComponents[root].push_back(i);
	}

    cout<< "#Connected Components: " << approxConnectedComponents.size()<<endl;
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
            dataArr[c] = vec2D[p.second[c]];
        };

		// generate a 2D matrix filled with all pair comparision results
        for (int i =0; i<numComponents; i++) {
            distmat[i][i] = true;
            for (int j = i+1; j < numComponents; j++)
            {
                if (isLinkageOk(dataArr[i], dataArr[j], blockingDistanceThreshold)) {
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
    cout<< "Total Nodes: "<< totalNodes << " unique records: " << totalUniqueRecords << endl;
    cout<< "Total Components: "<< finalConnectedComponents.size()<<endl;
}


int main(int argc, char** argv) {
    MPI_Init(NULL, NULL);
    int rank;
    int world;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world);
    printf("Hello: rank %d, world: %d\n",rank, world);

    numThreads = world;

    string filePath = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/";
    // string filePath = "/home/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/";
    string fileName = argv[1];
    filePath = filePath + argv[1];
    getFormattedDataFromCSV(filePath);
	totalRecords = vec2D.size();
    MPI_Barrier(MPI_COMM_WORLD);
	clock_t currTS_p0	= clock();
    double currWallT_p0 = getWallTime();
    if (rank == ROOT_ID) {
        getCombinedData();
        // Sort the Combined Data
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
        findBlockBoundaries();
        sortByBlockSizes();
        findBlockAssignments();
        double loadBalancing_p4_t	= (double)(clock() - currTS_p4) / CLOCKS_PER_SEC;
        cout<< "Get Load Balancing Time "<< loadBalancing_p4_t << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
	uf.setVariable(totalRecords);

	// TAGS:
	// 0 - assigned number of blocks
	// 1 - individual block sizes
	// 2 - individual block membersList
	// 3 - unionFind arr
    if(rank == ROOT_ID) {
		// for(int i=0; i<numthreads-1; i++){
		// 	MPI_Send(&totalRecords, 1, MPI_INT, i+1, 3, MPI_COMM_WORLD);
		// 	uf.setVariable(totalUniqueRecords);
		// }
        for(int i=0; i<numThreads-1; i++){
			int assignedBlocks = assignedBlocklists[i].size();
            MPI_Send(&assignedBlocks, 1, MPI_INT, i+1, 0, MPI_COMM_WORLD);
			for(int j=0; j<assignedBlocklists[i].size(); j++){
				MPI_Send(&assignedBlocklists[i][j].second, 1, MPI_INT, i+1, 1, MPI_COMM_WORLD);
			}
			for(int j=0; j<assignedBlocklists[i].size(); j++){
				int recordsIndicesForIthBlock[assignedBlocklists[i][j].second];
				for(int k=0; k < assignedBlocklists[i][j].second; k++){
					recordsIndicesForIthBlock[k] = uniqueRecords[assignedBlocklists[i][j].first+k].first;
				}
				MPI_Send(&recordsIndicesForIthBlock, assignedBlocklists[i][j].second, MPI_INT, i+1, 2, MPI_COMM_WORLD);
			}
        }
    } else {
		int assignedBlocksCount = 0;

		MPI_Recv(&assignedBlocksCount, 1, MPI_INT, 0, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        cout<< "Rank: " << rank << " Was assigned " << assignedBlocksCount << " Number of Blocks" << endl;
		assignedRecordIndices.resize(assignedBlocksCount);
		int assignedBlockSizes[assignedBlocksCount];
		for(int j=0; j<assignedBlocksCount; j++){
			int blockSize = 0;
			MPI_Recv(&blockSize, 1, MPI_INT, 0, 1,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			assignedBlockSizes[j] = blockSize;
		}
		int totalRecs = 0;
		for(int j=0; j<assignedBlocksCount; j++){
			// Todo remove redundency 
			int blockRecordIndices[assignedBlockSizes[j]];
			MPI_Recv(&blockRecordIndices, assignedBlockSizes[j], MPI_INT, 0, 2,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			vector<int> temp;
			for(int k=0; k<assignedBlockSizes[j]; k++){
				temp.push_back(blockRecordIndices[k]);
			}
			assignedRecordIndices.push_back(temp);
			totalRecs += assignedBlockSizes[j];
			doAssignedRecordComp(rank);
		}
		cout<< "Rank: " << rank << " received " << totalRecs << " Number of Records and found edges: " << totalEdges << endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if(rank == ROOT_ID) {
		for(int i=1; i<numThreads; i++){
			int tempArr[totalRecords];
            MPI_Recv(&tempArr, totalRecords, MPI_INT, i, 4,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			UnionFind ufTemp;
			ufTemp.setVariable(totalRecords);
			vector<int> tempVec(totalRecords);
			for(int j=0; j<totalRecords; j++){
				tempVec[j] = tempArr[j];
			}
			ufTemp.setParentArr(tempVec);
			mergeEdges(ufTemp);
			cout<< "Merged edges from rank: " << i << endl; 
		}
		cout<< "Merging Done" << endl;
		findConnComp();
		cout<< "finding connected components done" << endl;
		// findFinalConnectedComp();
		// cout<< "finding final connected components done" << endl;
		double final_t	= (double)(clock() - currTS_p0) / CLOCKS_PER_SEC;
        cout<< "Record Linkage time "<< final_t << endl;
	} else {
		int tempVec[totalRecords];
		for(int i=0; i<totalRecords; i++) {
			tempVec[i] = uf.parentInd[i];
		}
		MPI_Send(&tempVec, totalRecords, MPI_INT, 0, 4, MPI_COMM_WORLD);
	} 
	
    MPI_Finalize();
}