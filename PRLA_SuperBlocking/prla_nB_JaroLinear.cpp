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

double distanceThreshold = 0.95;
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
int numThreads = 6;
long long int totalCompRequired;

vector<string> recordVector;
int stringMaxLen = 0;
const int alphabetSize = 256;

int** jwMatS1;
int** jwMatS2;
int** s1CharCount;
int** s2CharCount;
int** s1MatchAtPos;
int** s2MatchAtPos;

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
vector<pair<int, int>> boundaryArr;
vector<vector<pair<int, int>>> assignedBlocklists;
vector<vector<int> > edgeArr;

priority_queue<pair<double, pair<int,int>>> pq;


class UnionFind {
  public:
    int numSets;
	vector<int> parentInd;
	vector<int> setSize;

    UnionFind() { // Constructor with parameters
		numSets = 0;
    }

	void setVariable(int numRecords) {
		this->numSets = numRecords;
      	this->parentInd.resize(numRecords);
		this->setSize.resize(numRecords);
		for (int i = 0; i < numRecords; i++)
		{
			this->parentInd[i] = i;
			this->setSize[i] = -1;
		}
	}

	int find(int recID) {
		int root = recID;
    	while (root != parentInd[root]){
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
		// make larger root point to smaller one
		if (setSize[i] > setSize[j]) { 
			parentInd[i] = j; 
			setSize[j] += setSize[i]; 
		}
		else { 
			parentInd[j] = i; 
			setSize[i] += setSize[j]; 
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
};

vector<UnionFind> uf;

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

int getKmerCount() {
	long long int totalKmerCount = 0;
	string blockingStr;
	for (int i = 0; i < totalUniqueRecords; i++) {
		totalKmerCount += vec1D[uniqueRecords[i].first*attributes + blockingAttrIndex].size();
	}
	return totalKmerCount;
}

double jaroDistanceLinear(string& s1, string& s2, int id){

	int len1 = s1.length();
	int len2 = s2.length();

	if (len1 == 0 || len2 == 0)
		return 0.0;

	int range = floor(max(len1, len2) / 2) - 1;

	int match = 0;
	int ptr1, ptr2;
	int count1, count2;
	
	// // initialize 
	// for(int i=0; i<len1; i++) {
	// 	s1CharCount[id][s1[i]] = 0;
	// }
	// for(int i=0; i<len2; i++) {
	// 	s2CharCount[id][s2[i]] = 0;
	// }

	for(int i=0; i<len1; i++) {
		jwMatS1[id][s1[i] * stringMaxLen + s1CharCount[id][s1[i]]] = i;
		s1CharCount[id][s1[i]]++;
	}

	for(int i=0; i<len2; i++) {
		jwMatS2[id][s2[i] * stringMaxLen + s2CharCount[id][s2[i]]] = i;
		s2CharCount[id][s2[i]]++;
	}

	for(int i=0; i<len1; i++) {
		count1 = s1CharCount[id][s1[i]];
		count2 = s2CharCount[id][s1[i]];

		// Do linear pass to find matchs
		ptr1 = 0;
		ptr2 = 0;
		int offset = s1[i] * stringMaxLen;
		while((ptr1<count1) && (ptr2<count2))
		{
			if (abs(jwMatS1[id][offset + ptr1] - jwMatS2[id][offset + ptr2]) <= range)
			{
				match++;
				s1MatchAtPos[id][jwMatS1[id][offset + ptr1]] = 1;
				s2MatchAtPos[id][jwMatS2[id][offset + ptr2]] = 1;
				ptr1++;
				ptr2++;
			} else {
				if (jwMatS2[id][offset + ptr2] < jwMatS1[id][offset + ptr1])
				{
					ptr2++;
				} else
				{
					ptr1++;
				}
			}
		}
		s1CharCount[id][s1[i]] = 0;
	}

	for(int i=0; i<len2; i++) {
		s2CharCount[id][s2[i]] = 0;
	}

	// Linear pass for transposition calculation
	int point = 0;
	double t = 0.0;
	for (int i = 0; i < len1; i++) {
		if (s1MatchAtPos[id][i])
		{
			while (s2MatchAtPos[id][point] == 0)
				point++;
			s2MatchAtPos[id][point] = 0;
			if (s1[i] != s2[point++])
				t++;
			
			s1MatchAtPos[id][i] = 0;
		}
	}

	t = t / 2.0;

	return (((double)match) / ((double)len1) + ((double)match) / ((double)len2) + ((double)match - t) / ((double)match)) / 3.0;
}


bool isLinkageOk(vector<string> &a, vector<string> &b, double threshold, int id)
{
	double blockField_dist;
	for(int i=1; i<attributes; i++){
		blockField_dist = jaroDistanceLinear(a[i], b[i], id);
		if(blockField_dist < threshold) return false;
	}
	return true;
}

double getLinkageValue(vector<string> &a, vector<string> &b, double threshold){
	double blockField_dist;
	double cumulativeDist;
	for(int i=1; i<attributes; i++){
		blockField_dist = jaroDistanceLinear(a[i], b[i], 0);
		if(blockField_dist < threshold) return 0.0;
		cumulativeDist+=blockField_dist;
	}
	return cumulativeDist / ((double)attributes-1.0) ;
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
			if(result[i].size() > stringMaxLen) {
				stringMaxLen = result[i].size();
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

void writeApproximateConnectedComponentToFile(string& result_file_name) {
	ofstream out_file;
    out_file.open(result_file_name);

	for (auto const& p : approxConnectedComponents) {
        for (int i=0; i<p.second.size(); i++) {
			for (int j=0; j<exactmatches[p.second[i]].size(); j++) {
				out_file<< vec2D[uniqueRecords[p.second[i]].first][0] << ",";
			}
        }
        out_file<< "\n";
	}
	out_file.close();
}

void writeFinalConnectedComponentToFile(string& result_file_name) {
	ofstream out_file;
    out_file.open(result_file_name);

	for(int i = 0; i< finalConnectedComponents.size(); i++) {
		for(int j = 0; j< finalConnectedComponents[i].size(); j++) {
            for(int k=0; k< exactmatches[finalConnectedComponents[i][j]].size(); k++) {
                out_file << vec2D[uniqueRecords[finalConnectedComponents[i][j]].first][0] << ",";
            }
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

// Do exact clustering from lexically sorted vector of int,string pair
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
        // cout<< uniqueRecords[i].first << "\t" << uniqueRecords[i].second << endl;
    }
}

// Blocking Functions
void getBlockingIDArray() {
	int totalKMers = getKmerCount();
	blockingIDList.resize(totalKMers);
	int ind = 0;
	int blockID = 0;
	int indATUnique = 0;
	string blockingStr;
	for (int i = 0; i < totalUniqueRecords; i++) {
		indATUnique = i;
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
			blockingIDList[ind].second = indATUnique;
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
			totalCompRequired = totalCompRequired + pow(range,2);
			curBlockInd++;
			startInd = i;
		}
	}
	// Enter last Block info
	range = numRecords-startInd;
	totalCompRequired = totalCompRequired + pow(range,2);
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
		// cout<< "Thread: "<< i << " was assigned: " << curAssignmentSize << " comparisions where threshold is: " << threshold << endl;
	}
	// cout<< "Start ind: " << startInd << " last ind " << lastInd << endl; 
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
	
    cout<< "#Connected Components: " << approxConnectedComponents.size()<<endl;
}

void findFinalConnectedComp() {
    int totalNodes = 0;
    int pairsAccessed = 0;
    for (auto const& p : approxConnectedComponents) {
        pairsAccessed++;
        int numComponents = p.second.size();
        totalNodes+=numComponents;
		if(numComponents<3) {
			vector<int> connectedComponent;
			for(int i=0; i < numComponents; i++) {
				connectedComponent.push_back(p.second[i]);
			}
			finalConnectedComponents.push_back(connectedComponent);
		} else {
			double distmat[numComponents*numComponents];
			memset(distmat, 0.0, numComponents*numComponents * sizeof(double));
			pq = priority_queue<pair<double, pair<int,int>>>();
			pair<double, pair<int,int> > edge;
			vector<vector<string>> dataArr(numComponents); // to make cache-efficient, keep records in a row
			// Copy Data in cluster
			for(int c=0; c<p.second.size(); c++) {
				dataArr[c] = vec2D[uniqueRecords[p.second[c]].first];
			};

			// generate a 2D matrix filled with all pair comparision results
			for (int i =0; i<numComponents; i++) {
				for (int j = i+1; j < numComponents; j++)
				{
					double similarity = getLinkageValue(dataArr[i], dataArr[j], distanceThreshold);
					if (similarity) {
						edge.first = similarity;
						edge.second.first = i;
						edge.second.second = j;
						pq.push(edge);
						distmat[i*numComponents + j] = similarity;
						distmat[j*numComponents + i] = similarity;
					}
				}
			}

			map<int, vector<int>> clusters;
			int assignedCluster[numComponents];
			for(int i=0; i<numComponents; i++) {
				assignedCluster[i] = i;
				clusters[i].push_back(i);
			}
			int clusterID = numComponents;
			vector<int> cid_first_records;
			vector<int> cid_second_records;
			while(!pq.empty()){
				edge = pq.top();
				int cid_first = assignedCluster[edge.second.first];
				int cid_second = assignedCluster[edge.second.second];
				cid_first_records = clusters[cid_first];
				cid_second_records = clusters[cid_second];
				bool isCompatible = true;
				for(int i=0; i<cid_first_records.size(); i++) {
					for(int j = 0; j<cid_second_records.size(); j++) {
						if(distmat[cid_first_records[i] * numComponents + cid_second_records[j]] == 0.0) {
							isCompatible = false;
							break;
						}
					}
					if(!isCompatible) break;
				}
				if(isCompatible) {
					for(int i=0; i<cid_first_records.size(); i++) {
						assignedCluster[cid_first_records[i]] = clusterID;
						clusters[clusterID].push_back(cid_first_records[i]);
					}
					for(int j = 0; j<cid_second_records.size(); j++) {
						assignedCluster[cid_second_records[j]] = clusterID;
						clusters[clusterID].push_back(cid_second_records[j]);
					}
					clusterID++;
					clusters.erase(cid_first);
					clusters.erase(cid_second);
				}
				pq.pop();
			}

			for (auto const& q : clusters) {
				vector<int> connectedComponent;
				for(int i=0; i<q.second.size(); i++) {
					connectedComponent.push_back(p.second[q.second[i]]);
				}
				finalConnectedComponents.push_back(connectedComponent);
			}
		}
    }
    cout<< "Total Nodes: "<< totalNodes << " unique records: " << totalUniqueRecords << endl;
    cout<< "Total Components: "<< finalConnectedComponents.size()<<endl;
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

void getEdgesFromBlockedRecords(int id, vector<pair<int, vector<string>>> blockRecords) {
	for (int i = 0; i < blockRecords.size() - 1; i++)
	{
        for (int j = i+1; j < blockRecords.size(); j++)
        {
            int recID_i = blockRecords[i].first;
            int recID_j = blockRecords[j].first;
            if (!uf[id].isConnected(recID_i, recID_j)) {
				// cout<< "will check linkage now between: " << blockRecords[i].second[blockingAttrIndex] << " and " << blockRecords[j].second[blockingAttrIndex] << endl;
				if (isLinkageOk(blockRecords[i].second, blockRecords[j].second, distanceThreshold, id)) {
                    uf[id].weightedUnion(recID_i, recID_j);
                }
			}
        }  
	}
}

void doNormalThreadedBlocking(int tID) {
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
	// cout<< "Thread: " << threadID << " Running" << endl;
	doNormalThreadedBlocking(threadID);
	return 0;
}

int main(int argc, char** argv) {
    string filePath = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/";
    // string filePath = "/home/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/";
    string fileName = argv[1];
    filePath = filePath + argv[1];
    
	// Outputs
    string out_file_path = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/data/";
    // string out_file_path = "/home/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/data/";
	string out_name1 = out_file_path + "out_single_linkage_"+ fileName + "_JaroSimilarity_LINEAR_" + to_string(distanceThreshold);
	string out_name2 = out_file_path + "out_complete_linkage_"+ fileName + "_JaroSimilarity_LINEAR_" + to_string(distanceThreshold);
	string stat_file_name = "stat_"+ fileName + "_JaroSimilarity_LINEAR_" + to_string(distanceThreshold);
	
	
	getFormattedDataFromCSV(filePath);
	totalRecords = vec2D.size();
	getCombinedData();

	//Initialize


	jwMatS1 = new int*[numThreads];
	jwMatS2 = new int*[numThreads];
	s1CharCount = new int*[numThreads];
	s2CharCount = new int*[numThreads];
	s1MatchAtPos = new int*[numThreads];
	s2MatchAtPos = new int*[numThreads];

	for(int i = 0; i < numThreads; ++i) {
		jwMatS1[i] = new int[alphabetSize * stringMaxLen];
		jwMatS2[i] = new int[alphabetSize * stringMaxLen];
		s1CharCount[i] = new int[alphabetSize];
		s2CharCount[i] = new int[alphabetSize];
		s1MatchAtPos[i] = new int[stringMaxLen];
		s2MatchAtPos[i] = new int[stringMaxLen];
		memset(s1MatchAtPos[i], 0, stringMaxLen * sizeof(int));
		memset(s2MatchAtPos[i], 0, stringMaxLen * sizeof(int));
		memset(s1CharCount[i], 0, alphabetSize* sizeof(int));
		memset(s2CharCount[i], 0, alphabetSize* sizeof(int));
	}

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
		usleep(10);
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
	double edgeMergingDone_p6_t	= (double)(clock() - currTS_p6) / CLOCKS_PER_SEC;
	cout<< "Get edge Merging Time "<< edgeMergingDone_p6_t << endl;
	cout<< "Total Edges: " << totalUniqueRecords-uf[numThreads-1].getSetCount() << endl;

	// Find Connected components (Single Linkage)
	clock_t currTS_p7	= clock();
    findConnComp();
    double findComp_p7_t	= (double)(clock() - currTS_p7) / CLOCKS_PER_SEC;
    cout<< "Connected Comp Find Time "<< findComp_p7_t << endl;

	// Find Connected components (Complete Linkage)
    clock_t currTS_p8_t	= clock();
    findFinalConnectedComp();
    double findFinalComp_t	= (double)(clock() - currTS_p8_t) / CLOCKS_PER_SEC;
    cout<< "Final Connected Comps Find Time "<< findFinalComp_t << endl;

	// Total Time Required
    double total_t	= (double)(clock() - currTS_p0) / CLOCKS_PER_SEC;
	cout<< "Total processor run time "<< total_t << endl;
	double allDone_pX_Wt = getWallTime();
	cout<< "Get Total Wall Time "<< (double)(allDone_pX_Wt - currWallT_p0) << endl;

	writeApproximateConnectedComponentToFile(out_name1);
	writeFinalConnectedComponentToFile(out_name2);

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
	stat_file << "Total Wall Time Taken: " << (double)(allDone_pX_Wt - currWallT_p0) << " Seconds" << endl;
	stat_file.close();

    return 0;
}