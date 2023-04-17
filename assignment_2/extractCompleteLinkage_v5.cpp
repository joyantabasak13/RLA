// Single Linkage to complete linkage
//	try randomized selection of minimum edges
// Currently supports NC voter dataset
// try k NN global random


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

// #include <bits/stdc++.h>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

using namespace std;

int threshold = 0;
int totalRecords;
double avgDistanceThreshold = 0.8;
vector<double> attrDistThreshold{0.8,0.8,0.8,0.8};
int lenMax;
int attributes;
int numSources = 5;
int knn = 100;

vector<vector<string> > vec2D;
vector<vector<int> > cluster2D;
map<int, int> vecRecToClusterRec;
map<int, int> clusterRecToVecRec;
vector<vector<short>> distMat;
map<int, vector<pair<int,int>>> knnMap;
vector<priority_queue<pair<double, pair<int,int> > > > pq;



int getDistance(vector<string>& vec1, vector<string>& vec2);

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

class CompleteCluster {
	public:
		map<int, int> sourceToRec;
		map<int, int> recordToSource;
		int clusterID;
		int recCount = 0;

		bool addRecord(int recID, int sourceID) {
			if (sourceToRec.count(sourceID) == 0) {
				if(isCompatibleLink(recID)) {
					sourceToRec[sourceID] = recID;
					recordToSource[recID] = sourceID;
					recCount++;
					return true;
				}
			}
			return false;
		}

		void setCluster(int ind) {
			clusterID = ind;
		}

		bool isCompatibleLink(int recID) {
			bool isCompatible = true;
			for (auto const& p : recordToSource){
				int dist = getDistance(vec2D[recID], vec2D[p.first]);
				if (dist<(avgDistanceThreshold*1000)) {
					return false;
				}
			}
			return true;
		}

		double getAverageRecordDistance(int recID) {
			int totDist = 0.0;;
			for (auto const& p : recordToSource){
				int dist = getDistance(vec2D[recID], vec2D[p.first]);
				if (dist<(avgDistanceThreshold*1000)) {
					return ((double)threshold);
				} else {
					totDist+=((double)dist);
				}
			}
			return (double)(((double)totDist)/((double)recordToSource.size()));
		}

		bool mergeCompleteClusters(CompleteCluster& cluster) {
			for (auto const& p : cluster.sourceToRec){
				if(this->sourceToRec.count(p.first) != 0) {
					return false;
				}
				if(isCompatibleLink(p.second)== false) {
					return false;
				}
			}
			for (auto const& p : cluster.sourceToRec){
				this->sourceToRec[p.first] = p.second;
				this->recordToSource[p.second] = p.first;
				recCount++;
			}
			return true;
		}

		int getRecCount() {
			return recCount;
		}

		double getClusterDistance(CompleteCluster& cluster){
			double totDist = 0.0;
			int comp = 0;
			for (auto const& p : cluster.recordToSource){
				for (auto const& q : this->recordToSource){
					int dist = getDistance(vec2D[p.first], vec2D[q.first]);
					totDist += ((double)dist);
					comp++;
				}
			}
			return (double)(((double)totDist)/((double)(comp)));
		}
};

vector<CompleteCluster> completeClusters;
vector<CompleteCluster> minCompleteClusters;
vector<CompleteCluster> maxCompleteClusters;

class Cluster {
	public:
		set<int> sources;
		set<int> recIDs;
		bool isCompatible(const set<int>& srcs) {
			for(auto src : srcs) {
				if(sources.count(src)>0) {
					return false;
				}
			}
			return true;
		}

		bool mergeCluster(const Cluster& cluster) {
			bool compatible = isCompatible(cluster.sources);
			if (compatible) {
				for(auto src : cluster.sources) {
					this->sources.insert(src);
				}
				for(auto recID : cluster.recIDs) {
					this->recIDs.insert(recID);
				}
				return true;
			}
			return false;
		}

		double getDistance(const Cluster& cluster){
			if(!isCompatible(cluster.sources)) {
				return ((double)threshold);
			}
			// double totalDist = 0.0;
			//double minDist = ((double)threshold);
			double maxDist = 1000.0;
			for(auto recID1 : cluster.recIDs) {
				for(auto recID2 : this->recIDs) {
					//totalDist += ((double)distMat[recID1][recID2]);
					double dist = ((double)distMat[recID1][recID2]);
					if (dist < avgDistanceThreshold*1000) {
						return ((double)threshold);
					} else if (maxDist > dist) {
						maxDist = dist;
					}
					// if(totalDist>threshold) {
					// 	return ((double)threshold);
					// }
				}
			}
			// double avgDist = (double)((double)totalDist)/((double)(cluster.recIDs.size()*this->recIDs.size()));
			// return avgDist;
			return maxDist;
		}
};

map<int, Cluster> clusterIdToCluster;

int jaro_distance(string& s1, string& s2)
{
	// Length of two strings
	int len1 = s1.length();
	int len2 = s2.length();

	if (len1 == 0 || len2 == 0)
		return 0.0;

	// Maximum distance upto which matching is allowed
	int max_dist = floor(max(len1, len2) / 2) - 1;

	// Count of matches
	int match = 0;

	// Hash for matches
	std::vector<int> hash_s1(s1.length(),0);
	std::vector<int> hash_s2(s2.length(),0);

	// Traverse through the first string
	for (int i = 0; i < len1; i++) {

		// Check if there is any matches
		for (int j = max(0, i - max_dist);
			j < min(len2, i + max_dist + 1); j++)
			// If there is a match
			if (s1[i] == s2[j] && hash_s2[j] == 0) {
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
		if (hash_s1[i]) {

			// Find the next matched character
			// in second string
			while (hash_s2[point] == 0)
				point++;

			if (s1[i] != s2[point++])
				t++;
		}

	t /= 2;

	// Return the Jaro Similarity
	double jaroSimVal = (((double)match) / ((double)len1)
			+ ((double)match) / ((double)len2)
			+ ((double)match - t) / ((double)match))
		/ 3.0;
	
	jaroSimVal = jaroSimVal* 100.0;
	return (int)jaroSimVal;
}

// I/O Functions
void writeCompleteClusters(string& recID_file) {
	ofstream out_file;
    out_file.open(recID_file);
	for (int i = 0; i< completeClusters.size(); i++) {
		for (auto const& p : completeClusters[i].recordToSource){
			out_file<< vec2D[p.first][0] << ",";
			// out_file<< p.first << ",";
        }
        out_file<< "\n";
	}
	out_file.close();
}

void getData(string& filePath) {
    string line;
    ifstream records(filePath);

    int ind = 0;
    while (getline (records, line)) {
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
    cout<< "Records read: " << vec2D.size() << endl;
	cout<< "Attributes: "<<attributes << endl;
}

int getDistance(vector<string>& vec1, vector<string>& vec2){
	int dist = 0;
	double avgDist = 0.0;
	for(int i=1; i<attributes-1; i++){
		dist = jaro_distance(vec1[i], vec2[i]);
		if(dist<attrDistThreshold[i-1]) {
			return threshold;
		}
		avgDist +=((double)dist);
	}
	avgDist = avgDist / ((double)attrDistThreshold.size());
	if(avgDist < avgDistanceThreshold) {
		return threshold;
	}
	return (int)(avgDist);
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

void getSingleLinkageClusters(string& filePath) {
	string line;
    ifstream clusters(filePath);

    int rec = 0;
    while (getline (clusters, line)) {
        vector<string> vec;
		vector<int> cluster;
        boost::split(vec, line, boost::is_any_of(","));
		for(int i = 0; i < vec.size()-1; i++) {
			cluster.push_back(stoi(vec[i]));
			rec++;
		}
        cluster2D.push_back(cluster);
    }
    clusters.close();
	cout<< "Total Single linkage Clusters: " << cluster2D.size() << endl;
	cout<< "Total records found in Clusters: " << rec << endl;
}

void getCombinedData(vector<pair<int, string> >& combinedData, int sl_clusterID) {
	string strSample(50, 'a');
	int clusterSize = cluster2D[sl_clusterID].size();
	combinedData.resize(clusterSize);
	int max = 0;
	for (int i = 0; i < clusterSize; i++) {
		pair<int, string> p;
		p.first = cluster2D[sl_clusterID][i];
		for(int j = 1; j<attributes; j++ ) {
			p.second = p.second + vec2D[p.first][j];
		}
		combinedData[i]=p;
		if (max<p.second.size()) {
			max = p.second.size();
		}
	}
	lenMax = max;

    for (int i = 0; i < clusterSize; ++i) {
		int lenDiff		= lenMax - combinedData[i].second.size();
		if(lenDiff > 0)
			combinedData[i].second	+= strSample.substr(0, lenDiff);
	}
}

void extractRecordConflictStatus(vector<pair<int, string> >& combinedData, vector<int>& conflictRecords, vector<int>& noConflictRecords) {
	vector<bool> conflictStatus;
	conflictStatus.resize(combinedData.size(), false);
	for (int i = 0; i < combinedData.size()-1; i++) {
		if(combinedData[i].second.compare(combinedData[i + 1].second) == 0) {
			conflictStatus[i] = true;
			conflictStatus[i+1] = true;
		}
	}
	for (int i = 0; i< combinedData.size(); i++) {
		if(conflictStatus[i] == true) {
			conflictRecords.push_back(combinedData[i].first);
		} else {
			noConflictRecords.push_back(combinedData[i].first);
		}
	}
// 	cout << "total Conflicted Records: " << conflictRecords.size() << endl;
// 	cout << "total Non-Conflicted Records: " << noConflictRecords.size() << endl;
}

void makeConflictingRecordsSingletonClusters(vector<int>& conflictRecords){
	for(int i = 0; i<conflictRecords.size(); i++) {
		CompleteCluster completeCluster;
		int recID = conflictRecords[i];
		int sourceID = stoi(vec2D[recID][attributes-1]);
		bool isAdded = completeCluster.addRecord(recID,sourceID);
		if(isAdded) {
			completeClusters.push_back(completeCluster);
		} else {
			cout<< "Error: in making conflicting records Singleton clusters" << endl;
		}
	}
}

void getNoConflictCombinedData(vector<pair<int, string> >& noConflictCombinedData, vector<int>& noConflictRecords) {
	string strSample(50, '0');
	noConflictCombinedData.resize(noConflictRecords.size());
	int max = 0;
	for (int i = 0; i< noConflictRecords.size(); i++) {
		pair<int, string> p;
		p.first = noConflictRecords[i];
		for(int j = 1; j<attributes-1; j++ ) {
			p.second = p.second + vec2D[noConflictRecords[i]][j];
		}
		noConflictCombinedData[i]=p;
		if (max<p.second.size()) {
			max = p.second.size();
		}
	}
	lenMax = max;
    for (int i = 0; i < noConflictCombinedData.size(); ++i) {
		int lenDiff		= lenMax - noConflictCombinedData[i].second.length();
		if(lenDiff > 0)
			noConflictCombinedData[i].second	+= strSample.substr(0, lenDiff);
	}
}

void getExactMatches(vector<pair<int, string> >& noConflictCombinedData, map<int, vector<int> >& exactMatches) {
	vector<int> tempVec;
	vecRecToClusterRec.clear();
	clusterRecToVecRec.clear();
	clusterIdToCluster.clear();
	int count = 0;
	if(noConflictCombinedData.size() == 0) {
		return;
	}

	tempVec.push_back(noConflictCombinedData[0].first);

	for (int i = 1; i < noConflictCombinedData.size(); ++i) {
		if(noConflictCombinedData[i].second.compare(noConflictCombinedData[i - 1].second) == 0)
			tempVec.push_back(noConflictCombinedData[i].first);
		else {
			if(!exactMatches.count(tempVec[0])) {
				if (tempVec.size() == numSources) {
					CompleteCluster cluster;
					for(int k=0; k<numSources; k++) {
						bool isAdded = cluster.addRecord(tempVec[k], stoi(vec2D[tempVec[k]][attributes-1]));
						if(!isAdded) {
							cout<< "ERROR: FullSize cluster extraction error" << endl;
						}
					}
					completeClusters.push_back(cluster);
				} else {
					exactMatches[tempVec[0]] = tempVec;
					vecRecToClusterRec[tempVec[0]] = count;
					clusterRecToVecRec[count] = tempVec[0];
					Cluster cluster;
					cluster.recIDs.insert(count);
					for(int l =0; l<tempVec.size(); l++){
						cluster.sources.insert(stoi(vec2D[tempVec[l]][attributes-1]));
					}
					clusterIdToCluster[count] = cluster;
					count++;

				}	
			} else {
				cout<< "ERROR: Element already EXISTS !! " << endl;
			}
			tempVec.clear();
			tempVec.push_back(noConflictCombinedData[i].first);
		}
	}
	if(!exactMatches.count(tempVec[0])) {
		if (tempVec.size() == numSources) {
			CompleteCluster cluster;
			for(int k=0; k<numSources; k++) {
				bool isAdded = cluster.addRecord(tempVec[k], stoi(vec2D[tempVec[k]][attributes-1]));
				if(!isAdded) {
					cout<< "ERROR: FullSize cluster extraction error" << endl;
				}
			}
			completeClusters.push_back(cluster);
		} else {
			exactMatches[tempVec[0]] = tempVec;
			vecRecToClusterRec[tempVec[0]] = count;
			clusterRecToVecRec[count] = tempVec[0];
			Cluster cluster;
			cluster.recIDs.insert(count);
			for(int l =0; l<tempVec.size(); l++){
				cluster.sources.insert(stoi(vec2D[tempVec[l]][attributes-1]));
			}
			clusterIdToCluster[count] = cluster;
			count++;
		}
	}
	cout << "total exact clusters: " << exactMatches.size() << " Labled " << count << " records" << endl;
}

void getDistMat(){
	// cout<< "Started distMat calc" << endl;
	pq = vector<priority_queue<pair<int, pair<int,int> >>>();

	distMat.clear();
	int n=clusterRecToVecRec.size();
	pq.resize(n);
	vector<vector<string>> dataArr;
	dataArr.resize(n);
	cout<< endl;
	for(int i=0; i<n; i++) {
		dataArr[i] = vec2D[clusterRecToVecRec[i]];
	}
	distMat.resize(n);

	for (int i = 0; i < n; ++i) {
		distMat[i].resize(n, threshold);
	}
	for (int i = 0; i < n; ++i) {
		double minDist = (double)threshold-1.0;
		int minIndex = -1;
		for(int j = i+1; j<n; j++) {
			int dist = getDistance(dataArr[i], dataArr[j]);
			distMat[i][j] = dist; 
			distMat[j][i] = dist;
			// for nearest neighbours keep a pq. Later Extract the mins and choose randomly
			if(dist > avgDistanceThreshold) { 
				minDist = (double)dist;
				minIndex = j;
				pair<double, pair<int,int> > edge;
				edge.first = minDist;
				edge.second.first = i;
				edge.second.second = minIndex;
				pq[i].push(edge);
				// pq[minIndex].push(edge);
			}
		}
	}
	for(int i=0; i<n; i++) {
		vector<pair<double, pair<int,int>>>  kMinEdges;
		// Extract the k mins
		while (!pq[i].empty())
		{
			kMinEdges.push_back(pq[i].top());
			// surpressed k neighboughr logic
			// if(kMinEdges.size() >= knn) {
			// 	break;
			// }
			pq[i].pop();
		}

		for(auto const& edge: kMinEdges){
			knnMap[(int)(edge.first)].push_back(edge.second);
		}
	}
}

void extractCompleteClusters(map<int, vector<int> >& exactMatches) {
	int clusterID = clusterIdToCluster.size();
	int maxKey = 0;
	int edgeTot = 0;
	// shuffle all edges
	for(auto const& p: knnMap ){
		//random shuffling
		auto rng = std::default_random_engine {};
		std::shuffle(std::begin(knnMap[p.first]), std::end(knnMap[p.first]), rng);
		if(p.first > maxKey) {
			maxKey = p.first;
		}
		edgeTot += knnMap[p.first].size();
	}
 
	while(edgeTot > 0) {
		// update maxKey in each iteration if maxKey val edge exists
		if(knnMap[maxKey].size() < 1) {
			maxKey = 0;
			for(auto const& p: knnMap ){
				if (knnMap[p.first].size() > 0) {
					if (p.first> maxKey) {
						maxKey = p.first;
					}
				}
			}
		}
		// Do random selection
		int randMinIndex = rand() % knnMap[maxKey].size();
		pair<int,int> edge;
		edge = knnMap[maxKey][randMinIndex];
		int u = edge.first;
		int v = edge.second;
		

		// remove edge from map
		knnMap[maxKey][randMinIndex] = knnMap[maxKey][knnMap[maxKey].size()-1];
		knnMap[maxKey].pop_back();
		edgeTot--;

		if(clusterIdToCluster.count(u) > 0 && clusterIdToCluster.count(v) > 0) {
			Cluster cluster = clusterIdToCluster[u];
			bool isAdded;
			isAdded = cluster.mergeCluster(clusterIdToCluster[v]);
			if(isAdded) {
				clusterIdToCluster.erase(u);
				clusterIdToCluster.erase(v);
				if(cluster.sources.size() == numSources) {
					// Add to complete clusters
					CompleteCluster com_cluster;
					for(auto recID : cluster.recIDs) {
						for(int k=0; k<exactMatches[clusterRecToVecRec[recID]].size(); k++) {
							int orig_recID = exactMatches[clusterRecToVecRec[recID]][k];
							int orig_source = stoi(vec2D[orig_recID][attributes-1]);
							bool isAddedToComplete = com_cluster.addRecord(orig_recID, orig_source);
							if (!isAddedToComplete) {
								cout<< "ERROR: IN extractCompleteCluster, merge cluster has issues" << endl;
							}
						}
					}
					completeClusters.push_back(com_cluster);
				} else {
					// cout<< "Merged cluster is not Complete!" << endl;
					// Traverse the map and find the nearest neighbour index
					double minDist = ((double)threshold);
					int minIndex = -1;
					priority_queue<pair<double, pair<int,int> > > nnPQ;
					vector<pair<double, pair<int,int> > > edgesToAdd;
					for (const auto &p : clusterIdToCluster){
						double dist = cluster.getDistance(p.second);
						if (dist >= minDist)
						{
							minDist = (double)dist;
							minIndex = p.first;
							pair<double, pair<int,int> > edge;
							edge.first = (double)minDist;
							edge.second.first = clusterID;
							edge.second.second = minIndex;
							nnPQ.push(edge);
						}
					}
					while(!nnPQ.empty()) {
						edgesToAdd.push_back(nnPQ.top());
						if(edgesToAdd.size() >= knn) {
							break;
						}
						nnPQ.pop();
					}
					for(auto const& edge: edgesToAdd){
						int edgeWeight = edge.first;
						if (edgeWeight <= 0) {
							cout<< "ERROR: Edgeweight Zero or Negative" << endl;
						}
						knnMap[edgeWeight].push_back(edge.second);
						edgeTot++;
						if (edgeWeight< minKey) {
							minKey = edgeWeight;
						}
					}

					// cout<< "New merged cluster " << clusterID << " is in Play" << endl;
					clusterIdToCluster[clusterID] = cluster;
					clusterID++;
				}
			}
		}
	}
	// ADD all remaining clusters to complete clusters
	for (const auto &p : clusterIdToCluster) {
		CompleteCluster com_cluster;
		// cout<< "Cluster " << p.first << " is a remaining cluster" << endl;
		for(auto recID : p.second.recIDs) {
			for(int k=0; k<exactMatches[clusterRecToVecRec[recID]].size(); k++) {
				int orig_recID = exactMatches[clusterRecToVecRec[recID]][k];
				int orig_source = stoi(vec2D[orig_recID][attributes-1]);
				// cout<< "Addind rec " << recID << " from source " << orig_source << "to complete cluster" << endl;
				bool isAddedToComplete = com_cluster.addRecord(orig_recID, orig_source);
				if (!isAddedToComplete) {
					cout<< "ERROR: Merging remaining cluster, record: " << orig_recID << " From source " << orig_source << " Has conflict" << endl;
				} else {
					// cout<< "Added Succesfully" << endl;
				}
			}
		}
		completeClusters.push_back(com_cluster);
	}
}

void getAllCompleteClusters(string& filePath) {
	for(int i=0; i<cluster2D.size(); i++) {
	// for(int i=0; i<100; i++) {
		// move trivial size 1 clusters to complete clusters
		if (cluster2D[i].size() == 1)
		{
			CompleteCluster completeCluster;
			int recID = cluster2D[i][0];
			int sourceID = stoi(vec2D[recID][attributes-1]);
			completeCluster.addRecord(recID,sourceID);
			completeClusters.push_back(completeCluster);
			continue;
		}

		// copy data
		cout<< "Processing cluster: " << i << " of size: " << cluster2D[i].size() << endl;
		vector<pair<int, string> > combinedData;
		vector<int> conflictingRecords;
		vector<int> nonConflictingRecords;
		vector<vector<int> > fullExactMatches;
		vector<vector<int> > candidateExactmatches;
		map<int, vector<int> > exactMatches;

		getCombinedData(combinedData, i);
		radixSort(combinedData);
		extractRecordConflictStatus(combinedData, conflictingRecords, nonConflictingRecords);

		// remove same source exact matches and make them complete clusters
		makeConflictingRecordsSingletonClusters(conflictingRecords);

		vector<pair<int, string> > noConflictCombinedData;
		getNoConflictCombinedData(noConflictCombinedData, nonConflictingRecords);
		radixSort(noConflictCombinedData);
		getExactMatches(noConflictCombinedData, exactMatches);
		getDistMat();
		// cout<< "No error in DistMat for cluster:" << i << " Pq size: "<< pq.size() << endl;
		// extract complete clusters
		extractCompleteClusters(exactMatches);
	}
	writeCompleteClusters(filePath);
	cout<< "Total Clusters " << completeClusters.size() << endl;
}


int main(int argc, char** argv) {

	// IO PATHS
	// string filePath = "/home/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/";
    // string singleLinkageFilePath = "/home/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/data/";
	// string completeLinkageFilePath = "/home/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/data/";

	string filePath = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/";
    string singleLinkageFilePath = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/data/";
	string completeLinkageFilePath = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/data/";

	filePath = filePath + argv[1];
	string singleLinkageFileName = argv[2];
	singleLinkageFilePath = singleLinkageFilePath + argv[2];
	string completeLinkageClusterFile = "kNN_Global_Linked_NewImp" + singleLinkageFileName;
	completeLinkageFilePath = completeLinkageFilePath + completeLinkageClusterFile;


	getData(filePath);
	totalRecords = vec2D.size();
	getSingleLinkageClusters(singleLinkageFilePath);

	// Extract Complete linkage clusters from Single linkage clusters
	clock_t currTS_p0	= clock();
    getAllCompleteClusters(completeLinkageFilePath);
    double sorting_p0_t	= (double)(clock() - currTS_p0) / CLOCKS_PER_SEC;
    cout<< "Complete Linkage Time "<< sorting_p0_t << endl;

    return 0;
}