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

int threshold = 99;
int totalRecords;
int cumulativeDistanceThreshold = 5;
vector<int> attrDistThreshold{3,2,3,2};
int lenMax;
int attributes;
int numSources = 5;
int knn = 5;

vector<vector<string> > vec2D;
vector<vector<int> > cluster2D;
map<int, int> vecRecToClusterRec;
map<int, int> clusterRecToVecRec;
vector<vector<short>> distMat;
priority_queue<pair<double, pair<int,int> > > pq;



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
				if (dist>cumulativeDistanceThreshold) {
					return false;
				}
			}
			return true;
		}

		double getAverageRecordDistance(int recID) {
			int totDist = 0.0;;
			for (auto const& p : recordToSource){
				int dist = getDistance(vec2D[recID], vec2D[p.first]);
				if (dist>cumulativeDistanceThreshold) {
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
			double maxDist = 0.0;
			for(auto recID1 : cluster.recIDs) {
				for(auto recID2 : this->recIDs) {
					//totalDist += ((double)distMat[recID1][recID2]);
					double dist = ((double)distMat[recID1][recID2]);
					if (dist >= threshold -1) {
						return ((double)threshold);
					} else if (maxDist < dist) {
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
	for(int i=1; i<attributes-1; i++){
		dist = dist + calculateBasicED(vec1[i], vec2[i], attrDistThreshold[i-1]);
		if(dist>cumulativeDistanceThreshold) {
			return threshold;
		}
	}
	return dist;
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
	pq = priority_queue<pair<double, pair<int,int> > >();
	distMat.clear();
	int n=clusterRecToVecRec.size();
	vector<vector<string>> dataArr;
	dataArr.resize(n);
	cout<< endl;
	for(int i=0; i<n; i++) {
		dataArr[i] = vec2D[clusterRecToVecRec[i]];
		// for(int j = 0; j<attributes;j++) {
		// 	cout<< dataArr[i][j] << "," ;
		// }
		// cout<< endl;
	}
	distMat.resize(n);

	for (int i = 0; i < n; ++i) {
		distMat[i].resize(n, threshold);
	}
	for (int i = 0; i < n-1; ++i) {
		double minDist = (double)threshold-1.0;
		int minIndex = -1;
		priority_queue<pair<double, pair<int,int> > > nnPQ_i;
		for(int j = i+1; j<n; j++) {
			int dist = getDistance(dataArr[i], dataArr[j]);
			distMat[i][j] = dist; 
			distMat[j][i] = dist;
			// for nearest neighbours keep a pq. Later Extract the mins and choose randomly
			if((double)dist <= minDist) { 
				minDist = (double)dist;
				minIndex = j;
				pair<double, pair<int,int> > edge;
				edge.first = (double)((double)minDist)*(-1.0);
				edge.second.first = i;
				edge.second.second = minIndex;
				// cout<< "Added edge u:"<<i<< " v:"<<minIndex << " value: " << edge.first << endl;
				nnPQ_i.push(edge);
			}
		}

		vector<pair<double, pair<int,int>>>  kMinEdges;
		// Extract the k mins
		while (!nnPQ_i.empty())
		{
			kMinEdges.push_back(nnPQ_i.top());
			if(kMinEdges.size() >= knn) {
				break;
			}
			nnPQ_i.pop();
		}
		//random shuffling
		auto rng = std::default_random_engine {};
		std::shuffle(std::begin(kMinEdges), std::end(kMinEdges), rng);

		for(auto edge: kMinEdges){
			pq.push(edge);
		}
	}
}

void extractCompleteClusters(map<int, vector<int> >& exactMatches) {
	int clusterID = clusterIdToCluster.size();
	while(!pq.empty()) {
		vector<pair<double, pair<int,int> > > edgesToAdd;
		// pq has the nearest neighbours
		int u = pq.top().second.first;
		int v = pq.top().second.second;
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
					double minDist = ((double)threshold)-1.00;
					int minIndex = -1;
					priority_queue<pair<double, pair<int,int> > > nnPQ;
					for (const auto &p : clusterIdToCluster){
						double dist = cluster.getDistance(p.second);
						if (dist <= minDist)
						{
							minDist = (double)dist;
							minIndex = p.first;
							pair<double, pair<int,int> > edge;
							edge.first = (double)((double)minDist)*(-1.0);
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
					// cout<< "New merged cluster " << clusterID << " is in Play" << endl;
					clusterIdToCluster[clusterID] = cluster;
					clusterID++;
				}
			}
		}
		pq.pop();
		double topDist = pq.top().first;
		// cout<< "TopDist: " << topDist << endl;
		while(topDist>=pq.top().first && pq.empty() == false) {
			// cout<< "Current top: " << pq.top().first << " PQ Size: " << pq.size() << endl;
			edgesToAdd.push_back(pq.top());
			pq.pop();
		}
		// Shuffle 
		auto rng = std::default_random_engine {};
		std::shuffle(std::begin(edgesToAdd), std::end(edgesToAdd), rng);

		for(auto edge: edgesToAdd) {
			pq.push(edge);
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
	string filePath = "/home/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/";
    string singleLinkageFilePath = "/home/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/data/";
	string completeLinkageFilePath = "/home/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/data/";

	// string filePath = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/";
    // string singleLinkageFilePath = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/data/";
	// string completeLinkageFilePath = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/data/";

	filePath = filePath + argv[1];
	string singleLinkageFileName = argv[2];
	singleLinkageFilePath = singleLinkageFilePath + argv[2];
	string completeLinkageClusterFile = "kNN_Global_Linked_" + singleLinkageFileName;
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