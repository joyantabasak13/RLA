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
int matArr[50][50] = {0};

vector<string> vec1D;
vector<vector<string> > vec2D;
vector<vector<vector<string> > > sameEntities;
vector<pair<int, string> > combinedData;
vector<vector<vector<int> > > distances;
vector<vector<pair<int,int>>> distCompPairs;
vector<vector<int>> singleLinkageDistance;


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
					if(i>1 && j>1) {
						if((str1[i-1]==str2[j-2]) && (str1[i-2] == str2[j-1])) {
							matArr[i][j] = min(matArr[i][j], matArr[i-2][j-2]+1);
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
					else {
						matArr[i][j] = min(min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
						if((i-2>=0) && (j + (i - dist) - 2 >=0)) {
							if (((int)s1[i - 2] == (int)s2[j + (i - dist) - 1]) && ((int)s1[i - 1] == (int)s2[j + (i - dist) - 2])) {
								matArr[i][j] = min(matArr[i][j], matArr[i-2][j+2]+1);
							}
						}
					}
							
				}
				else
				{
					if ((int)s1[i - 1] == (int)s2[j - (dist - i) - 1])
						matArr[i][j] = min(matArr[i - 1][j], matArr[i][j - 1] + 1);
					else {
						matArr[i][j] = min(matArr[i - 1][j] + 1, matArr[i][j - 1] + 1);
						if((i-2>=0) && (j + (i - dist) - 2 >=0)) {
							if (((int)s1[i - 2] == (int)s2[j + (i - dist) - 1]) && ((int)s1[i - 1] == (int)s2[j + (i - dist) - 2])) {
								matArr[i][j] = min(matArr[i][j], matArr[i-2][j+2]+1);
							}
						}
					}
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
					else {
						matArr[i][j] = min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1);
						if((i-2>=0) && (j + (i - dist) - 2 >=0)) {
							if (((int)s1[i - 2] == (int)s2[j + (i - dist) - 1]) && ((int)s1[i - 1] == (int)s2[j + (i - dist) - 2])) {
								matArr[i][j] = min(matArr[i][j], matArr[i-2][j+2]+1);
							}
						}
					}
				}
				else if (j != (col - 1))
				{
					if ((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j] = min(min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
					else {
						matArr[i][j] = min(min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
						if((i-2>=0) && (j + (i - dist) - 2 >=0)) {
							if (((int)s1[i - 2] == (int)s2[j + (i - dist) - 1]) && ((int)s1[i - 1] == (int)s2[j + (i - dist) - 2])) {
								matArr[i][j] = min(matArr[i][j], matArr[i-2][j+2]+1);
							}
						}
					}
				}
				else
				{
					if ((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j] = min(matArr[i - 1][j], matArr[i][j - 1] + 1);
					else {
						matArr[i][j] = min(matArr[i - 1][j] + 1, matArr[i][j - 1] + 1);
						if((i-2>=0) && (j + (i - dist) - 2 >=0)) {
							if (((int)s1[i - 2] == (int)s2[j + (i - dist) - 1]) && ((int)s1[i - 1] == (int)s2[j + (i - dist) - 2])) {
								matArr[i][j] = min(matArr[i][j], matArr[i-2][j+2]+1);
							}
						}
					}
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
					else {
						matArr[i][j] = min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1);
						if((i-2>=0) && (j + (i - dist) - 2 >=0)) {
							if (((int)s1[i - 2] == (int)s2[j + (i - dist) - 1]) && ((int)s1[i - 1] == (int)s2[j + (i - dist) - 2])) {
								matArr[i][j] = min(matArr[i][j], matArr[i-2][j+2]+1);
							}
						}
					}
				}
				else
				{
					if ((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j] = min(min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
					else {
						matArr[i][j] = min(min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
						if((i-2>=0) && (j + (i - dist) - 2 >=0)) {
							if (((int)s1[i - 2] == (int)s2[j + (i - dist) - 1]) && ((int)s1[i - 1] == (int)s2[j + (i - dist) - 2])) {
								matArr[i][j] = min(matArr[i][j], matArr[i-2][j+2]+1);
							}
						}
					}
				}
				if ((j == diagonal) && (matArr[i][j] > dist))
					return threshold + 1;
			}
		}
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
		for(int i = 0; i<result.size(); i++){
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

void getDistances() {
	for(int i = 0; i< sameEntities.size(); i++) {
		bool printThis = false;
	
		vector<vector<int>> intraClusterDists;
		vector<pair<int, int>> intraClusterPairs;
		for (int j = 0; j < sameEntities[i].size()-1; j++)
		{
			for (int k = j+1 ; k < sameEntities[i].size(); k++)
			{
				vector<int> pairwiseAttrDists;
				pair<int,int> pairDist;
				pairDist.first = j;
				pairDist.second = k;
				pairwiseAttrDists.resize(attributes-1, 0);
				for (int attr = 1; attr < attributes-1; attr++)
				{
					int attrDist = calculateBasicED(sameEntities[i][j][attr], sameEntities[i][k][attr], 10);
					pairwiseAttrDists[attr-1] = attrDist;
				}
				intraClusterDists.push_back(pairwiseAttrDists);
				intraClusterPairs.push_back(pairDist);
			}
		}
		distances.push_back(intraClusterDists);
		distCompPairs.push_back(intraClusterPairs);
	}
}

void getMaxDistanceInEachAttr() {
	vector<int> maxDistanceSummary;
	maxDistanceSummary.resize(attributes-1);
	for (int i=0; i<attributes-1; i++) {
		maxDistanceSummary[i] = 0;
	}
	int maxCumulativeDistance = 0;
	for (int i=0; i<distances.size(); i++) {
		for (int j=0; j<distances[i].size(); i++) {
			int cumulativeDistance = 0;
			for (int k = 0; k<attributes-1; k++) {
				cumulativeDistance += distances[i][j][k];
				if (distances[i][j][k] > maxDistanceSummary[k])
				{
					maxDistanceSummary[k] = distances[i][j][k];
				}
			}
			if (cumulativeDistance > maxCumulativeDistance) {
				maxCumulativeDistance = cumulativeDistance;
			}
		}
	}
	cout<< "Max Cumulative distance in a Group: " << maxCumulativeDistance << endl;
	for (int i = 0; i < attributes-1; i++)
	{
		cout<< "MaxDist in Atrribute: "<< i << " is " << maxDistanceSummary[i] << endl;
	}	
}

void getAllDistanceInEachAttr() {
	vector<map<int, int>> allDistanceSummary;
	allDistanceSummary.resize(attributes);

	for (int i=0; i<distances.size(); i++) {
		for (int j=0; j<distances[i].size(); i++) {
			int cumulativeDistance = 0;
			for (int k = 0; k<attributes-1; k++) {
				cumulativeDistance += distances[i][j][k];
				if (!allDistanceSummary[k].count(distances[i][j][k])) {
					allDistanceSummary[k][distances[i][j][k]] = 0;
				}
				allDistanceSummary[k][distances[i][j][k]]++;
			}
			if (!allDistanceSummary[attributes-1].count(cumulativeDistance)) {
					allDistanceSummary[attributes-1][cumulativeDistance] = 0;
				}
			allDistanceSummary[attributes-1][cumulativeDistance]++;
		}
	}
	for(int i=0; i<attributes ; i++) {
			cout<< "COLUMN " << i << " Distances " << endl;
		for(auto it = allDistanceSummary[i].cbegin(); it != allDistanceSummary[i].cend(); ++it)
		{
			std::cout << "Distance: " << it->first << " has pair " << it->second << endl;
		}
	}	
}

int getCumulativeDist(vector<int>& distanceVec){
	int cumDist = 0;
	for(int i = 0; i< attributes -1; i++) {
		cumDist += distanceVec[i];
	}
	return cumDist;
}

int getMaxCutMinSpanTree(int clusterNum){
	int clusterSize = sameEntities[clusterNum].size();
	UnionFind uf;
	uf.setVariable(clusterSize);
	priority_queue<pair<int, pair<int, pair<int,int>>>> pq;
	int maxCutEdge = -1;
	// Find the minimum spanning tree
	// Push all edges to a priority queue
	for (int j=0; j<distances[clusterNum].size(); j++) {
		int u = distCompPairs[clusterNum][j].first;
		int v = distCompPairs[clusterNum][j].second;
		int edgeWeight = getCumulativeDist(distances[clusterNum][j]);
		pair<int, pair<int, pair<int,int>>> edge;
		edge.first = edgeWeight*(-1);
		edge.second.first = j;
		edge.second.second.first = u;
		edge.second.second.second = v;
		pq.push(edge);
	}
	// cout<< "Got " << pq.size() << " Edges" << endl;
	// Do Kruskal
	while(pq.empty()==false){
		int u = pq.top().second.second.first;
		int v = pq.top().second.second.second;
		int ind = pq.top().second.first;
		int root_u = uf.find(u);
		int root_v = uf.find(v);
		if (root_u != root_v) {
			uf.weightedUnion(root_u, root_v);
			maxCutEdge = ind;
		}
		if (uf.getRootVal(root_u) + clusterSize == 0) {
			break;
		}
		pq.pop();
	}
	// if (maxCutEdge >= 0) {
	// 	// cout<< "MaxCut Size: " << distances[clusterNum][maxCutEdge].size() << endl;
	// } else {
	// 	// cout<< "Cluster has no edge!" << endl;
	// }
	
	return maxCutEdge;
}

void getSingleLinakgeSummary(vector<vector<int>> SLDistances) {
	vector<map<int, int>> SLDistanceSummary;
	SLDistanceSummary.resize(attributes);

	for (int i=0; i<SLDistances.size(); i++) {
		int cumulativeDistance = 0;
		for (int k = 0; k<attributes-1; k++) {
			cumulativeDistance += SLDistances[i][k];
			if (!SLDistanceSummary[k].count(SLDistances[i][k])) {
				SLDistanceSummary[k][SLDistances[i][k]] = 0;
			}
			SLDistanceSummary[k][SLDistances[i][k]]++;
		}
		if (!SLDistanceSummary[attributes-1].count(cumulativeDistance)) {
				SLDistanceSummary[attributes-1][cumulativeDistance] = 0;
			}
		SLDistanceSummary[attributes-1][cumulativeDistance]++;
	}

	for(int i=0; i<attributes ; i++) {
			cout<< "COLUMN " << i << " Distances " << endl;
		for(auto it = SLDistanceSummary[i].cbegin(); it != SLDistanceSummary[i].cend(); ++it)
		{
			std::cout << "Distance: " << it->first << " has pair " << it->second << endl;
		}
	}	
}

void checkThresholdValidity(vector<int> attrThresholds, int cumulativeThreshold) {
	int validCount = 0;
	int invalidCount = 0;
	int numRecordsInValidCounts = 0;

	for (int i=0; i<singleLinkageDistance.size(); i++) {
		int cumulativeDistance = 0;
		bool isValid = true;
		for (int k = 0; k<attributes-2; k++) {
			cumulativeDistance += singleLinkageDistance[i][k];
			if ((singleLinkageDistance[i][k] > attrThresholds[k])||(cumulativeDistance > cumulativeThreshold)) {
				// cout<< "Error for clusterID " << sameEntities[i][0][0] << " With distance: " <<  singleLinkageDistance[i][k] << " in attr: " << k << endl;
				// cout<< "CumulativeDistance is: " << cumulativeDistance << endl;
				isValid = false;
				invalidCount++;
				break;
			}
		}
		if (isValid) {
				validCount++;
				numRecordsInValidCounts += sameEntities[i].size();
		}
	}
	cout<< "Total Clusters: " << singleLinkageDistance.size() << endl;
	cout<< "Valid Clusters: " << validCount << endl;
	cout<< "Invalid Clusters: " << invalidCount << endl;
	cout<< "Records in valid clusters: " << numRecordsInValidCounts << endl;
}

void getSingleLinakgeDistance() {
	singleLinkageDistance.resize(distances.size());
	for (int i=0; i<distances.size(); i++) {
		singleLinkageDistance[i].resize(attributes-1);
	}
	cout<< "Initialized singleLinakge vectors" << endl;

	for (int i=0; i<distances.size(); i++) {

		vector<int> maxEdgeWeights;
		int maxCutEdge;
		// cout<< "Getting the MaxCutMinSpaningTree for cluster: " << i << " That has " << sameEntities[i].size() << " records" << endl;
		maxCutEdge = getMaxCutMinSpanTree(i);
		if (maxCutEdge >= 0) {
			maxEdgeWeights = distances[i][getMaxCutMinSpanTree(i)];
		} else {
			maxEdgeWeights.resize(attributes-1);
			for (int k=0; k < attributes-1; k++) {
				maxEdgeWeights[k] = 0;
			}
		}
		// cout<< "Got maxCutMinTree for cluster: " << i << endl; 
		// cout<< "maxEdgeWeights size: " << maxEdgeWeights.size() << endl;
		// cout<< "And total attributes: " << attributes << endl;
		for (int j = 0; j < attributes-1; j++)
		{
			singleLinkageDistance[i][j] = maxEdgeWeights[j];
		}
		
	}
	
	int max = 0;
	int maxInd = 0;
	for(int i = 0; i<singleLinkageDistance.size(); i++) {
		int cumDist = getCumulativeDist(singleLinkageDistance[i]);
		if (cumDist>max) {
			max = cumDist;
			maxInd = i;
		}
	}
	cout<< "Max Dist: " << max << " AT cluster: " << maxInd << " With id: " << sameEntities[maxInd][0][0] << endl;
	for (int i = 0; i < singleLinkageDistance[maxInd].size(); i++)
	{
		cout<< "Attribute: " << i << " Dist: " << singleLinkageDistance[maxInd][i] << endl;
	}
	// getSingleLinakgeSummary(singleLinkageDistance);
}

bool checkIfConnected(int clusterID, vector<int>& attrThresholds, int cumulativeThreshold) {
	// cout<< "start" << endl;
	UnionFind uf;
	uf.setVariable(sameEntities[clusterID].size());

	for(int i = 0; i<distances[clusterID].size(); i++) {
		bool isOK = true;
		int thisCumDist = 0;
			
		for(int k=0; k<attrThresholds.size(); k++){
			if (distances[clusterID][i][k] > attrThresholds[k]) {
				isOK = false;
				break;
			}
			thisCumDist += distances[clusterID][i][k];
		}
		if(isOK) {
			if (thisCumDist <= cumulativeThreshold) {
				uf.weightedUnion(distCompPairs[clusterID][i].first, distCompPairs[clusterID][i].second);
			}
		}
		
	}

	int numCluster = uf.getSetCount();
	if(numCluster == 1) {
		return true;
	} else {
		return false;
	}
}

pair<double, double> calculateThesholdCoverage(vector<int>& attrThresholds, int cumulativeThreshold){
	int coveredClusters = 0;
	int coveredRecords = 0;
	for(int i=0; i<sameEntities.size(); i++) {
		bool isConnected = checkIfConnected(i, attrThresholds, cumulativeThreshold);
		// bool isConnected = true;
		if (isConnected) {
			coveredClusters++;
			coveredRecords += sameEntities[i].size();
		}
	}
	cout<< "Total Clusters: " << sameEntities.size() << endl;
	cout<< "Covered Clusters: " << coveredClusters << endl;
	cout<< "Covered Record: " << coveredRecords << " %Total: " << ((double)coveredRecords)/((double)vec2D.size());
	cout<< "Coverage: " << (double)(((double)coveredClusters) / ((double)sameEntities.size())) << endl;
	pair<double, double> p;
	p.first = ((double)coveredRecords)/((double)vec2D.size());
	p.second = (double)(((double)coveredClusters) / ((double)sameEntities.size()));
	return p;
}

int main(int argc, char** argv) {

	// IO PATHS
    string filePath = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/";
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

	clock_t currTS_p1	= clock();
    getSameEntities();
	double exactClustering_p1_t	= (double)(clock() - currTS_p1) / CLOCKS_PER_SEC;
    cout<< "Same Entity Extraction time: "<< exactClustering_p1_t << endl;
	
	clock_t currTS_p2	= clock();
	getDistances();
	double distanceCalculation_p2_t	= (double)(clock() - currTS_p2) / CLOCKS_PER_SEC;
    cout<< "Intra Cluster Distance Calculation time: "<< distanceCalculation_p2_t << endl;

	clock_t currTS_p3	= clock();
	// getAllDistanceInEachAttr();
	// getSingleLinakgeDistance();
	// double distanceCalculation_p3_t	= (double)(clock() - currTS_p3) / CLOCKS_PER_SEC;
    // cout<< "MaxCutMinSpan tree calculation time: "<< distanceCalculation_p3_t << endl;

	//Add the theshold amd cumulative distance
	int limit = 4;
	int cumulativeLimit = 7;

	for (int i = 0; i < limit; i++){
		for (int j=0; j< limit; j++) {
			for (int k=0; k<limit; k++) {
				for(int m = 0; m < limit; m++) {
					for (int n = 0; n < cumulativeLimit; n++){
						int cumulativeThreshold = n;
						vector<int> attrThreshold{i,j,k,m};
						cout<< "i: " << i << " j: " << j << " k: " << k << " m: " << m << " n: " << n << endl;
						clock_t start_t	= clock();
						pair<double, double> result = calculateThesholdCoverage(attrThreshold, cumulativeThreshold);
						if(result.first > 0.999 && result.second > 0.999) {
							return 0;
						}
						double singleIterTime	= (double)(clock() - start_t) / CLOCKS_PER_SEC;
    					cout<< "time "<< singleIterTime << endl;
					}	
				}
			}
		}
	}
	
	// vector<int> attrThresholds{3,3,3,3};
	// int cumulativeThreshold = 4;
	// // checkThresholdValidity(attrThresholds, cumulativeThreshold);
	// calculateThesholdCoverage(attrThresholds, cumulativeThreshold);
	double totalTime	= (double)(clock() - currTS_p0) / CLOCKS_PER_SEC;
    cout<< "Total time "<< totalTime << endl;


	string x = "XYZ10";
	boost::to_lower(x);
	cout<< x << endl;
	string y = "";
	cout<< "Distance is " << calculateBasicED(x, y, threshold) << endl; 

    return 0;
}