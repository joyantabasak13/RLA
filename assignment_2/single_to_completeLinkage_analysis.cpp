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

// #include <bits/stdc++.h>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

using namespace std;

int threshold = 99;
int totalRecords;
int lenMax;
int attributes;
int blockFieldIndex = 4;
int blockingDistanceThreshold = 2;
int singleNonBlockingDistanceThreshold = 5;
int cumulativeNonBlockingDistanceThreshold = 10;

vector<vector<string> > vec2D;
vector<vector<int> > cluster2D;
map<int, vector<int>> completeClusters; 

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

// I/O Functions
int getDistance(vector<string>& vec1, vector<string>& vec2){
	int dist = 0;
	for(int i=1; i<attributes-1; i++){
		dist = dist + calculateBasicED(vec1[i], vec2[i], threshold);
	}
	return dist;
}

void getData(string& filePath) {
    string line;
    ifstream records(filePath);

    int ind = 0;
    while (getline (records, line)) {
        vector<string> result;
        boost::split(result, line, boost::is_any_of(","));
        vector<string> vec;
        vec.push_back(result[0]);
		for(int i = 1; i<result.size(); i++){
			boost::to_lower(result[i]);
			vec.push_back(result[i]);
		}
        vec2D.push_back(vec);
    }
    records.close();

    attributes = vec2D[0].size();
    cout<< "Records read: " << vec2D.size() << endl;
	cout<< "Attributes: "<<attributes << endl;

}

void getSingleLinakgeClusters(string& filePath) {
	string line;
    ifstream clusters(filePath);

    int rec = 0;
    while (getline (clusters, line)) {
        vector<string> vec;
		vector<int> cluster;
        boost::split(vec, line, boost::is_any_of(","));
		for(int i = 0; i < vec.size()-1; i++) {
			// cout<< vec[i] << endl;
			cluster.push_back(stoi(vec[i]));
			rec++;
		}
        cluster2D.push_back(cluster);
    }
    clusters.close();
	cout<< "Total Single linkage Clusters: " << cluster2D.size() << endl;
	cout<< "Total records found in Clusters: " << rec << endl;
}

void sourceConsistentCompleteLinkage(){
	for(int i=0; i<cluster2D.size(); i++){
		int numComponents = cluster2D[i].size();
		// cout<< "Cluster: " << i << " with Components: " << numComponents << endl;
		vector<vector<string>> records(numComponents);
		UnionFind uf_finalClusters;
		uf_finalClusters.setVariable(numComponents);
		vector<set<string>> clusterSources(numComponents);
		// copy records
		for(int j=0; j < numComponents; j++){
			records[j] = vec2D[cluster2D[i][j]];
			clusterSources[j].insert(records[j][attributes-1]);
		}
		// sort records by source
		for (int j = 0; j < numComponents; j++)
		{
			for (int k = 0 ; k < numComponents - j - 1; k++) {
				if (records[k][attributes-1] > records[k+1][attributes-1])
				{
					// cout<< "K " << k << " : " << records[k][attributes-1]<< endl << " and K+1 " << k+1 << " : " << records[k+1][attributes-1] << endl;

					records[k+1].swap(records[k]);
				}
			}
		}

		// compute all pair distance
		priority_queue<pair<int, pair<int,int> > > pq;
		for (int j = 0; j < numComponents; j++)
		{
			for (int k = j + 1 ; k < numComponents; k++) {
				if(records[j][attributes-1] != records[k][attributes-1]) {
					int dist = getDistance(records[j],records[k]);
					pair<int, int> uv;
					uv.first = j;
					uv.second = k;
					pair<int, pair<int,int>> edge;
					edge.first = dist;
					edge.second = uv;
					pq.push(edge);
				}
			}
		}
		while(pq.empty()==false) {
			int u = pq.top().second.first;
			int x = pq.top().second.second;
			int root_u = uf_finalClusters.find(u);
			int root_v = uf_finalClusters.find(x);
			if (root_u != root_v) {
				bool noCommon = true;
				for (auto v : clusterSources[root_u]) {
					if (clusterSources[root_v].count(v) > 0) {
						noCommon = false;
					}
				}
				if (noCommon)
				{
					uf_finalClusters.weightedUnion(root_u, root_v);
					if (root_u == uf_finalClusters.find(root_u))
					{
						for (auto v : clusterSources[root_v]) {
							clusterSources[root_u].insert(v);
						}
					} else {
						for (auto u : clusterSources[root_u]) {
							clusterSources[root_v].insert(u);
						}
					}
				}
				
			}

			// cout<< "Edge weight: " << pq.top().first << " Edge u: " << pq.top().second.first << " Edge v: " << pq.top().second.second << endl;
			pq.pop();
		}

		for(int j=0; j<numComponents; j++) {
			int root = uf_finalClusters.find(j);
			completeClusters[cluster2D[i][root]].push_back(cluster2D[i][j]);
		}

		if(i%100000 == 0) {
			cout<< "Processed " << i << " Single Linkage Clusters" << endl;
		}

	}
}

void writeCompleteConnectedComponentToFile(string& result_file_name) {
	ofstream out_file;
    out_file.open(result_file_name);
	for (auto const& p : completeClusters) {
        for (int i=0; i<p.second.size(); i++) {
			out_file<< vec2D[p.second[i]][0] << ",";
        }
        out_file<< "\n";
	}
	out_file.close();
}

void getCompleteClusters(string& fileName) {
	sourceConsistentCompleteLinkage();
	writeCompleteConnectedComponentToFile(fileName);
	cout<< "Total Clusters " << completeClusters.size() << endl;
}


int main(int argc, char** argv) {

	// IO PATHS
    string filePath = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/";
    string fileName = argv[1];
	string singleLinkageClusterFile = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/Server_results/genRLA_NC/out_SB_RLA_SingleLinkage_RecInd_NO_DEDUP_NC_voterData_5M_Source_Annotated.csv_pGEN_NC_lastName_6_threads_dist_1_9";
    string completeLinkageClusterFile = singleLinkageClusterFile + "_COMPLETE";
	filePath = filePath + argv[1];
	getData(filePath);
	totalRecords = vec2D.size();
	getSingleLinakgeClusters(singleLinkageClusterFile);

	// // Outputs
    // string out_file_path = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/data/";
	// string out_file_name = "CompleteLinkage_NC_CLIPLike_Clustering";

	// Sort the Combined Data
	clock_t currTS_p0	= clock();
    getCompleteClusters(completeLinkageClusterFile);
    double sorting_p0_t	= (double)(clock() - currTS_p0) / CLOCKS_PER_SEC;
    cout<< "Complete Linkage Time "<< sorting_p0_t << endl;

    return 0;
}