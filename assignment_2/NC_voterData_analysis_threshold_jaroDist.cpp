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

int totalRecords;
int lenMax;
int attributes;

vector<vector<string> > vec2D;
vector<vector<vector<string> > > sameEntities;
vector<pair<int, string> > combinedData;
vector<vector<vector<double> > > distances;
vector<vector<pair<int,int>>> distCompPairs;


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


double jaro_distance(string& s1, string& s2)
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
	return (((double)match) / ((double)len1)
			+ ((double)match) / ((double)len2)
			+ ((double)match - t) / ((double)match))
		/ 3.0;
}

double calculateJaroWinklerDist(string& str1,string& str2) {
	int min_string_len = (str1.size() - str2.size());

	if (min_string_len < 0){
		min_string_len = -min_string_len;
	}
	if (min_string_len > 1) {
		return 0.20;
	}

	double jaro_dist = jaro_distance(str1, str2);

	if (jaro_dist > 0.0) {
		int prefix = 0;

		for (int i = 0; i < min(str1.length(), str2.length()); i++) {
			// If the characters match
			if (str1[i] == str2[i])
				prefix++;

			// Else break
			else
				break;
		}

		// Maximum of 4 characters are allowed in prefix
		prefix = min(4, prefix);

		// Calculate jaro winkler Similarity
		jaro_dist += 0.1 * prefix * (1 - jaro_dist);
	}
	return jaro_dist;

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

void getDistances() {
	for(int i = 0; i< sameEntities.size(); i++) {
		bool printThis = false;
	
		vector<vector<double>> intraClusterDists;
		vector<pair<int, int>> intraClusterPairs;
		for (int j = 0; j < sameEntities[i].size()-1; j++)
		{
			for (int k = j+1 ; k < sameEntities[i].size(); k++)
			{
				vector<double> pairwiseAttrDists;
				pair<int,int> pairDist;
				pairDist.first = j;
				pairDist.second = k;
				pairwiseAttrDists.resize(attributes-1, 0);
				for (int attr = 1; attr < attributes-1; attr++)
				{
					// double attrDist = jaro_distance(sameEntities[i][j][attr], sameEntities[i][k][attr]);
					double attrDist = calculateJaroWinklerDist(sameEntities[i][j][attr], sameEntities[i][k][attr]);
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

bool checkIfConnected(int clusterID, vector<double>& attrThresholds, double avgThreshold) {
	// cout<< "start" << endl;
	UnionFind uf;
	uf.setVariable(sameEntities[clusterID].size());

	for(int i = 0; i<distances[clusterID].size(); i++) {
		bool isOK = true;
		double thisTotSim = 0.0;
			
		for(int k=0; k<attrThresholds.size(); k++){
			if (distances[clusterID][i][k] < attrThresholds[k]) {
				isOK = false;
				break;
			}
			thisTotSim += distances[clusterID][i][k];
		}
		if(isOK) {
			if ((double)(thisTotSim/((double)attrThresholds.size())) > avgThreshold) {
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

pair<double,double> calculateThesholdCoverage(vector<double>& attrThresholds, double avgThreshold){
	int coveredClusters = 0;
	int coveredRecords = 0;
	for(int i=0; i<sameEntities.size(); i++) {
		bool isConnected = checkIfConnected(i, attrThresholds, avgThreshold);
		if (isConnected) {
			coveredClusters++;
			coveredRecords += sameEntities[i].size();
		}
	}
	cout<< "Total Clusters: " << sameEntities.size() << endl;
	cout<< "Covered Clusters: " << coveredClusters << endl;
	cout<< "Covered Record: " << coveredRecords << " %Total: " << ((double)coveredRecords)/((double)vec2D.size());
	cout<< "Coverage: " << (double)(((double)coveredClusters) / ((double)sameEntities.size())) << endl;
	pair<double, double> covs;
	covs.first = ((double)coveredRecords)/((double)vec2D.size());
	covs.second = (double)(((double)coveredClusters) / ((double)sameEntities.size()));
	return covs;
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

	clock_t currTS_p1	= clock();
    getSameEntities();
	double exactClustering_p1_t	= (double)(clock() - currTS_p1) / CLOCKS_PER_SEC;
    cout<< "Same Entity Extraction time: "<< exactClustering_p1_t << endl;
	
	clock_t currTS_p2	= clock();
	getDistances();
	double distanceCalculation_p2_t	= (double)(clock() - currTS_p2) / CLOCKS_PER_SEC;
    cout<< "Intra Cluster Distance Calculation time: "<< distanceCalculation_p2_t << endl;

	clock_t currTS_p3	= clock();

	//Add the theshold amd cumulative distance
	double limit = 0.74;
	double avgLimit = 0.74;
	vector<double> bestThresholds;
	double bestAvgThreshold;
	double bestCoverage = 0.995;

	for (double i = .90; i > limit; i-=0.05){
		for (double j = .90; j > limit; j-=0.05) {
			for (double k = .90; k > limit; k-=0.05) {
				for(double m = .90; m > limit; m-=0.05) {
					for (double n = .90; n > avgLimit; n-=0.05){
						double avgThreshold = n;
						vector<double> attrThreshold{i,j,k,m};
						cout<< endl;
						cout<< "i: " << i << " j: " << j << " k: " << k << " m: " << m << " n: " << n << endl;
						clock_t start_t	= clock();
						pair<double,double> curCov = calculateThesholdCoverage(attrThreshold, avgThreshold);
						if((curCov.first > bestCoverage) && (curCov.second > bestCoverage)){
							double singleIterTime	= (double)(clock() - start_t) / CLOCKS_PER_SEC;
    						cout<< "time "<< singleIterTime << endl;
							return 0;
						}
						double singleIterTime	= (double)(clock() - start_t) / CLOCKS_PER_SEC;
    					cout<< "time "<< singleIterTime << endl;
					}	
				}
			}
		}
	}

	cout<< "Best Coverage: " << bestCoverage << endl;
	cout<< "Best Avg Threshold: " << bestAvgThreshold << endl;
	cout<< "Best Thresholds: " << endl;
	for (size_t i = 0; i < bestThresholds.size(); i++)
	{
		cout<< "Attribute: " << i << " has threshold: " << bestThresholds[i] << endl;
	}
	
	double totalTime	= (double)(clock() - currTS_p0) / CLOCKS_PER_SEC;
    cout<< "Total time "<< totalTime << endl;

	string x = "farmville";
	string y = "farmVelel";
	boost::to_lower(x);
	boost::to_lower(y);

	cout<< "Distance is " << jaro_distance(x, y) << endl; 

    return 0;
}