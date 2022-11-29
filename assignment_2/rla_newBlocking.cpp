#include <iostream>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <set>
#include <tuple>
#include <utility>
#include <algorithm>
#include <map>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

using namespace std;

int threshold = 99;
int totalRecords = 50000;
int lenMax;
int totalUniqueRecords;
int attributes;
int extraEdges = 0;

vector<int> blockfieldArr;
// vector<string> uniqueblockfieldArr;
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

vector<vector<int> > edgeArr;
set<pair<int, int> > set_of_edges;

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

		//Remove digits and alphanumerics from name
		auto last = std::remove_if(result[1].begin(), result[1].end(), [](auto ch) {
        								return ::isdigit(ch) || ::ispunct(ch) || ::iswpunct(ch);
    								});
		result[1].erase(last, result[1].end()); //Remove junk left by remove_if() at the end of iterator
        boost::to_lower(result[1]);
		vec.push_back(result[1]);
		vec.push_back(result[2]);
		vec.push_back(result[3]);
        vec2D.push_back(vec);
    }
    records.close();
    // vector<vector<string>> temp2D;

    // temp2D.resize(vec2D[0].size(), vec2D.size());
    // temp2D.resize(vec2D[0].size());
    // for (int i = 0; i < vec2D[0].size(); ++i) temp2D[i].resize(vec2D.size());

    vec1D.resize(vec2D[0].size()*vec2D.size());


    // for (size_t i = 0; i < vec2D.size(); i++)
    // {
    //     for(size_t j = 0; j< vec2D[0].size(); j++) {
    //         temp2D[j][i] = vec2D[i][j];
    //     }
    // }
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
	cout << endl << "total exact clusters: " << totalUniqueRecords << endl;
}

void getUniqueEntries() {
	uniqueRecords.resize(totalUniqueRecords);

	for (size_t i = 0; i < totalUniqueRecords; i++)
    {
        uniqueRecords[i] = combinedData[exactmatches[i][0]];
        // cout<< uniqueRecords[i].first << "\t" << uniqueRecords[i].second << endl;
    }
}

void doSortedComp() {
	headlessCopies.resize(2*totalUniqueRecords);
	for(int i=0; i< uniqueRecords.size(); i++) { 
		headlessCopies[i].first = i;
		headlessCopies[i].second = uniqueRecords[i].second;
		headlessCopies[totalUniqueRecords+i].first = i;
		headlessCopies[totalUniqueRecords+i].second = uniqueRecords[i].second.substr(1,uniqueRecords[i].second.size()-1) + '0';
		// cout<< uniqueRecords[i].second << endl;
		// cout<< headlessCopies[totalUniqueRecords+i].second << endl;
	}
	radixSort(headlessCopies);
	for (int i = 1; i < headlessCopies.size(); i++) {
		// cout<< headlessCopies[i].second << endl;
		if (headlessCopies[i-1].second.compare(headlessCopies[i].second) == 0) {
			// cout<< headlessCopies[i-1].first << " " << headlessCopies[i-1].second << endl;
			// cout<< headlessCopies[i].first << " " << headlessCopies[i].second << endl;
			pair<int, int> edge_pair;
            int i_th_record_id = headlessCopies[i-1].first;
            int j_th_record_id = headlessCopies[i].first;
            if (i_th_record_id < j_th_record_id) {
                edge_pair.first = i_th_record_id;
                edge_pair.second = j_th_record_id;
            } else {
                edge_pair.first = j_th_record_id;
                edge_pair.second = i_th_record_id;
            }
			
            if (!set_of_edges.count(edge_pair)) {
                set_of_edges.insert(edge_pair);
				extraEdges++;
            }
		}
	}
	cout<< "Edges Added: "<< extraEdges << endl;
}

void doNormalBlocking() {
	int base = 26;
	int kmer = 3;
    int blockID = 0;
	int total_blocked = 0;
    int unique_blocked = 0;
    int total_str_size = 0;
	int blockTotal = pow(base,kmer);
	block_list.resize(blockTotal);

	cout<< "Total unique records: " << totalUniqueRecords << endl;

	for (int i = 0; i < totalUniqueRecords; i++)
	{
		string blockingStr = vec1D[uniqueRecords[i].first*attributes + 1];
        // string blockingStr = vec2D[uniqueRecords[i].first][1];
        total_str_size += blockingStr.size();
        // cout<< i << "\t" << blockingStr << "\t" << blockingStr.size() << endl;
		for (int j = 0; j < blockingStr.size() - kmer + 1 ; ++j)
		{
			blockID = 0;

			for (int k = 0; k < kmer; ++k)
			{
				blockID += ((int)blockingStr[j+k] - 97) * pow(base,k);
			}
			if(!block_list[blockID].count(i)){
                block_list[blockID].insert(i);
                unique_blocked++;
            }
			total_blocked++;
		}
	}
	cout<< "Total blocked: " << total_blocked << endl;
    cout<< "Uniquely blocked: " << unique_blocked << endl;
    cout<< "Total string size covered: " << total_str_size << endl;
    cout<< "Expected blocks:" << total_str_size - (kmer-1)*totalUniqueRecords << endl;
}

void doSuperBlocking() {
	int base = 26;
	int kmer = 3;
    int blockID = 0;
	int total_blocked = 0;
    int unique_blocked = 0;
    int total_str_size = 0;
	int perAplhaBlocks = pow(base,kmer);
	int blockTotal = base * perAplhaBlocks;
	block_list.resize(blockTotal);

	cout<< "Total unique records: " << totalUniqueRecords << endl;

	for (int i = 0; i < totalUniqueRecords; i++)
	{
		string blockingStr = vec1D[uniqueRecords[i].first*attributes + 1];
        // string blockingStr = vec2D[uniqueRecords[i].first][1];
        total_str_size += blockingStr.size();
        // cout<< i << "\t" << blockingStr << "\t" << blockingStr.size() << endl;
		for (int j = 0; j < blockingStr.size() - kmer + 1 ; ++j)
		{
			blockID = 0;

			for (int k = 0; k < kmer; ++k)
			{
				blockID += ((int)blockingStr[j+k] - 97) * pow(base,k);
			}
			blockID = (blockingStr[0]-97)*perAplhaBlocks + blockID;
			if(!block_list[blockID].count(i)){
                block_list[blockID].insert(i);
                unique_blocked++;
            }
			total_blocked++;
		}
	}
	cout<< "Total blocked: " << total_blocked << endl;
    cout<< "Uniquely blocked: " << unique_blocked << endl;
    cout<< "Total string size covered: " << total_str_size << endl;
    cout<< "Expected blocks:" << total_str_size - (kmer-1)*totalUniqueRecords << endl;
}

bool isLinkageOk(vector<string> &a, vector<string> &b, int threshold)
{
    //int dist = 0;
	int name_dist = calculateBasicED(a[1], b[1], 1);
    //dist +=name_dist;
    if (name_dist <= threshold) {
        int dod_dist = calculateBasicED(a[2], b[2], 1);
        //dist+=dod_dist;   
        if (dod_dist <= threshold) {
            int dob_dist = calculateBasicED(a[3], b[3], 1);
            //dist+=dod_dist;
            // if(dist==0) {
            //     //Self edge?
            //     return false;
            // }
            if (dob_dist <= threshold) {
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    } else {
        return false;
    }
}

void generateEdgilist(set<int>& blockRowArr)
{
	int blockItemTotal	= blockRowArr.size();

	vector<vector<string>> dataArr(blockItemTotal); // to make cache-efficient, keep records in a row
	vector<int> uniqueIndices(blockItemTotal);
    int l = 0;
    BOOST_FOREACH(int p, blockRowArr) {
        dataArr[l] = vec2D[uniqueRecords[p].first];
        uniqueIndices[l] = p;
        l++;
    };

	for (int i = 0; i < blockItemTotal-1; i++)
	{
        for (int j = i+1; j < blockItemTotal; j++)
        {
            pair<int, int> edge_pair;
            int i_th_record_id = uniqueIndices[i];
            int j_th_record_id = uniqueIndices[j];
            if (i_th_record_id < j_th_record_id) {
                edge_pair.first = i_th_record_id;
                edge_pair.second = j_th_record_id;
            } else {
                edge_pair.first = j_th_record_id;
                edge_pair.second = i_th_record_id;
            }

            if (!set_of_edges.count(edge_pair)) {
                if (isLinkageOk(dataArr[i], dataArr[j], 1)) {
                    set_of_edges.insert(edge_pair);
                }
            }
        }
        
	}

	dataArr.clear();
    uniqueIndices.clear();
}

void printEdges() {
    ofstream log_file;
	log_file.open("log_file");
    BOOST_FOREACH(pair p, set_of_edges) {
        string U	= uniqueRecords[p.first].second;
		string V	= uniqueRecords[p.second].second;
        log_file<< "U: "<< p.first << " V: "<< p.second <<endl;
        log_file<< "U: "<< U <<endl;
        log_file<< "V: "<< V << endl; 
        log_file<< endl;
    }
    log_file.close();
}

void createClusterEdgeList()
{
	cout << "createClusterEdgeList" << endl;
	int	blockTotal	= block_list.size();
	for (int i = 0; i < blockTotal; ++i)
	{
		if (block_list[i].size() > 0)
		{
			generateEdgilist(block_list[i]);
		}
	}
}

// find root of a point in components
int findRoot(int pointID, vector<int> &parentArr)
{
	if (parentArr[pointID] != pointID)
		parentArr[pointID]	= findRoot(parentArr[pointID], parentArr);

	return parentArr[pointID];
}

//  find clusters as connected components in a graph where edges are connection among records and vertices are record index
void findConnComp()
{
	int i, rootU, rootV, edgeTotal;
    vector<int> parentArr;

	for(i = 0; i < totalUniqueRecords; ++i)
	{
		parentArr.push_back(i);
	}

	edgeTotal	= set_of_edges.size();
	cout << "Number of edges " << edgeTotal  << endl;

    BOOST_FOREACH(pair p, set_of_edges) {
        rootU	= findRoot(p.first, parentArr);
		rootV	= findRoot(p.second, parentArr);

		if(rootU != rootV)
		{
            parentArr[rootV] = rootU;
		}
    }
    int componentsInClusters = 0;
    for(int i =0; i<parentArr.size(); i++) {
        int root;
        if (parentArr[i] == i) {
            // cout<< "i: " << i << " Parent: "<< parentArr[i] << " Val: " << uniqueRecords[i].first << " String "<< uniqueRecords[i].second<<endl;
            root = i;
        } else {
            root = findRoot(i, parentArr);
        }
        if (!approxConnectedComponents.count(root)) {
            vector<int> compononents;
            compononents.push_back(root);
            approxConnectedComponents[root] = compononents;
        } else {
            approxConnectedComponents[root].push_back(i);
            componentsInClusters++;
        }
    }
    cout<< "#Connected Components: " << approxConnectedComponents.size()<<endl;
    cout<< "#Total Non Root Nodes in graph: " << componentsInClusters << endl;
}

void printApproximateCluster() {
    int count = 0; 
    for (auto const& p : approxConnectedComponents) {
        for (int i=0; i<p.second.size(); i++) {
            cout<< uniqueRecords[p.second[i]].second << endl;
        }
        cout<< endl;
        count++;
        if (count>5) break;
    }
}

void findFinalConnectedComp(int intraCompDist) {
    int totalNodes = 0;
    int pairsAccessed = 0;
    for (auto const& p : approxConnectedComponents) {
        pairsAccessed++;
        int numComponents = p.second.size();
        totalNodes+=numComponents;
        bool distmat[numComponents][numComponents];
        vector<vector<string>> dataArr(numComponents); // to make cache-efficient, keep records in a row

        for(int c=0; c<p.second.size(); c++) {
            dataArr[c] = vec2D[uniqueRecords[p.second[c]].first];
        };

        for (int i =0; i<numComponents; i++) {
            distmat[i][i] = true;
            for (int j = i+1; j < numComponents; j++)
            {
                if (isLinkageOk(dataArr[i], dataArr[j], intraCompDist)) {
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

        //Printing
        // cout<< endl;
        // for(int i= 0; i< numComponents; i++) {
        //     for(int j= 0; j< numComponents; j++) {
        //         cout<< distmat[i][j] << "\t";
        //     }
        //     cout<< endl;
        // }

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

void printFinalConnectedClusters() {
    for(int i=0; i< 2; i++) {
        for(int j=0; j< finalConnectedComponents[i].size(); j++) {
            cout<< vec2D[uniqueRecords[finalConnectedComponents[i][j]].first][0] << endl;
            cout<< exactmatches[finalConnectedComponents[i][j]].size()<<endl;
            for(int k = 0; k<exactmatches[finalConnectedComponents[i][j]].size(); k++) {
                cout<< exactmatches[finalConnectedComponents[i][j]][k] <<endl;
            }
        }
        cout<< endl;
    }
}

void writeFinalConnectedComponentToFile(string& result_file_name) {
	ofstream out_file;
    out_file.open(result_file_name);

	for(int i = 0; i< finalConnectedComponents.size(); i++) {
        //cout<< "Cluster Size: "<< finalConnectedComponents[i].size() << endl;
		for(int j = 0; j< finalConnectedComponents[i].size(); j++) {
            //cout<< "Num Copies :" << exactmatches[finalConnectedComponents[i][j]].size() << endl;
            for(int k=0; k< exactmatches[finalConnectedComponents[i][j]].size(); k++) {
                //cout<< vec2D[exactmatches[finalConnectedComponents[i][j]][0]][0] << endl;
                out_file << vec2D[uniqueRecords[finalConnectedComponents[i][j]].first][0] << ",";
            }
		}
        out_file<< "\n";
	}
	out_file.close();
}


int main(int argc, char** argv) {
    // string filePath = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/";
    string filePath = "/home/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/";
    string fileName = argv[1];
    filePath = filePath + argv[1];
    getFormattedDataFromCSV(filePath);
	totalRecords = vec2D.size();
	getCombinedData();
    
    cout<< "Number of Records: " << vec2D.size() << endl;
	cout<< "Number of Combined data pairs: "<< combinedData.size() << endl;

    clock_t currTS1	= clock();
    radixSort(combinedData);
    double sorting_t	= (double)(clock() - currTS1) / CLOCKS_PER_SEC;
    cout<< "Sorting time "<< sorting_t << endl;
    // printSortedRecords();
    clock_t currTS3	= clock();
    getExactMatches();
    getUniqueEntries();
    double findingExact_t	= (double)(clock() - currTS3) / CLOCKS_PER_SEC;
    cout<< "My Exact Clustering time "<< findingExact_t << endl;

	clock_t currTS3_1	= clock();
	doSortedComp();
	double sortingComp_t	= (double)(clock() - currTS3_1) / CLOCKS_PER_SEC;
	cout<< "Sorting and Comparision time "<< sortingComp_t << endl;
    
    // for (size_t i = 0; i < uniqueRecords.size(); i++)
	// {
    //     cout<< uniqueRecords[i].first << " " << uniqueRecords[i].second <<endl;
    //     cout<< vec1D[uniqueRecords[i].first*attributes + 1] << endl;
	// }
    

    clock_t currTS4	= clock();
    doSuperBlocking();
    double blocking_t	= (double)(clock() - currTS4) / CLOCKS_PER_SEC;
    cout<< "Super Blocking time "<< blocking_t << endl;

    clock_t currTS5	= clock();
    createClusterEdgeList();
    double createEdges_t	= (double)(clock() - currTS5) / CLOCKS_PER_SEC;
    cout<< "Edge list creation Time "<< createEdges_t << endl;
    //printEdges();
    clock_t currTS6	= clock();
    findConnComp();
    double findComp_t	= (double)(clock() - currTS6) / CLOCKS_PER_SEC;
    cout<< "Connected Comp Find Time "<< findComp_t << endl;
    //printApproximateCluster();
    clock_t currTS7	= clock();
    findFinalConnectedComp(1);
    double findFinalComp_t	= (double)(clock() - currTS7) / CLOCKS_PER_SEC;
    cout<< "Final Connected Comps Find Time "<< findFinalComp_t << endl;
    //printFinalConnectedClusters();
    double total_t	= (double)(clock() - currTS1) / CLOCKS_PER_SEC;
    cout<< "Total run time "<< total_t << endl;

    //count number of pairs compared
	long long int total_comp = 0;
	for (size_t i = 0; i < block_list.size(); i++)
	{
		int cur_size = block_list[i].size();
		int cur_comp = (int)((cur_size * (cur_size - 1)) / 2);
		total_comp += cur_comp;
	}
	long long int tot_possible_com = (vec2D.size()*(vec2D.size() - 1))/2 ;
	
	cout << "Number of Possible comparison: " << tot_possible_com << endl;
	cout << "Total Comp. Done: "<< total_comp << endl;
	cout << "Reduction Ratio:" << ((long double)total_comp / (long double) tot_possible_com) << endl; 
	cout<< "Total Approximately Connected Components: " << approxConnectedComponents.size()<< endl;
	cout<< "Total Final Connected Components: " << finalConnectedComponents.size() << endl;

    // string out_file_path = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/data/";
    string out_file_path = "/home/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/data/";
	string out_name1 = out_file_path + "out_single_linkage_"+ fileName + "_super_blocking";
	string out_name2 = out_file_path + "out_complete_linkage_"+ fileName + "_super_blocking";
	string stat_file_name = "stat_"+ fileName + "_super_blocking";

	writeFinalConnectedComponentToFile(out_name2);

	string stat_file_path = out_file_path + stat_file_name;
    ofstream stat_file;
	stat_file.open(stat_file_path);
	stat_file << "DataSize: "<< vec2D.size() << endl;
	stat_file << "Number of Possible comparison: " << tot_possible_com << endl;
	stat_file << "Number of pairs compared: " << total_comp  << endl;
	stat_file << "Reduction Ratio:" << ((long double)total_comp / (long double) tot_possible_com) << endl;
	stat_file << "Number of Edges: "<< set_of_edges.size() << endl;
	stat_file << "Total Single Clusters: " << approxConnectedComponents.size()<< endl;
	stat_file << "Total Complete Clusters " << finalConnectedComponents.size() << endl;
	stat_file << "Total Time taken: " << total_t << " Seconds" << endl;
	stat_file << "Sorting time: " << sorting_t << " Seconds" << endl;
	stat_file << "Sorting with head removed and Comparison time: " << sortingComp_t << " Seconds" << endl;
	stat_file << "Finding Exact Clusters time: " << findingExact_t << " Seconds" << endl;
    stat_file << "Blocking time: " << blocking_t << " Seconds" << endl;
    stat_file << "Edge generation Time: " << createEdges_t << " Seconds" << endl;
    stat_file << "Find Approximate Clusters Time: " << findComp_t << " Seconds" << endl;
    stat_file << "Find Final Clusters Time: " << findFinalComp_t << " Seconds" << endl;
	stat_file.close();

    return 0;
}