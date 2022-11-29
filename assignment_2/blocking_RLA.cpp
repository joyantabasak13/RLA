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
// #include </usr/local/boost_1_80_0/boost/algorithm/string.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace std::chrono;

int threshold = 99;
int lenMax;
int totalRecords;
int totalUniqueRecords;

vector<int> blockfieldArr;
vector<string> uniqueblockfieldArr;
vector<vector<int> > exactmatches;
vector<vector<int> > block_list;
vector<vector<string> > vec2D;
vector<pair<int, string> > combinedData;
vector<pair<int, string> > uniqueRecords;

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

// Input: two vectors representing records, theshold
// Output: true if edit distance <= theshold, else false
bool isWithinThesholdEditDistance (vector<string> &vec1, vector<string> &vec2, int theshold) {
	bool isWithinTheshold = false;
	string first_element = vec1[1];
	string second_element = vec2[1];
	int edit_distance = 0;
	// Optimize more here in terms of vec-of-vec access
	int name_dist = calculateBasicED(first_element, second_element, theshold);
	if (name_dist <= theshold) {
		string first_dod = vec1[2];
		string second_dod = vec2[2];
		int dod_dist = calculateBasicED(first_dod, second_dod, theshold);
		if (dod_dist <= theshold) {
			string first_dob = vec1[3];
			string second_dob = vec2[3];
			int dob_dist = calculateBasicED(first_dob, second_dob, theshold);
			if (dob_dist <= theshold) {
				return true;
			}
		} 
	}
	return false;
}

// Combines the fields to a string and also keeps original vector index
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

void getBlockFieldData() {
	blockfieldArr.resize(vec2D.size());
	for (size_t i = 0; i < vec2D.size(); i++)
	{
		blockfieldArr[i] = vec2D[i][1].length();
	}
	
}
// reads data from one comma deliminated dataset (.CSV) file
// returns a vector of string vector which contains ssn, name, DoD, DoB 
vector<vector<string> > getFormattedDataFromCSV(string& file_path) {
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
    return vec2D;
}

// for graph datastucture
struct Edge {
    int src, dest;
    Edge(int my_src, int my_dest) {
        src = my_src;
        dest = my_dest;
    }
};

class Graph
{
public:

    vector<vector<int> > adjList;
    vector<vector<int> > connected_component_list;
	vector<vector<int> > final_cluster_list;
    Graph(vector<Edge> const &edges, int n)
    {
        adjList.resize(n);
		//undirected
        for (auto &edge: edges)
        {
            adjList[edge.src].push_back(edge.dest);
            adjList[edge.dest].push_back(edge.src);
        }
    }

	// Required for connected component extraction
    void DFSUtil(int vertex, vector<bool> &visited, vector<int> &connected_component)
    {
        // Mark the current node as visited
        visited[vertex] = true;
		connected_component.push_back(vertex);
        
        int i = 0;
		
		if (adjList[vertex].size() > 0) {
			for (i = 0; i < adjList[vertex].size(); i++) {
				if (visited[adjList[vertex][i]] == false) {
					DFSUtil(adjList[vertex][i], visited, connected_component);
				}
			}
		}
        
    }
	// finds connected components and keeps them in connected component lists
    void connectedComponents(long long int n)
    {
        // Mark all the vertices as not visited
        vector<bool> visited;
        for (int v = 0; v < n; v++)
            visited.push_back(false);

        for (int v = 0; v < n; v++) {
            if (visited[v] == false) {
                vector<int> connected_component;
                DFSUtil(v, visited, connected_component);
                connected_component_list.push_back(connected_component);
            }
        }
    }

    vector<vector<int> > getConnectedComponents() {
        return connected_component_list;
    }

	void generateFinalClusters(vector<vector<string> > &vec, int theshold) {
		for (int i = 0; i < connected_component_list.size(); i++)
		{
			bool mat[connected_component_list[i].size()][connected_component_list[i].size()];
			int arr[connected_component_list[i].size()];
			// Initialize
			for (int l = 0; l < connected_component_list[i].size(); l++) {
				arr[l] = l;
				for (int m = 0; m< connected_component_list[i].size(); m++) {
					mat[l][m] = false;
				}
			}
			for (int j = 0; j<connected_component_list[i].size(); j++) {
				for (int k = j + 1; k < connected_component_list[i].size(); k++) {
					bool isClose = isWithinThesholdEditDistance(vec[connected_component_list[i][j]], vec[connected_component_list[i][k]], theshold);
					mat[j][k] = isClose;
					mat[k][j] = isClose;
				}
			}
			// Extract paths
			int count = connected_component_list[i].size();
			
			vector<int> vec_path;
			for (size_t p = 0; p < connected_component_list[i].size(); p++)
			{
				if (arr[p] != -1) {
					arr[p] = -1;
					count--;
					vec_path.push_back(connected_component_list[i][p]);
					for (size_t q = 0; q< connected_component_list[i].size(); q++) {
						if (mat[p][q] == true) {
							arr[q] = -1;
							vec_path.push_back(connected_component_list[i][q]);
							count--;
						}
					}
					final_cluster_list.push_back(vec_path);
				}
			}
    	}
	}

	void printFinalConnectedComponents(vector<vector<string> > &vec) {
		cout<<"Final components of size "<< final_cluster_list.size()<<endl;
		for (int i = 0; i < final_cluster_list.size(); i++)
		{
			if(final_cluster_list[i].size() > 0) {
				//cout<<"has connected component"<<endl;
				for (int v: final_cluster_list[i]) {
				    cout << vec[v][0] << "--> ";
				}
				cout << endl;
			}
    	}
	}

	void printConnectedComponents(vector<vector<string> > &vec) {
		cout<<"Connected components of size "<< connected_component_list.size()<<endl;
		for (int i = 0; i < connected_component_list.size(); i++)
		{
			if(connected_component_list[i].size() > 0) {
				//cout<<"has connected component"<<endl;
				for (int v: connected_component_list[i]) {
				    cout << vec[v][0] << "--> ";
				}
				cout << endl;
			}
    	}
	}
};
 
void printGraph(Graph const &graph, int n)
{
    for (int i = 0; i < n; i++)
    {
        if(graph.adjList[i].size()>0) {
            cout << i << " ——> ";
            for (int v: graph.adjList[i]) {
                cout << v << ", ";
            }
        cout << endl;
        }
        
    }
}


void expand_connected_component_list(vector<vector<int> > &connected_component_list, vector<vector<string> > &exact_list, vector<vector<string> > &vec2D) {
	for(int i = 0; i< connected_component_list.size(); i++) {
		int current_list_size = connected_component_list[i].size();
		for(int j = 0; j< current_list_size; j++) {
			int start_ind = stoi(exact_list[connected_component_list[i][j]][4]);
			int end_ind = stoi(exact_list[connected_component_list[i][j]][5]);
			int extra_copies = end_ind - start_ind;
			for (size_t k = 0; k < extra_copies; k++)
			{
				connected_component_list[i].push_back(connected_component_list[i][j]);
			}
			// j = j + extra_copies;
		}
	}
}

void writeConnectedComponentToFile(vector<vector<int> > &connected_component_list, vector<vector<string> > &exact_list, vector<vector<string> > &vec2D, string& result_file_name) {
	ofstream out_file;
    out_file.open(result_file_name);
	expand_connected_component_list(connected_component_list, exact_list, vec2D );
	for(int i = 0; i< connected_component_list.size(); i++) {
		for(int j = 0; j< connected_component_list[i].size(); j++) {
			out_file << exact_list[connected_component_list[i][j]][0]; //here
			if (j == connected_component_list[i].size()-1)
			{
				out_file << "\n";
			} else {
				out_file << ",";
			}
			
		}
	}
	out_file.close();
}

// Input: vector of strings, NEEDs lenmax global var
// Output: Vector laxically sorted
void radixSort(vector<pair<int, string> > &combData){
	vector<pair<int, string>> tempArr(totalRecords);
	
	for (int i = lenMax - 1; i >= 0; --i) {
		vector<int> countArr(256, 0);
		
		for (int j = 0; j < totalRecords; ++j) {
			countArr[(combData[j].second)[i]]++;
		}
		
		for (int k = 1; k < 256; ++k)
			countArr[k]	+= countArr[k - 1];
		
		for (int j = totalRecords - 1; j >= 0; --j)
			tempArr[--countArr[(combData[j].second)[i]]]	= combData[j];
		
		for (int j = 0; j < totalRecords; ++j)
			combData[j]	= tempArr[j];
	}
}

// Do exact clustering from lexically sorted vector of int,string pair
void getExactMatches() {
	vector<int> tempVec;
	tempVec.push_back(combinedData[0].first);
	
	for (int i = 1; i < totalRecords; ++i) {
		if(combinedData[i].second.compare(combinedData[i - 1].second) == 0)
			tempVec.push_back(combinedData[i].first);
		else {
			exactmatches.push_back(tempVec);
			tempVec.clear();
			tempVec.push_back(combinedData[i].first);
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
    }
}

void getUniqueBlockFieldArr() {
	uniqueblockfieldArr.resize(totalUniqueRecords);
	for (size_t i = 0; i < totalUniqueRecords; i++)
	{
		uniqueblockfieldArr[i] = blockfieldArr[uniqueRecords[i].first];
	}
	
}


void doNormalBlocking() {
	int base = 26;
	int kmer = 3;
	//long long int total_blocked = 0;
	int blockTotal = pow(base,kmer);
	block_list.resize(blockTotal);

	// cout<< "Total unique clusters: " << totalUniqueRecords << endl;

	for (size_t i = 0; i < totalUniqueRecords; i++)
	{
		string blockingStr = vec2D[1][uniqueRecords[i].first];

		for (size_t j = 0; j < blockingStr.size() - kmer + 1 ; ++j)
		{ 
			int blockID = 0;

			for (size_t k = 0; k < kmer; ++k)
			{
				blockID += ((int)blockingStr[j+k] - 97) * pow(base,k);
			}
			block_list[blockID].push_back(i);
			//total_blocked++;
		}
	}
	//cout<< "Total blocked: " << total_blocked << endl;
}

void doSuperBlocking(vector<vector<int> > &block_list, vector<vector<string> > &exactmatches) {
	int base = 26;
	int kmer = 3;
	int blockTotal = pow(26,kmer); //number of blocks in a superblock
	// vector<vector<int> > block_list;
	block_list.resize(blockTotal*26);
	// int blocked_total = 0;

	for (size_t i = 0; i < exactmatches.size(); i++)
	{
		// blocked_total += exactmatches[i][1].size() - kmer + 1;
		for (size_t j = 0; j < exactmatches[i][1].size() - kmer + 1 ; j++)
		{
			int blockID = 0;
			for (size_t k = j; k < j+ kmer; k++)
			{
				blockID += (int)(exactmatches[i][1][k] - 97) * pow(base,k-j);
			}
			// cout<< blockID << endl;
			int superBlockId = exactmatches[i][1][0] - 97;
			block_list[blockID + superBlockId*blockTotal].push_back(i);
		}
	}
}

void doReverseSuperBlocking(vector<vector<int> > &block_list, vector<vector<string> > &exactmatches) {
	int base = 26;
	int kmer = 3;
	int blockTotal = pow(26,kmer); //number of blocks in a superblock
	block_list.resize(blockTotal*26);

	for (size_t i = 0; i < exactmatches.size(); i++)
	{
		// blocked_total += exactmatches[i][1].size() - kmer + 1;
		for (size_t j = exactmatches[i][1].size() - 1; j >= kmer - 1 ; j--)
		{
			int blockID = 0;
			for (size_t k = j; k > j - kmer; k--)
			{
				blockID += (int)(exactmatches[i][1][k] - 97) * pow(base,j-k);
			}
			// cout<< blockID << endl;
			int superBlockId = exactmatches[i][1][exactmatches[i][1].size() - 1] - 97;
			block_list[blockID + superBlockId*blockTotal].push_back(i);
		}
	}
}


int main(int argc, char** argv)
{
	string file_path = argv[1];
	ofstream stat_file;
	getFormattedDataFromCSV(file_path);
	totalRecords = vec2D.size();
	getCombinedData();
	getBlockFieldData();
	cout<< "Number of Records: " << vec2D.size() << endl;

	// for (size_t i = 0; i < 10; i++)
	// {
	// 	cout<< combinedData[i].second << endl;
	// }

	// Start counting time
	clock_t start_t = clock();

	// Sort the Strings
	radixSort(combinedData);
	cout<< "Sorting Done" << endl;

	clock_t sorting_t = clock();
	double sorting_duration_t	= (double)(clock() - start_t) / CLOCKS_PER_SEC;
    cout<< "Sorting time: "<< sorting_duration_t <<" seconds" << endl;

    // for (size_t i = 0; i < 10; i++) {
    //     cout<< combinedData[i].second << endl;
    // }
	clock_t exact_clustering_start_t = clock();
	getExactMatches();
	getUniqueEntries();
	
	double exact_clustering_duration_t	= (double)(clock() - exact_clustering_start_t) / CLOCKS_PER_SEC;
    cout<< "Exact Clustering time: "<< exact_clustering_duration_t <<" seconds" << endl;


	// for (size_t i = 0; i < exactmatches.size(); i++)
	// {
	// 	if ((exactmatches[i].size() != 4) && (exactmatches[i].size()!=1)) {
	// 		for(int j=0; j<exactmatches[i].size(); j++) {
	// 			cout<< vec2D[exactmatches[i][j]][0] << " " << vec2D[exactmatches[i][j]][1] << "  ";
	// 		}
	// 		cout<< endl;
	// 	}
	// }

	// Do usual Blocking
	clock_t blocking_start_t = clock();
	doNormalBlocking();
	double blocking_duration_t	= (double)(clock() - blocking_start_t) / CLOCKS_PER_SEC;
    cout<< "Normal Blocking time: "<< blocking_duration_t <<" seconds" << endl;

	// doSuperBlocking(block_list, exactmatches);
	// doReverseSuperBlocking(block_list, exactmatches);
	cout<< "Blocking Done" << endl;
	vector<Edge> edges;
	set<tuple<int, int> > set_of_edges;

	for( int i = 0; i< block_list.size(); i++) {
		for (int j = 0; j< block_list[i].size(); j++) {
			for (size_t k = j+1; k < block_list[i].size(); k++)
			{
				string first_element = uniqueRecords[block_list[i][j]].second;
				string second_element = uniqueRecords[block_list[i][k]].second;
				// Optimize more here in terms of vec-of-vec access
				int name_dist = calculateBasicED(first_element, second_element, 1);
				// if (name_dist <= 1) {
				// 	string first_dod = exactmatches[block_list[i][j]][2];
				// 	string second_dod = exactmatches[block_list[i][k]][2];
				// 	int dod_dist = calculateBasicED(first_dod, second_dod, 1);
				// 	if (dod_dist <= 1) {
				// 		string first_dob = exactmatches[block_list[i][j]][3];
				// 		string second_dob = exactmatches[block_list[i][k]][3];
				// 		int dob_dist = calculateBasicED(first_dob, second_dob, 1);
				// 		if (dob_dist > 1) {
				// 			edit_distance = threshold + 1; //Maybe you don't need to do this, dist func already returns this
				// 		}
				// 	} else {
				// 		edit_distance = threshold + 1;
				// 	}
				// } else {
				// 	edit_distance = threshold + 1;
				// }
				if (name_dist <= 1) {
					//cout<< "No Prob. i = "<< i << " j= "<< j << " k= "<< k <<endl;
					tuple<int, int> edge_tuple;
					int j_th_record_id = block_list[i][j];
					int k_th_record_id = block_list[i][k];
					if (j_th_record_id < k_th_record_id) {
						edge_tuple = make_tuple(j_th_record_id, k_th_record_id);
					} else {
						edge_tuple = make_tuple(k_th_record_id, j_th_record_id);
					}

					if (!set_of_edges.count(edge_tuple)) {
						// Edge edge(i, j);
						// edges.push_back(edge);
						set_of_edges.insert(edge_tuple);
					} 

				}
				
			}	
		}

	}

	for(auto e : set_of_edges) {
		int u = std::get<0>(e);
		int v = std::get<1>(e);
		// cout<< u << " " << v << endl;
		Edge edge(u, v);
		edges.push_back(edge);
	}
	
	// cout<<" Number of Edges: "<< edges.size() << endl;
	Graph graph(edges, totalUniqueRecords);

	// // printGraph(graph, exactmatches.size());

	graph.connectedComponents(totalUniqueRecords);
	// //graph.printConnectedComponents(exactmatches);

	double approx_clustering_duration_t	= (double)(clock() - start_t) / CLOCKS_PER_SEC;
    cout<< "Approx Clustering time: "<< approx_clustering_duration_t <<" seconds" << endl;
	return 0;
	// cout << "Time taken: (For single clustering) " << duration_single.count()<< " milli seconds" << endl;

	// // Get final Clusters
	// graph.generateFinalClusters(exactmatches, 1);
	// //graph.printFinalConnectedComponents(exactmatches);

	// auto stop = high_resolution_clock::now();
	// auto duration = duration_cast<std::chrono::milliseconds>(stop - start);

	// cout << "Time taken: " << duration.count()<< " milli seconds" << endl;

	// //count number of pairs compared
	// int total_comp = 0;
	// for (size_t i = 0; i < block_list.size(); i++)
	// {
	// 	int cur_size = block_list[i].size();
	// 	int cur_comp = (int)((cur_size * (cur_size - 1)) / 2);
	// 	total_comp += cur_comp;
	// }
	// long long int tot_possible_com = (int)((vec2D.size()*(vec2D.size() - 1))/2);
	
	// cout << "Number of Possible comparison: " << (int)((vec2D.size()*(vec2D.size() - 1))/2)<< endl;
	// cout << "Total comp: "<< total_comp << endl;
	// cout << "Reduction Ratio:" << double(total_comp)/double(tot_possible_com) << endl; 
	// cout<< "Total Approximately Connected Components: " << graph.connected_component_list.size()<< endl;
	// cout<< "Total Final Connected Components: " << graph.final_cluster_list.size() << endl;
	// string out_name1 = "out_single_linkage_"+ file_name + "_normal_blocking";
	// string out_name2 = "out_complete_linkage_"+ file_name + "_normal_blocking";
	// string stat_file_name = "stat_"+ file_name + "_normal_blocking";
	// writeConnectedComponentToFile(graph.connected_component_list, exactmatches ,vec2D, out_name1);
	// writeConnectedComponentToFile(graph.final_cluster_list, exactmatches ,vec2D, out_name2);
	// // string stat_file_path = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/data/"+stat_file_name;
	// stat_file.open(stat_file_name);
	// stat_file << "DataSize: "<< vec2D.size() << endl;
	// stat_file << "Number of Possible comparison: " << (int)((vec2D.size()*(vec2D.size() - 1))/2)<< endl;
	// stat_file << "Number of pairs compared: " << total_comp  << endl;
	// stat_file << "Reduction Ratio:" << float(total_comp)/float(tot_possible_com) << endl;
	// stat_file << "Number of Edges: "<< edges.size() << endl;
	// stat_file << "Total Single Clusters: " << graph.connected_component_list.size()<< endl;
	// stat_file << "Total Conplete Clusters " << graph.final_cluster_list.size() << endl;
	// stat_file << "Total Time taken: " << duration.count() << " miliSeconds" << endl;
	// stat_file << "Single Linkage Time Taken: " << duration_single.count() << " milliSeconds" << endl;
	// stat_file.close();
}

// g++ -std=c++17 -I boost_1_51_0 -o blocking_RLA blocking_RLA.cpp 

// g++ -std=c++17 -I /usr/local/boost_1_80_0 -o blocking_RLA blocking_RLA.cpp

// ./blocking_RLA ds7_1M
