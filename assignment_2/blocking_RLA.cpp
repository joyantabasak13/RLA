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
#include </usr/local/boost_1_80_0/boost/algorithm/string.hpp>

using namespace std;
using namespace std::chrono;

int threshold = 99;

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

// NOT USED ANYMORE
// reads data from two tab deliminated datasets and merges them
// returns a vector of string vector which contains ssn & name in lowercase without whitespace 
vector<vector<string> > getData(string file_path_1, string file_path_2) {
    string line;
    vector<vector<string> > vec2D;

    // Read from the text file
    ifstream records_1(file_path_1);
    ifstream records_2(file_path_2);

    int ind = 0;

    while (getline (records_1, line)) {
        vector<string> result;
        boost::split(result, line, boost::is_any_of("\t"));
        vector<string> vec;
        string name = result[1] + result[2];
        boost::to_lower(name);
        vec.push_back(result[0]);
        vec.push_back(name);
        vec.push_back(to_string(ind));
        vec.push_back(to_string(0));
        vec2D.push_back(vec);
        ind++;
    }

    while (getline (records_2, line)) {
        vector<string> result;
        boost::split(result, line, boost::is_any_of(","));
        vector<string> vec;
        string name = result[1] + result[2];
        boost::to_lower(name);
        vec.push_back(result[0]);
        vec.push_back(name);
        vec.push_back(to_string(ind));
        vec.push_back(to_string(0));
        vec2D.push_back(vec);
        ind++;
    }
    records_1.close();
    records_2.close();
    return vec2D;
}

// reads data from one comma deliminated dataset (.CSV) file
// returns a vector of string vector which contains ssn & name 
vector<vector<string> > getFormattedDataFromCSV(string& file_path) {
    string line;
    vector<vector<string> > vec2D;
    ifstream records(file_path);

    int ind = 0;
    while (getline (records, line)) {
        vector<string> result;
        boost::split(result, line, boost::is_any_of(","));
        vector<string> vec;
        vec.push_back(result[0]);

		//Remove digits and alphanumerics from name
		auto last = std::remove_if(result[1].begin(), result[1].end(), [](auto ch) {
        								return ::isdigit(ch) || ::ispunct(ch);
    								});
		result[1].erase(last, result[1].end()); //Remove junk left by remove_if() at the end of iterator
        boost::to_lower(result[1]);
		vec.push_back(result[1]);
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

	void printConnectedComponents() {
		cout<<"Connected components of size "<< connected_component_list.size()<<endl;
		for (int i = 0; i < connected_component_list.size(); i++)
		{
			if(connected_component_list[i].size() > 0) {
				cout<<"has connected component"<<endl;
				    for (int v: connected_component_list[i]) {
				        cout << v << "--> ";
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


int calculateAccuracy(vector<vector<int> > &connected_component_list, vector<vector<string> > &vec2D) {
	int correct_component = 0;
	for(int i = 0; i< connected_component_list.size(); i++) {
		string ssn = vec2D[connected_component_list[i][0]][0];
		bool isSame = true;
		for(int j = 0; j< connected_component_list[i].size(); j++) {
			if (ssn != vec2D[connected_component_list[i][j]][0]) {
				isSame = false;
			}
		}
		if (isSame) {
			correct_component++;
		}
	}
	return correct_component;
}

void expand_connected_component_list(vector<vector<int> > &connected_component_list, vector<vector<string> > &exact_list, vector<vector<string> > &vec2D) {
	for(int i = 0; i< connected_component_list.size(); i++) {
		int current_list_size = connected_component_list[i].size();
		for(int j = 0; j< current_list_size; j++) {
			int start_ind = stoi(exact_list[connected_component_list[i][j]][2]);
			int end_ind = stoi(exact_list[connected_component_list[i][j]][3]);
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

//Input: vector of strings, size of the vector
//Output: size of longest string in the vector
size_t getMax(vector<vector<string> > &vec, int n){
    size_t max = vec[0][1].size();
    for (int i = 1; i < n; i++){
        if (vec[i][1].size()>max)
            max = vec[i][1].size();
    }
    return max;
}

// Input: vector of strings, size of the vector, digit
// Output: Vector laxically sorted by the digit-th character
// Note: Only handles Alphabets
void countSort(vector<vector<string> > &vec, int size, size_t k){
    string *b = NULL; int *c = NULL; string *d = NULL;
    b = new string[size];
    c = new int[27];
	d = new string[size];

    for (int i = 0; i <27; i++){
        c[i] = 0;
    }

    for (int j = 0; j <size; j++){   
        c[k < vec[j][1].size() ? (int)(unsigned char)vec[j][1][k] - 97 + 1 : 0]++;   //vec[j][1] is the name string
    }

    for (int f = 1; f <27; f++){
        c[f] += c[f - 1];
    }

    for (int r = size - 1; r >= 0; r--){
		int ind = k < vec[r][1].size() ? (int)(unsigned char)vec[r][1][k] - 97 + 1 : 0;
		if(ind < 0 || ind > 27) {
			cout<< "Naughty Character is: " << vec[r][1][k] << endl;
			cout<< "Naughty String is: "<< vec[r][1] << endl;
		}
        b[c[k < vec[r][1].size() ? (int)(unsigned char)vec[r][1][k] - 97 + 1 : 0] - 1] = vec[r][1];
        d[c[k < vec[r][1].size() ? (int)(unsigned char)vec[r][1][k] - 97 + 1 : 0] - 1] = vec[r][0];
		c[k < vec[r][1].size() ? (int)(unsigned char)vec[r][1][k] -97 + 1 : 0]--;
    }

    for (int l = 0; l < size; l++){
        vec[l][1] = b[l];
		vec[l][0] = d[l];
    }

    // avold memory leak
    delete[] b;
    delete[] c;
	delete[] d;
}

// Input: referance to a vector of strings, number of first r strings we want to sort (generally r = vec.size())
// Output: laxically sorted vector of strings
void radixSort(vector<vector<string> > &vec, int r){
    size_t max = getMax(vec, r);
	cout << "max is: "<< max << endl;
    for (size_t digit = max; digit > 0; digit--){ // size_t is unsigned, so avoid using digit >= 0, which is always true
        // cout<< "Digit: "<< digit << endl;
		countSort(vec, r, digit - 1);
    }

}

// Do exact clustering from lexically sorted vector of strings
void getExactMatches(vector<vector<string> > &vec, vector<vector<string> > &matches) {
	vector<string> tempVec;
	tempVec.push_back(vec[0][0]);
	tempVec.push_back(vec[0][1]);
	tempVec.push_back(to_string(0)); // starting index
	for (size_t i = 1; i < vec.size(); i++)
	{
		if(vec[i][1] != vec[i-1][1]) {
			tempVec.push_back(to_string(i-1)); //ending index
			matches.push_back(tempVec);
			tempVec.clear();
			tempVec.push_back(vec[i][0]);
			tempVec.push_back(vec[i][1]);
			tempVec.push_back(to_string(i));
		}
	}
	tempVec.push_back(to_string(vec.size()-1));
	matches.push_back(tempVec);
}

void doNormalBlocking(vector<vector<int> > &block_list, vector<vector<string> > &exactmatches) {
	int base = 26;
	int kmer = 3;
	int blockTotal = pow(26,kmer);
	// vector<vector<int> > block_list;
	block_list.resize(blockTotal);
	// int blocked_total = 0;
	//cout<< "Block resized to: "<< block_list.size() << endl;
	for (size_t i = 0; i < exactmatches.size(); i++)
	{
		// cout<< "String size: "<< exactmatches[i][1].size()<<endl;
		// blocked_total += exactmatches[i][1].size() - kmer + 1;
		for (size_t j = 0; j < exactmatches[i][1].size() - kmer + 1 ; j++)
		{
			int blockID = 0;
			for (size_t k = j; k < j+ kmer; k++)
			{
				// cout<< "Here k: " << k << endl;
				blockID += (int)(exactmatches[i][1][k] - 97) * pow(base,k-j);
			}
			// cout<<"i: "<< i << " j: "<< j << " BlockId "<< blockID << endl;
			block_list[blockID].push_back(i);
		}
	}
	// cout<< "Blocks assigned" << endl;
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

int main(int argc, char** argv)
{
	string file_name = argv[1];
    // Read Data
    //string file_path_1 = "/Users/joyanta/Documents/Research/Record\ Linkage/codes/my\ codes/19500671/ds1.1.1";
    //string file_path_2 = "/Users/joyanta/Documents/Research/Record\ Linkage/codes/my\ codes/19500671/ds1.1.2";
    string file_path = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/data/"+file_name;
	ofstream stat_file;
	vector<vector<string> > vec2D = getFormattedDataFromCSV(file_path);
	// for (size_t i = 0; i < vec2D.size(); i++)
	// {
	// 	cout<< vec2D[i][0] << " " << vec2D[i][1] << endl;
	// }
	
	// Start counting time
	auto start = high_resolution_clock::now();

	// Sort the Strings
	radixSort(vec2D, vec2D.size());

	auto sorting_time = high_resolution_clock::now();
	auto sorting_duration = duration_cast<std::chrono::milliseconds>(sorting_time - start);
	cout << "Time taken For sorting: " << sorting_duration.count()<< " milliSeconds" << endl;

    cout<< "after sorting:";
    for (size_t i = 0; i < vec2D.size(); i++) {
        cout<< vec2D[i][0] << " " << vec2D[i][1] << endl;
    }

	vector<vector<string> > exactmatches;

	getExactMatches(vec2D, exactmatches);

	cout << "Exact Clustering size: " << exactmatches.size() << endl;

	for (size_t i = 0; i < exactmatches.size(); i++)
	{
		cout<< exactmatches[i][0] << " " << exactmatches[i][1] << " " << exactmatches[i][2] << " " << exactmatches[i][3] << endl; 
	}

	// Do usual Blocking
	vector<vector<int> > block_list;
	doNormalBlocking(block_list, exactmatches);
	//doSuperBlocking(block_list, exactmatches);
	cout<< "Blocking Done" << endl;
	vector<Edge> edges;
	set<tuple<int, int> > set_of_edges;

	for( int i = 0; i< block_list.size(); i++) {
		for (int j = 0; j< block_list[i].size(); j++) {
			for (size_t k = j+1; k < block_list[i].size(); k++)
			{
				if (j!=k) {
					string first_element = exactmatches[block_list[i][j]][1];
					string second_element = exactmatches[block_list[i][k]][1];
					int edit_distance = calculateBasicED(first_element, second_element, 1);
					if (edit_distance <= 1) {
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
						// else {
						// 	cout << "i: "<<i<<" j: "<< j << " Exists!"<<endl;
						// }
					}
				}
			}	
		}
		// cout<< i<<" th Block Done!"<< endl;
	}

	for(auto e : set_of_edges) {
		int u = std::get<0>(e);
		int v = std::get<1>(e);
		// cout<< u << " " << v << endl;
		Edge edge(u, v);
		edges.push_back(edge);
	}
	
	cout<<" Number of Edges: "<< edges.size() << endl;
	Graph graph(edges, exactmatches.size());

	// printGraph(graph, exactmatches.size());

	graph.connectedComponents(exactmatches.size());

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<std::chrono::milliseconds>(stop - start);

	cout << "Time taken: " << duration.count()<< " milli seconds" << endl;

	//count number of pairs compared
	int total_comp = 0;
	for (size_t i = 0; i < block_list.size(); i++)
	{
		int cur_size = block_list[i].size();
		int cur_comp = (int)((cur_size * (cur_size - 1)) / 2);
		total_comp += cur_comp;
	}
	cout << "Total comp: "<< total_comp << endl;
	cout<<"Total Connected Components: " << graph.connected_component_list.size()<< endl;
	string out_name = "out_"+ file_name + "_super_blocking";
	string stat_file_name = "stat_"+ file_name + "_super_blocking";
	writeConnectedComponentToFile(graph.connected_component_list, exactmatches ,vec2D, out_name);
	string stat_file_path = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/data/"+stat_file_name;
	stat_file.open(stat_file_path);
	stat_file << "DataSize: "<< vec2D.size() << endl;
	stat_file << "Number of Possible comparison: " << (int)((vec2D.size()*(vec2D.size() - 1))/2)<< endl;
	stat_file << "Number of pairs compared: " << total_comp  << endl;
	stat_file << "Number of Edges: "<< edges.size() << endl;
	stat_file << "Total Connected Components: " << graph.connected_component_list.size()<< endl;
	stat_file << "Time taken: " << duration.count() << " miliSeconds" << endl;
	stat_file.close();
}

// g++ -std=c++17 -I /usr/local/boost_1_80_0 -o blocking_RLA blocking_RLA.cpp 
