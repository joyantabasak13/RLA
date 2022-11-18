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
// #include </usr/local/boost_1_80_0/boost/algorithm/string.hpp>
#include <boost/algorithm/string.hpp>

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

// reads data from one comma deliminated dataset (.CSV) file
// returns a vector of string vector which contains ssn, name, DoD, DoB 
vector<vector<string> > getFormattedDataFromCSV(string& file_path) {
    string line;
    vector<vector<string> > vec2D;
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
void countSort(vector<vector<string> > &vec, int size, size_t k, int col){
    int *c = NULL; 
	string *name = NULL; 
	string *ssn = NULL;
	string *dod = NULL;
	string *dob = NULL;
	int alphabet_size = 0;
	int mode = 0; // Mode = 1 for english alphabets, Mode = 2 for numericals(date)
	if (col == 1) {
		alphabet_size = 27;
		mode = 1;
	} else {
		alphabet_size = 11;
		mode = 2;
	}
	c = new int[alphabet_size];
    name = new string[size];
	ssn = new string[size];
	dod = new string[size];
	dob = new string[size];

    for (int i = 0; i <alphabet_size; i++){
        c[i] = 0;
    }
	
    for (int j = 0; j <size; j++){
		if (mode == 1) {
			c[k < vec[j][col].size() ? (int)(unsigned char)vec[j][col][k] - 97 + 1 : 0]++;   //vec[j][1] is the name string
		} else if (mode == 2) {
			c[k < vec[j][col].size() ? (int)(unsigned char)vec[j][col][k] - 48 + 1 : 0]++;	// For dates
		}
    }

    for (int f = 1; f <alphabet_size; f++){
        c[f] += c[f - 1];
    }
	
    for (int r = size - 1; r >= 0; r--){
		int ind = 0;
		if (mode == 1) {
			ind = k < vec[r][col].size() ? (int)(unsigned char)vec[r][col][k] - 97 + 1 : 0;
		} else if (mode == 2) {
			ind = k < vec[r][col].size() ? (int)(unsigned char)vec[r][col][k] - 48 + 1 : 0;
		}
		// Only fire when erroreneous character is present
		// if(ind < 0 || ind > alphabet_size) {
		// 	cout<< "Naughty Character is: " << vec[r][col][k] << endl;
		// 	cout<< "Naughty String is: "<< vec[r][col] << endl;
		// }
		// cout<< "Accessing out-of-array? " << ind << endl;
		// TO-DO : Need a better mechnism than keeping 5 arrays. Keep a vec of vec instead ?
		if (mode == 1) {
			name[c[k < vec[r][col].size() ? (int)(unsigned char)vec[r][col][k] - 97 + 1 : 0] - 1] = vec[r][1];
        	ssn[c[k < vec[r][col].size() ? (int)(unsigned char)vec[r][col][k] - 97 + 1 : 0] - 1] = vec[r][0];
			dod[c[k < vec[r][col].size() ? (int)(unsigned char)vec[r][col][k] - 97 + 1 : 0] - 1] = vec[r][2];
			dob[c[k < vec[r][col].size() ? (int)(unsigned char)vec[r][col][k] - 97 + 1 : 0] - 1] = vec[r][3];
			c[k < vec[r][col].size() ? (int)(unsigned char)vec[r][col][k] - 97 + 1 : 0]--;
		} else if (mode == 2) {
			name[c[k < vec[r][col].size() ? (int)(unsigned char)vec[r][col][k] - 48 + 1 : 0] - 1] = vec[r][1];
        	ssn[c[k < vec[r][col].size() ? (int)(unsigned char)vec[r][col][k] - 48 + 1 : 0] - 1] = vec[r][0];
			dod[c[k < vec[r][col].size() ? (int)(unsigned char)vec[r][col][k] - 48 + 1 : 0] - 1] = vec[r][2];
			dob[c[k < vec[r][col].size() ? (int)(unsigned char)vec[r][col][k] - 48 + 1 : 0] - 1] = vec[r][3];
			c[k < vec[r][col].size() ? (int)(unsigned char)vec[r][col][k] - 48 + 1 : 0]--;
		}
    }

    for (int l = 0; l < size; l++){
		vec[l][0] = ssn[l];
        vec[l][1] = name[l];
		vec[l][2] = dod[l];
		vec[l][3] = dob[l];
		
    }

    // avold memory leak
	delete[] c;
    delete[] ssn;
	delete[] name;
	delete[] dod;
	delete[] dob;
}


// Input: referance to a vector of strings, number of first r strings we want to sort (generally r = vec.size()), column of vec we want to sort
// Output: laxically sorted vector of strings
void radixSort(vector<vector<string> > &vec, int r, int col){
	size_t max = 0;
	if (col == 1) {
		// cout<<"getting max" << endl;
		max = getMax(vec, r);
		// cout<< "got max" << endl;
	} else {
		max = 8;	// Dates are all 8 digit numbers
	}			
	cout << "max is: "<< max << endl;
    for (size_t digit = max; digit > 0; digit--){ // size_t is unsigned, so avoid using digit >= 0, which is always true
        // cout<< "Digit: "<< digit << endl;
		countSort(vec, r, digit - 1, col);
    }
}

// Do exact clustering from lexically sorted vector of strings
void getExactMatches(vector<vector<string> > &vec, vector<vector<string> > &matches) {
	vector<string> tempVec;
	tempVec.push_back(vec[0][0]);
	tempVec.push_back(vec[0][1]);
	tempVec.push_back(vec[0][2]);
	tempVec.push_back(vec[0][3]);
	tempVec.push_back(to_string(0)); // starting index
	for (size_t i = 1; i < vec.size(); i++)
	{
		if ((vec[i][1] != vec[i-1][1])|| (vec[i][2] != vec[i-1][2]) || (vec[i][3] != vec[i-1][3])) {
			tempVec.push_back(to_string(i-1)); //ending index
			matches.push_back(tempVec);
			tempVec.clear();
			tempVec.push_back(vec[i][0]);
			tempVec.push_back(vec[i][1]);
			tempVec.push_back(vec[i][2]);
			tempVec.push_back(vec[i][3]);
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
	string file_name = argv[1];
    // Read Data
	// For server
    // string file_path = "/home/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/"+file_name;
	// For my laptop
	string file_path = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/"+file_name;
	ofstream stat_file;
	vector<vector<string> > vec2D = getFormattedDataFromCSV(file_path);
	cout<< "Number of Records: " << vec2D.size() << endl;
	// for (size_t i = 0; i < vec2D.size(); i++)
	// {
	// 	cout<< vec2D[i][0] << " " << vec2D[i][1] << " " << vec2D[i][2] << " " << vec2D[i][3] << endl;
	// }

	// Start counting time
	auto start = high_resolution_clock::now();

	// Sort the Strings
	radixSort(vec2D, vec2D.size(), 3); // By Date of Birth
	cout<< "DOB Done" << endl;
	radixSort(vec2D, vec2D.size(), 2); // By Date of Death
	cout<< "DoD done" << endl;
	radixSort(vec2D, vec2D.size(), 1); // By name 
	cout<< "Name Done" << endl;

	auto sorting_time = high_resolution_clock::now();
	auto sorting_duration = duration_cast<std::chrono::milliseconds>(sorting_time - start);
	cout << "Time taken For sorting: " << sorting_duration.count()<< " milliSeconds" << endl;

    cout<< "after sorting:";
    // // for (size_t i = 0; i < vec2D.size(); i++) {
    // //     cout<< vec2D[i][0] << " " << vec2D[i][1] << " " << vec2D[i][2] << " " << vec2D[i][3] << endl;
    // // }

	// vector<vector<string> > exactmatches;

	// getExactMatches(vec2D, exactmatches);

	// cout << "Exact Clustering size: " << exactmatches.size() << endl;

	// // for (size_t i = 0; i < exactmatches.size(); i++)
	// // {
	// // 	cout<< exactmatches[i][0] << " " << exactmatches[i][1] << " " << exactmatches[i][2] << " " << exactmatches[i][3] << " " << exactmatches[i][4] << " " << exactmatches[i][5] << endl; 
	// // }

	// // for (size_t i = 0; i < exactmatches.size(); i++)
	// // {
	// // 	if((exactmatches[i][4] != exactmatches[i][5]) && ((stoi(exactmatches[i][5]) - stoi(exactmatches[i][4])) != 3)) {
	// // 		cout<< exactmatches[i][0] << " " << exactmatches[i][1] << " " << exactmatches[i][2] << " " << exactmatches[i][3] << " " << exactmatches[i][4] << " " << exactmatches[i][5] << endl; 
	// // 	}	
	// // }

	// // Do usual Blocking
	// vector<vector<int> > block_list;
	// doNormalBlocking(block_list, exactmatches);
	// // doSuperBlocking(block_list, exactmatches);
	// // doReverseSuperBlocking(block_list, exactmatches);
	// cout<< "Blocking Done" << endl;
	// vector<Edge> edges;
	// set<tuple<int, int> > set_of_edges;

	// for( int i = 0; i< block_list.size(); i++) {
	// 	for (int j = 0; j< block_list[i].size(); j++) {
	// 		for (size_t k = j+1; k < block_list[i].size(); k++)
	// 		{
	// 			if (j!=k) {
	// 				string first_element = exactmatches[block_list[i][j]][1];
	// 				string second_element = exactmatches[block_list[i][k]][1];
	// 				int edit_distance = 0;
	// 				// Optimize more here in terms of vec-of-vec access
	// 				int name_dist = calculateBasicED(first_element, second_element, 1);
	// 				if (name_dist <= 1) {
	// 					string first_dod = exactmatches[block_list[i][j]][2];
	// 					string second_dod = exactmatches[block_list[i][k]][2];
	// 					int dod_dist = calculateBasicED(first_dod, second_dod, 1);
	// 					if (dod_dist <= 1) {
	// 						string first_dob = exactmatches[block_list[i][j]][3];
	// 						string second_dob = exactmatches[block_list[i][k]][3];
	// 						int dob_dist = calculateBasicED(first_dob, second_dob, 1);
	// 						if (dob_dist > 1) {
	// 							edit_distance = threshold + 1; //Maybe you don't need to do this, dist func already returns this
	// 						}
	// 					} else {
	// 						edit_distance = threshold + 1;
	// 					}
	// 				} else {
	// 					edit_distance = threshold + 1;
	// 				}
	// 				if (edit_distance <= 1) {
	// 					//cout<< "No Prob. i = "<< i << " j= "<< j << " k= "<< k <<endl;
	// 					tuple<int, int> edge_tuple;
	// 					int j_th_record_id = block_list[i][j];
	// 					int k_th_record_id = block_list[i][k];
	// 					if (j_th_record_id < k_th_record_id) {
	// 						edge_tuple = make_tuple(j_th_record_id, k_th_record_id);
	// 					} else {
	// 						edge_tuple = make_tuple(k_th_record_id, j_th_record_id);
	// 					}

	// 					if (!set_of_edges.count(edge_tuple)) {
	// 						// Edge edge(i, j);
	// 						// edges.push_back(edge);
	// 						set_of_edges.insert(edge_tuple);
	// 					} 
	// 					// else {
	// 					// 	cout << "i: "<<i<<" j: "<< j << " Exists!"<<endl;
	// 					// }
	// 				}
	// 			}
	// 		}	
	// 	}
	// 	// cout<< i<<" th Block Done!"<< endl;
	// }

	// for(auto e : set_of_edges) {
	// 	int u = std::get<0>(e);
	// 	int v = std::get<1>(e);
	// 	// cout<< u << " " << v << endl;
	// 	Edge edge(u, v);
	// 	edges.push_back(edge);
	// }
	
	// cout<<" Number of Edges: "<< edges.size() << endl;
	// Graph graph(edges, exactmatches.size());

	// // printGraph(graph, exactmatches.size());

	// graph.connectedComponents(exactmatches.size());
	// //graph.printConnectedComponents(exactmatches);
	// auto stop_single = high_resolution_clock::now();
	// auto duration_single = duration_cast<std::chrono::milliseconds>(stop_single - start);

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
