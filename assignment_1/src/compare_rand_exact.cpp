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
#include </opt/homebrew/Cellar/boost/1.79.0_2/include/boost/algorithm/string.hpp>

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
        vec.push_back(result[1]);
        vec.push_back(to_string(ind));
        vec.push_back(to_string(0));
        vec2D.push_back(vec);
        ind++;
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

void writeConnectedComponentToFile(vector<vector<int> > &connected_component_list, vector<vector<string> > &vec2D, string& result_file_name) {
	ofstream out_file;
    out_file.open(result_file_name);
	for(int i = 0; i< connected_component_list.size(); i++) {
		for(int j = 0; j< connected_component_list[i].size(); j++) {
			out_file << vec2D[connected_component_list[i][j]][0];
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

int main(int argc, char** argv)
{
	string file_name = argv[1];
    // Read Data
    //string file_path_1 = "/Users/joyanta/Documents/Research/Record\ Linkage/codes/my\ codes/19500671/ds1.1.1";
    //string file_path_2 = "/Users/joyanta/Documents/Research/Record\ Linkage/codes/my\ codes/19500671/ds1.1.2";
    string file_path = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/data/"+file_name;
	vector<vector<string> > vec2D = getFormattedDataFromCSV(file_path);

    cout<< "Vector Size: "<< vec2D.size()<<endl;
    // set seed
    srand(1);
	long long int tot_record = vec2D.size();
	long long int tot_pair = tot_record*(tot_record - 1);
	cout<< tot_record*(tot_record - 1) << " " << tot_pair << endl;
	vector<Edge> edges;
	set<tuple<int, int> > set_of_edges;

	for( int i = 0; i<tot_pair; i++) {
		int ind1 = rand() % tot_record;
		int ind2 = rand() % tot_record;
		while(ind1 == ind2) {
			ind2 = rand() % tot_record;
		}
		int edit_distance = calculateBasicED(vec2D[ind1][1], vec2D[ind2][1], 1);
		if (edit_distance <= 1) {
					
			tuple<int, int> edge_tuple;
			if (ind1<ind2) {
				edge_tuple = make_tuple(ind1,ind2);
			} else {
				edge_tuple = make_tuple(ind2,ind1);
			}

			if (!set_of_edges.count(edge_tuple)) {
				Edge edge(ind1, ind2);
				edges.push_back(edge);
				set_of_edges.insert(edge_tuple);
			}
		}
		if (i%100000 == 0)
		{
			cout<< "Number of pairs done" << i << endl;
		}	
	}
	cout<<" Number of Edges: "<< edges.size() << endl;
	cout<<" Number of Tuples: "<< set_of_edges.size() << endl;

}

// g++ -I /opt/homebrew/Cellar/boost/1.79.0_2/include -o data_processing data_processing.cpp 
