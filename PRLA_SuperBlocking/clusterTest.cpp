#include <mpi.h>
#include <stdio.h>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <iostream>
#include <fstream>


using namespace std;


int main(int argc, char** argv) {
    MPI_Init(NULL, NULL);
    int rank;
    int world;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world);
    string filePath = "/home/job22011/RecordLinkage/data/ds1_50k";
    cout<< "World: " << world << endl;
    string line;
    vector<vector<string> > vec2D;
    ifstream records(filePath);
    int attributes;
    int ind = 0;
    while (getline (records, line)) {
        vector<string> result;
        boost::split(result, line, boost::is_any_of(","));
        vector<string> vec;
		for(int i=0; i<result.size(); i++) {
			auto last = std::remove_if(result[i].begin(), result[i].end(), [](auto ch) {
        								return ::ispunct(ch) || ::iswpunct(ch);
    								});
			result[i].erase(last, result[i].end()); //Remove junk left by remove_if() at the end of iterator
			boost::to_lower(result[i]);
			vec.push_back(result[i]);
		}
        vec2D.push_back(vec);
    }
    records.close();
    attributes = vec2D[0].size();
	cout<< "Attributes: "<<attributes << " rank: " << rank << endl;

    MPI_Finalize();
    return 0;
}

