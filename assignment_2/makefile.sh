g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -o pGenRLA_v2 pGenRLA_v2.cpp

./pGenRLA_v2 NC_VoterData_5M.csv

# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -o NC_voterData_analysis NC_voterData_analysis.cpp

# ./NC_voterData_analysis nc_test_5k.csv


# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -o prla_normalBlocking_unionFind prla_normalBlocking_unionFind.cpp

# ./prla_normalBlocking_unionFind simulated_records_original_1000_records.csv

# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -o prla_superBlocking prla_superBlocking.cpp

# ./prla_superBlocking ds1_50k

# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -o rla_newBlocking rla_newBlocking.cpp

# ./rla_newBlocking ds1_50k


# ./rla_newBlocking NC_voter_data_5M.csv

# ./rla_newBlocking ds1_50k

# ./rla_newBlocking ds2_100k

# ./rla_newBlocking ds3_200k

# ./rla_newBlocking ds4_400k

# ./rla_newBlocking ds5_600k

# ./rla_newBlocking ds6_800k

# ./rla_newParallelNormalBlocking ds7_1M