# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -o rla_newParallelNormalBlocking rla_newParallelNormalBlocking.cpp

# ./rla_newParallelNormalBlocking nc_test_5k.csv

g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -o rla_newParallelNormalBlocking rla_newParallelNormalBlocking.cpp

./rla_newParallelNormalBlocking ds7_1M

# ./rla_newBlocking NC_voter_data_5M.csv

# ./rla_newBlocking ds1_50k

# ./rla_newBlocking ds2_100k

# ./rla_newBlocking ds3_200k

# ./rla_newBlocking ds4_400k

# ./rla_newBlocking ds5_600k

# ./rla_newBlocking ds6_800k

# ./rla_newParallelNormalBlocking ds7_1M