# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -o blocking_RLA blocking_RLA.cpp

# ./rla_newParallelNormalBlocking /Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/ds1_50k

g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -o rla_newParallelNormalBlocking rla_newParallelNormalBlocking.cpp

# ./rla_newParallelNormalBlocking ds1_50k

# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -o rla_newBlocking rla_newBlocking.cpp

# ./rla_newBlocking ds1_50k

# ./rla_newBlocking ds2_100k

# ./rla_newBlocking ds3_200k

# ./rla_newBlocking ds4_400k

# ./rla_newBlocking ds5_600k

# ./rla_newBlocking ds6_800k

./rla_newParallelNormalBlocking ds7_1M