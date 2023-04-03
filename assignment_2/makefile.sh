# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -I /usr/local/  -o single_to_completeLinkage_analysis single_to_completeLinkage_analysis.cpp

# ./single_to_completeLinkage_analysis NC_voterData_5M_Source_Annotated.csv

g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -o pGenRLA_v4_NC pGenRLA_v4_NC.cpp

# ./pGenRLA_v4_NC nc_test_5k.csv

# ./pGenRLA_v4_NC nc_TEST.csv

./pGenRLA_v4_NC NC_voterData_5M_Source_Annotated.csv


# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -o NC_voterData_analysis_threshold NC_voterData_analysis_threshold.cpp

# ./NC_voterData_analysis_threshold nc_test_5k.csv

# ./NC_voterData_analysis_threshold NC_voterData_5M_Source_Annotated.csv



# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -o NC_voterData_analysis_BlockingKey NC_voterData_analysis_BlockingKey.cpp

# ./NC_voterData_analysis_BlockingKey nc_test_5k.csv

# ./NC_voterData_analysis_BlockingKey NC_voterData_5M_Source_Annotated.csv


# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -o pGenRLA_v3_NC pGenRLA_v3_NC.cpp

# # ./pGenRLA_v3_NC nc_test_5k.csv

# ./pGenRLA_v3_NC NC_voterData_5M_Source_Annotated.csv

# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -o pGenRLA_v2_NC pGenRLA_v2_NC.cpp

# ./pGenRLA_v2_NC NC_voterData_5M_Source_Annotated.csv

# ./pGenRLA_v2 nc_test_5k.csv

# ./pGenRLA_v2 NC_VoterData_5M.csv

# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -I /usr/local/ -o NC_voterData_analysis NC_voterData_analysis.cpp

# ./NC_voterData_analysis NC_voterData_5M_Source_Annotated.csv

# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -o prla_normalBlocking_v2 prla_normalBlocking_v2.cpp

# ./prla_normalBlocking_v2 simulated_records_UnionFindTest_80000_records.csv

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