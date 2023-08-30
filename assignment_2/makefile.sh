# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -I /usr/local/  -o NC_voterData_analysis_threshold_jaroDist NC_voterData_analysis_threshold_jaroDist.cpp
# ./NC_voterData_analysis_threshold_jaroDist NC_voterData_5M_Source_Annotated.csv
# ./NC_voterData_analysis_threshold_jaroDist NC_5Percent

# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -I /usr/local/  -o NC_voterData_analysis_threshold NC_voterData_analysis_threshold.cpp
# ./NC_voterData_analysis_threshold_jaroDist NC_voterData_5M_Source_Annotated.csv
# ./NC_voterData_analysis_threshold NC_5Percent

g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -I /usr/local/  -o linearJW linearJW.cpp

./linearJW f_name_alpha.csv
# ./linearJW fullName_alpha.csv
# ./linearJW fullName_suburb_alpha.csv
# ./linearJW strSize10_alpha.csv
# ./linearJW strSize26_alpha.csv
# ./linearJW strSize36.csv

# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -I /usr/local/  -o ukkonenLDEditDistance ukkonenLDEditDistance.cpp
# ./ukkonenLDEditDistance

# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -o pGenRLA_v6_NC_jaro pGenRLA_v6_NC_jaro.cpp
# ./pGenRLA_v6_NC_jaro nc_test_5k.csv
# ./pGenRLA_v6_NC_jaro NC_voterData_5M_Source_Annotated.csv

# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -o pGenRLA_v5_NC pGenRLA_v5_NC.cpp
# ./pGenRLA_v5_NC NC_voterData_5M_Source_Annotated.csv


# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -I /usr/local/  -o extractCompleteLinkage_v5 extractCompleteLinkage_v5.cpp
# ./extractCompleteLinkage_v5 nc_TEST.csv ncTEST.csv
# ./extractCompleteLinkage_v5 NC_voterData_5M_Source_Annotated.csv ALL_RECID_SingleLinkage_jaroWinkler_80

# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -I /usr/local/  -o extractCompleteLinkage_v6 extractCompleteLinkage_v6.cpp
# ./extractCompleteLinkage_v6 NC_voterData_5M_Source_Annotated.csv ALL_RECID_SingleLinkage_NC_voterData_5M_Source_Annotated.csv_Lev

# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -o pGenRLA_v4_NC pGenRLA_v4_NC.cpp

# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -o pGenRLA_v5_NC pGenRLA_v5_NC.cpp
# ./pGenRLA_v5_NC nc_test_5k.csv
# ./pGenRLA_v5_NC nc_TEST.csv
# ./pGenRLA_v5_NC NC_voterData_5M_Source_Annotated.csv



# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -o NC_voterData_analysis_BlockingKey NC_voterData_analysis_BlockingKey.cpp
# ./NC_voterData_analysis_BlockingKey NC_1Percent
# ./NC_voterData_analysis_BlockingKey NC_5Percent
# ./NC_voterData_analysis_BlockingKey NC_10Percent
# ./NC_voterData_analysis_BlockingKey nc_test_5k.csv
#  ./NC_voterData_analysis_BlockingKey NC_voterData_5M_Source_Annotated.csv


# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -o pGenRLA_v3_NC pGenRLA_v3_NC.cpp
# ./pGenRLA_v3_NC nc_test_5k.csv
# ./pGenRLA_v3_NC NC_voterData_5M_Source_Annotated.csv


# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -o prla_normalBlocking_v2 prla_normalBlocking_v2.cpp
# ./prla_normalBlocking_v2 ds7_1M

# ./prla_normalBlocking_v2 simulated_records_UnionFindTest_80000_records.csv

# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -o prla_superBlocking_v3 prla_superBlocking_v3.cpp
# ./prla_superBlocking_v3 ds7_1M

# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -o rla_newBlocking rla_newBlocking.cpp
# ./rla_newBlocking ds1_50k
# ./rla_newBlocking NC_voter_data_5M.csv

# ./rla_newBlocking ds1_50k
# ./rla_newBlocking ds2_100k
# ./rla_newBlocking ds3_200k
# ./rla_newBlocking ds4_400k
# ./rla_newBlocking ds5_600k
# ./rla_newBlocking ds6_800k
# ./rla_newBlocking ds7_1M