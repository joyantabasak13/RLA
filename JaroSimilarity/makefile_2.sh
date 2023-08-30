# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -I /usr/local/  -o eColiGenomeComp eColiGenomeComp.cpp
# ./eColiGenomeComp eColi_genes_300.txt
# ./eColiGenomeComp eColi_genes_500.txt
# ./eColiGenomeComp eColi_genes_1000.txt

#./makefile.sh |& tee -a output_set_2.txt

g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -I /usr/local/  -o randomStringComp randomStringComp.cpp

./randomStringComp randomString_25_alpha_5_set_2.txt
./randomStringComp randomString_25_alpha_10_set_2.txt
./randomStringComp randomString_25_alpha_15_set_2.txt
./randomStringComp randomString_25_alpha_20_set_2.txt
./randomStringComp randomString_25_alpha_26_set_2.txt

./randomStringComp randomString_50_alpha_5_set_2.txt
./randomStringComp randomString_50_alpha_10_set_2.txt
./randomStringComp randomString_50_alpha_15_set_2.txt
./randomStringComp randomString_50_alpha_20_set_2.txt
./randomStringComp randomString_50_alpha_26_set_2.txt

./randomStringComp randomString_75_alpha_5_set_2.txt
./randomStringComp randomString_75_alpha_10_set_2.txt
./randomStringComp randomString_75_alpha_15_set_2.txt
./randomStringComp randomString_75_alpha_20_set_2.txt
./randomStringComp randomString_75_alpha_26_set_2.txt

./randomStringComp randomString_100_alpha_5_set_2.txt
./randomStringComp randomString_100_alpha_10_set_2.txt
./randomStringComp randomString_100_alpha_15_set_2.txt
./randomStringComp randomString_100_alpha_20_set_2.txt
./randomStringComp randomString_100_alpha_26_set_2.txt

./randomStringComp randomString_125_alpha_5_set_2.txt
./randomStringComp randomString_125_alpha_10_set_2.txt
./randomStringComp randomString_125_alpha_15_set_2.txt
./randomStringComp randomString_125_alpha_20_set_2.txt
./randomStringComp randomString_125_alpha_26_set_2.txt

./randomStringComp randomString_150_alpha_5_set_2.txt
./randomStringComp randomString_150_alpha_10_set_2.txt
./randomStringComp randomString_150_alpha_15_set_2.txt
./randomStringComp randomString_150_alpha_20_set_2.txt
./randomStringComp randomString_150_alpha_26_set_2.txt


# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -I /usr/local/  -o prla_nB_JaroLinear prla_nB_JaroLinear.cpp
# ./prla_nB_JaroLinear ds7_1M

# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -I /usr/local/  -o prla_NC_jaro_linear prla_NC_jaro_linear.cpp
# ./prla_NC_jaro_linear NC_voterData_5M_Source_Annotated.csv

# g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -I /usr/local/  -o extractCL_JaroQuad extractCL_JaroQuad.cpp
# ./extractCL_JaroQuad NC_voterData_5M_Source_Annotated.csv ALL_RECID_SingleLinkage_quad_jaro_




# Datasets:

# Alphanuemeric Strings
# ./executable f_name_alpha.csv
# ./executable fullName_alpha.csv
# ./executable fullName_suburb_alpha.csv
# ./executable strSize10_alpha.csv
# ./executable strSize26_alpha.csv
# ./executable strSize36.csv

# CT Health Data
# ./executable ds1_50k
# ./executable ds2_100k
# ./executable ds3_200k
# ./executable ds4_400k
# ./executable ds5_600k
# ./executable ds6_800k
# ./executable ds7_1M

# NC Voter data 

# nc_test_5k.csv
# NC_voterData_5M_Source_Annotated.csv

# Biological
# eColi_genes

# RandomStrings

# randomString_size26_alpha5.txt
# randomString_size26_alpha10.txt
# randomString_size26_alpha15.txt
# randomString_size26_alpha20.txt
# randomString_size26_alpha26.txt
# randomString_size20_alpha5.txt
# randomString_size20_alpha10.txt
# randomString_size20_alpha15.txt
# randomString_size20_alpha20.txt
# randomString_size20_alpha26.txt
# randomString_size15_alpha5.txt
# randomString_size15_alpha10.txt
# randomString_size15_alpha15.txt
# randomString_size15_alpha20.txt
# randomString_size15_alpha26.txt
# randomString_size10_alpha5.txt
# randomString_size10_alpha10.txt
# randomString_size10_alpha15.txt
# randomString_size10_alpha20.txt
# randomString_size10_alpha26.txt
# randomString_size5_alpha5.txt
# randomString_size5_alpha10.txt
# randomString_size5_alpha15.txt
# randomString_size5_alpha20.txt
# randomString_size5_alpha26.txt