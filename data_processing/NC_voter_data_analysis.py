import math

import numpy as np
import pandas as pd
import csv
import random
from math import comb
import re


def get_data_clusters(path, total_cluster_vec):
    counter = 0
    with open(path, "r", encoding="utf8") as file:
        csv_reader = csv.reader(file, delimiter=",")
        is_first_row = True
        data_counter = 0
        for row in csv_reader:
            if is_first_row:
                print(row)
                is_first_row = False
            else:
                data_counter +=1
                cluster_vec = []
                is_first_element = True
                for field in row:
                    if is_first_element:
                        cluster_vec.append(int(field))
                        is_first_element = False
                    else:
                        field = re.sub('[^A-Za-z0-9]+', '', field)
                        cluster_vec.append(field)
                total_cluster_vec.append(cluster_vec)
                counter += 1
                # if counter >= 5000:
                #     break
    print(f"Data points read: {data_counter}")
    return total_cluster_vec


# Takes input vector of clusters
# Outputs (i) a dictionary where recid are keys and a tuple is the index of the records of that key
# Outputs (ii) a dictionary where recid are keys and a value is the count of records found with that recid

def get_ssn_info(cluster_vec):
    recid_exact_indices = {}
    recid_counts = {}
    for i in range(0, len(cluster_vec)):
        if cluster_vec[i][0] in recid_exact_indices:
            recid_counts[cluster_vec[i][0]] += 1
            recid_exact_indices[cluster_vec[i][0]].append(i)
        else:
            recid_exact_indices[cluster_vec[i][0]] = []
            # recid_exact_indices[cluster_vec[i][0]].append(i)
            recid_counts[cluster_vec[i][0]] = 1
    return recid_exact_indices, recid_counts


def get_cluster_size_stat(recid_counts):
    t_counts = {}
    for key,val in recid_counts.items():
        if val in t_counts:
            t_counts[val] += 1
        else:
            t_counts[val] = 1
    return t_counts


def get_all_records_with_ssn_count(count, tot_rec_dict, count_dict):
    clusters = []
    for key,val in count_dict.items():
        if val == count:
            for indices in tot_rec_dict[key]:
                clusters.append(cluster_members_vec[indices])
    return clusters


def check_mismatch_and_count(sorted_data, counter):
    first_name_miss = 0
    last_name_miss = 0
    address_miss = 0
    zipcode_miss = 0
    for i in range(0, len(sorted_data)-1):
        if sorted_data[i][0] == sorted_data[i+1][0]:
            if sorted_data[i][1] != sorted_data[i+1][1]:
                counter +=1
                first_name_miss += 1
                # print("First Name mismatch")
                # print(sorted_data[i])
                # print(sorted_data[i+1])
                # print()
            if sorted_data[i][2] != sorted_data[i+1][2]:
                counter += 1
                last_name_miss += 1
                # print("last Name mismatch")
                # print(sorted_data[i])
                # print(sorted_data[i + 1])
                # print()
            if sorted_data[i][3] != sorted_data[i + 1][3]:
                counter += 1
                address_miss += 1
                # print("Address mismatch")
                # print(sorted_data[i])
                # print(sorted_data[i + 1])
                # print()
            if sorted_data[i][4] != sorted_data[i + 1][4]:
                counter += 1
                zipcode_miss += 1
                # print("Post Code mismatch")
                # print(sorted_data[i])
                # print(sorted_data[i + 1])
                # print()
    print(f"First Name Miss: {first_name_miss}")
    print(f"Last Name Miss: {last_name_miss}")
    print(f"Address Miss: {address_miss}")
    print(f"Zip Code Miss: {zipcode_miss}")
    print(f"Total mismatch: {counter}")

# Read Data
file_path_1 = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/NC_voter_data_5M/vds_1_1M.csv"
file_path_2 = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/NC_voter_data_5M/vds_2_1M.csv"
file_path_3 = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/NC_voter_data_5M/vds_3_1M.csv"
file_path_4 = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/NC_voter_data_5M/vds_4_1M.csv"
file_path_5 = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/NC_voter_data_5M/vds_5_1M.csv"

cluster_members_vec = []
cluster_members_vec = get_data_clusters(file_path_1, cluster_members_vec)
print("first file read")
cluster_members_vec = get_data_clusters(file_path_2, cluster_members_vec)
print("second file read")
cluster_members_vec = get_data_clusters(file_path_3, cluster_members_vec)
print("third file read")
cluster_members_vec = get_data_clusters(file_path_4, cluster_members_vec)
print("fourth file read")
cluster_members_vec = get_data_clusters(file_path_5, cluster_members_vec)
print("fifth file read")
# cluster_1 = cluster_members_vec[: len(cluster_members_vec) // 2]
# cluster_2 = cluster_members_vec[len(cluster_members_vec) // 2:]
#
# with open("NC_voter_data_1_1M.csv", "w", newline="") as f:
#     writer = csv.writer(f)
#     writer.writerows(cluster_1)
#
# with open("NC_voter_data_2_1M.csv", "w", newline="") as f:
#     writer = csv.writer(f)
#     writer.writerows(cluster_2)

cluster_members_vec.sort(key=lambda x: int(x[0]))
print("List Sorted")
recid_tot_dict, recid_count_dict = get_ssn_info(cluster_members_vec)
print("Dictionaries built")
selected_clusters = get_all_records_with_ssn_count(5, recid_tot_dict, recid_count_dict)
print(selected_clusters)
random.shuffle(selected_clusters)
print(selected_clusters)

# recid_sizes_dict = get_cluster_size_stat(recid_count_dict)
# for key,val in recid_sizes_dict.items():
#     print(f"{val} ids have {key} records")
# mismatch_counter = 0
# check_mismatch_and_count(cluster_members_vec, mismatch_counter)
print(f"Total Records: {len(cluster_members_vec)}")
print(f"Total Selected Records: {len(selected_clusters)}")

with open("NC_voter_data_mixed_1M.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(selected_clusters)

cluster_1 = selected_clusters[: len(selected_clusters) // 2]
cluster_2 = selected_clusters[len(selected_clusters) // 2:]

with open("NC_voter_data_mixed_1_1M.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(cluster_1)

with open("NC_voter_data_mixed_2_1M.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(cluster_2)