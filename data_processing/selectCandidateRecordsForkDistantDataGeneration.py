import math
import numpy as np
import pandas as pd
import csv
import random
from math import comb


def get_my_clusters(path):
    total_cluster_vec = []
    with open(path, "r", encoding="utf8") as file:
        csv_reader = csv.reader(file, delimiter=",")
        for row in csv_reader:
            cluster_vec = []
            if len(row) > 0:
                for uid in row:
                    if len(uid) > 1:
                        cluster_vec.append(uid)
            cluster_vec.sort()
            total_cluster_vec.append(cluster_vec)
    return total_cluster_vec


def get_data(path):
    ssn_dict = {}
    with open(path, "r", encoding="utf8") as first_file:
        tsv_reader = csv.reader(first_file, delimiter=",")
        for row in tsv_reader:
            ssn = int(row[0])
            if ssn in ssn_dict:
                continue
            else:
                ssn_dict[ssn] = row
    return ssn_dict

path_clusters = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/Server_results/debug/out_single_linkage_ds7_1M_parallel_normal_blocking_fullBlocking_unionFind_6_threads"
path_data = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/ds7_1M"

my_clusters = get_my_clusters(path_clusters)
records_ssn_dict = get_data(path_data)

selected_records = []

for x in my_clusters:
    ssn = int(x[0])
    if ssn in records_ssn_dict:
        selected_records.append(records_ssn_dict[ssn])
    else:
        print(f"SSN {ssn} Not Found!")

df = pd.DataFrame(selected_records)
df.to_csv("unique_records.csv", index=False, header=False)
print(len(selected_records))
