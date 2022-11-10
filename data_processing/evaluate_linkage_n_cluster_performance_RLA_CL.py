import math

import numpy as np
import pandas as pd
import csv
import random
from math import comb


def get_data(path, max_size):
    total_component = 0
    correct_component_size = [0] * max_size
    incorrect_component_size = [0] * max_size
    is_new_cluster = False
    cur_ssn = ''
    cur_comp_size = 0
    is_wrong = False
    with open(path, "r", encoding="utf8") as file:
        csv_reader = csv.reader(file, delimiter="\t")
        for row in csv_reader:
            if len(row) == 1:
                total_component = total_component + 1
                is_new_cluster = True
                is_wrong = False
                cur_comp_size = 0
                cur_ssn = ''
            if len(row) > 1:
                if is_new_cluster:
                    cur_ssn = row[0]
                    cur_comp_size = 1
                    is_new_cluster = False
                else:
                    cur_comp_size = cur_comp_size + 1
                    if row[0] != cur_ssn:
                        is_wrong = True
            if len(row) < 1:
                if is_wrong == False:
                    if cur_comp_size > 0:
                        correct_component_size[cur_comp_size - 1] = correct_component_size[cur_comp_size - 1] + 1
                else:
                    if cur_comp_size > 0:
                        continue
                        #incorrect_component_size[cur_comp_size - 1] = incorrect_component_size[cur_comp_size - 1] + 1
                cur_comp_size = 0
    total_component = total_component - 1
    print(total_component)
    print(correct_component_size)
    print(incorrect_component_size)

    return total_component, correct_component_size, incorrect_component_size


def get_cluster_performance(total_component, correct_component_size, incorrect_component_size, expected_size, max_size):
    print("Total Components: " + str(total_component))
    for i in range(expected_size):
        frac = correct_component_size[i] / total_component
        frac = frac * 100.0
        print("Correct Comp of Size" + str(i + 1) + " : " + str(correct_component_size[i]) + ' (' + "{:.2f}".format(
            frac) + '%)')
        # print("Frac of   Correct Comp of Size" + str(i+1) + " : " + str(frac))
    # for i in range(max_size):
    #     frac = incorrect_component_size[i] / total_component
    #     frac = frac * 100.0
    #     print("Incorrect Comp of Size" + str(i + 1) + " : " + str(incorrect_component_size[i]) + ' (' + "{:.2f}".format(
    #         frac) + '%)')
        # print("Frac of   InCorrect Comp of Size" + str(i+1) + " : " + str(frac))
    return


def get_clusters(path):
    total_cluster_vec = []
    cluster_vec = []
    with open(path, "r", encoding="utf8") as file:
        csv_reader = csv.reader(file, delimiter="\t")
        for row in csv_reader:
            if len(row) == 1:
                cluster_vec = []
            if len(row) > 1:
                cur_ssn = row[0]
                cluster_vec.append(cur_ssn)
            if len(row) < 1:
                if len(cluster_vec) != 0:
                    cluster_vec.sort()
                    total_cluster_vec.append(cluster_vec)
                    cluster_vec = []
    # ruru
    # for i in total_cluster_vec:
    #     print(i)
    return total_cluster_vec


def get_clusterwise_ssn_count(total_clusters):
    cluster_ssn_count_vec = []
    for i in range(len(total_clusters)):
        for j in range(len(total_clusters[i])):
            curSSN = total_clusters[i][j]
            if j == 0:
                on_going_ssn = total_clusters[i][j]
                curSSN = on_going_ssn
                curr_ssn_count = 0

            if curSSN == on_going_ssn:
                curr_ssn_count = curr_ssn_count + 1
            if (curSSN != on_going_ssn) or (j == len(total_clusters[i]) - 1):
                cluster_ssn_count_vec.append(curr_ssn_count)
                # print(f"Cluster: {i} Has {on_going_ssn} times {curr_ssn_count}")
                curr_ssn_count = 1
                on_going_ssn = curSSN
    return cluster_ssn_count_vec


def get_linkage_tp(clusters_ssn_count):
    tp_count = 0
    for i in range(len(clusters_ssn_count)):
        tp_count = tp_count + math.comb(clusters_ssn_count[i], 2)
    return tp_count


def get_groups_in_each_cluster(clusters):
    cluster_ssn_groups = []
    for i in range(len(clusters)):
        cur_group = []
        for j in range(len(clusters[i])):
            curSSN = clusters[i][j]
            if j == 0:
                on_going_ssn = clusters[i][j]
                curSSN = on_going_ssn
                curr_ssn_count = 0

            if curSSN == on_going_ssn:
                curr_ssn_count = curr_ssn_count + 1
            if (curSSN != on_going_ssn) or (j == len(clusters[i]) - 1):
                cur_group.append(curr_ssn_count)
                # print(f"Cluster: {i} Has {on_going_ssn} times {curr_ssn_count}")
                curr_ssn_count = 1
                on_going_ssn = curSSN
        cluster_ssn_groups.append(cur_group)
    return cluster_ssn_groups


def get_linakge_fp(cluster_ssn_group_list):
    tot_fp = 0
    for i in range(len(cluster_ssn_group_list)):
        cluster_size = 0
        cluster_links = 0
        for j in cluster_ssn_group_list[i]:
            cluster_size = cluster_size + j
            cluster_links = cluster_links + math.comb(j,2)
        cur_fp = math.comb(cluster_size, 2) - cluster_links
        tot_fp = tot_fp + cur_fp
    return tot_fp


def get_ssn_groups(clusters):
    cluster_ssn_groups = {}
    for i in range(len(clusters)):
        for j in range(len(clusters[i])):
            curSSN = clusters[i][j]
            if not curSSN in cluster_ssn_groups:
                cluster_ssn_groups[curSSN] = []
            if j == 0:
                on_going_ssn = clusters[i][j]
                curSSN = on_going_ssn
                curr_ssn_count = 0

            if curSSN == on_going_ssn:
                curr_ssn_count = curr_ssn_count + 1
            if (curSSN != on_going_ssn) or (j == len(clusters[i]) - 1):
                cluster_ssn_groups[on_going_ssn].append(curr_ssn_count)
                # print(f"Cluster: {i} Has {on_going_ssn} times {curr_ssn_count}")
                curr_ssn_count = 1
                on_going_ssn = curSSN
    return cluster_ssn_groups


def get_linkage_fn(size, group_dict):
    tot_fn = 0
    tot_links = math.comb(size, 2)
    for x in group_dict:
        cur_links = 0
        for i in group_dict[x]:
            cur_links = cur_links + math.comb(i,2)
        cur_fn = tot_links - cur_links
        tot_fn = tot_fn + cur_fn
        # print(f"for {x} with {group_dict[x]} cur_links = {cur_links} and Curfn={cur_fn}")
    return tot_fn
### main ###

#file_path = "/Users/joyanta/Documents/Research/Record_Linkage/codes/rla_cl_exact/obj/ClOutE20k_base36_overflow.txtOutFinal.txt"
file_path = "/Users/joyanta/Downloads/output_edit_ds2.1.txtOutSingle"

### Calculate Cluster Accuracy
max_component_size = 100
tot_comp, correct_comp, incorrect_comp = get_data(file_path, max_component_size)
number_of_records = 100000
expected_size = 5
get_cluster_performance(tot_comp, correct_comp, incorrect_comp, expected_size, max_component_size)
cluster_members_vec = get_clusters(file_path)
clusters_ssn_count = get_clusterwise_ssn_count(cluster_members_vec)
linkage_TP = get_linkage_tp(clusters_ssn_count)
ssn_cluster_groups = get_groups_in_each_cluster(cluster_members_vec)
linkage_FP = get_linakge_fp(ssn_cluster_groups)
ssn_group_dictionary = get_ssn_groups(cluster_members_vec)
linkage_FN = get_linkage_fn(expected_size, ssn_group_dictionary)
linkage_TN = math.comb(number_of_records, 2) - linkage_TP - linkage_FP - linkage_FN
linkage_accuracy = (linkage_TP + linkage_TN) / (linkage_TP+ linkage_TN + linkage_FP + linkage_FN)
linkage_f1_score = (linkage_TP*2) / (linkage_TP*2 + linkage_FP + linkage_FN)
print(f"Linkage accuracy {linkage_accuracy}, Linkage f1 score {linkage_f1_score}")
print(f"Linkage TP = {linkage_TP}, Linkage_FP = {linkage_FP}, Linkage_FN = {linkage_FN}, Linkage TN = {linkage_TN}")
