import numpy as np
import pandas as pd
import csv
import random


# def get_max_component_size(path):
#     max_size = 0
#     with open(path, "r", encoding="utf8") as file:
#         csv_reader = csv.reader(file, delimiter=",")
#         for row in csv_reader:
#             if max_size < len(row):
#                 max_size = len(row)
#     return max_size


def get_data(path, max_size):
    total_component = 0
    correct_component_size = [0, 0, 0, 0, 0]
    incorrect_component_size = [0]*max_size
    is_new_cluster = False
    cur_ssn = ''
    cur_comp_size = 0
    is_wrong = False
    with open(path, "r", encoding="utf8") as file:
        csv_reader = csv.reader(file, delimiter="\t")
        for row in csv_reader:
            if len(row)==1:
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
                    if cur_comp_size>0:
                        correct_component_size[cur_comp_size - 1] = correct_component_size[cur_comp_size - 1] + 1
                else:
                    if cur_comp_size>0:
                        incorrect_component_size[cur_comp_size - 1] = incorrect_component_size[cur_comp_size - 1] + 1
                cur_comp_size = 0
    total_component = total_component - 1
    print(total_component)
    print(correct_component_size)
    print(incorrect_component_size)

    return total_component, correct_component_size, incorrect_component_size


def get_performance(total_component, correct_component_size, incorrect_component_size, expected_size, max_size):
    print("Total Components: "+str(total_component))
    for i in range(expected_size):
        frac = correct_component_size[i] / total_component
        frac = frac * 100.0
        print("Correct Comp of Size" + str(i+1) + " : " + str(correct_component_size[i]) + ' ('+"{:.2f}".format(frac)+'%)')
        # print("Frac of   Correct Comp of Size" + str(i+1) + " : " + str(frac))
    for i in range(max_size):
        frac = incorrect_component_size[i] / total_component
        frac = frac * 100.0
        print("Incorrect Comp of Size" + str(i+1) + " : " + str(incorrect_component_size[i])+ ' ('+"{:.2f}".format(frac)+'%)')
        # print("Frac of   InCorrect Comp of Size" + str(i+1) + " : " + str(frac))
    return


### main ###

# filename_prefix = "out_ds10000_"
# f_exact_name = "out_ds10000_exact"
# filename_suffix = ['1.000000', '0.950000', '0.900000', '0.850000', '0.800000', '0.750000']
#
# filename = f_exact_name
# print("Performance of "+filename)
file_path = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/results/rla_cl_clusters/10k/ClOutE10k_50Block_50rand.txt"
max_component_size = 50
tot_comp, correct_comp, incorrect_comp = get_data(file_path, max_component_size)
expected_size = 5
get_performance(tot_comp, correct_comp, incorrect_comp, expected_size, max_component_size)

# for x in filename_suffix:
#     filename = filename_prefix+x
#     print("Performance of "+filename)
#     file_path = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/results/"+filename
#     max_component_size = get_max_component_size(file_path)
#     tot_comp, correct_comp, incorrect_comp = get_data(file_path, max_component_size)
#     expected_size = 5
#     get_performance(tot_comp, correct_comp, incorrect_comp, expected_size, max_component_size)