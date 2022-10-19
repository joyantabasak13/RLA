import numpy as np
import pandas as pd
import csv
import random


def get_max_component_size(path):
    max_size = 0
    with open(path, "r", encoding="utf8") as file:
        csv_reader = csv.reader(file, delimiter=",")
        for row in csv_reader:
            if max_size < len(row):
                max_size = len(row)
    return max_size


def get_data(path, max_size):
    total_component = 0
    correct_component_size = [0, 0, 0, 0, 0]
    incorrect_component_size = [0]*max_size
    with open(path, "r", encoding="utf8") as file:
        csv_reader = csv.reader(file, delimiter=",")
        for row in csv_reader:
            total_component = total_component + 1
            cur_comp_size = len(row)
            cur_ssn = row[0]
            isWrong = False
            for x in row:
                if cur_ssn != x:
                    isWrong = True
            if isWrong:
                incorrect_component_size[cur_comp_size - 1] = incorrect_component_size[cur_comp_size - 1] + 1
            else:
                correct_component_size[cur_comp_size - 1] = correct_component_size[cur_comp_size - 1] + 1
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

filename_prefix = "out_ds10000_"
f_exact_name = "out_ds200_exact"
filename_suffix = ['1.000000', '0.950000', '0.900000', '0.850000', '0.800000', '0.750000']

filename = f_exact_name
print("Performance of "+filename)
file_path = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/results/"+filename
max_component_size = get_max_component_size(file_path)
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
