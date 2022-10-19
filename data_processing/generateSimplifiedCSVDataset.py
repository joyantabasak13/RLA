import numpy as np
import pandas as pd
import csv
import random


def get_data_and_write_csv(path1, path2, filename):
    ssn_list = []
    name_list = []
    with open(path1, "r", encoding="utf8") as first_file:
        tsv_reader = csv.reader(first_file, delimiter="\t")
        for row in tsv_reader:
            ssn_list.append(row[0])
            my_str = row[1] + row[2]
            my_str = my_str.replace(" ", "")
            my_str = my_str.lower()
            name_list.append(my_str)
    df1 = pd.DataFrame(list(zip(ssn_list, name_list)),
                       columns=['SSN', 'Name'])
    filename1 = filename + "_" + str(1)
    df1.to_csv(filename1, index=False, header=False)
    ssn_list = []
    name_list = []
    with open(path2, "r", encoding="utf8") as second_file:
        tsv_reader = csv.reader(second_file, delimiter="\t")
        for row in tsv_reader:
            ssn_list.append(row[0])
            my_str = row[1] + row[2]
            my_str = my_str.replace(" ", "")
            my_str = my_str.lower()
            name_list.append(my_str)

    df2 = pd.DataFrame(list(zip(ssn_list, name_list)),
                       columns=['SSN', 'Name'])
    filename2 = filename + "_" + str(2)
    df2.to_csv(filename2, index=False, header=False)


###### Main Func ######

file_path_1 = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/19500671/ds11.1.1"
file_path_2 = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/19500671/ds11.1.2"
file_path_3 = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/data_1M"

get_data_and_write_csv(file_path_1, file_path_2, file_path_3)
