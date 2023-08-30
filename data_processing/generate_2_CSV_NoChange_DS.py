import numpy as np
import pandas as pd
import csv
import random
import re


def get_data_and_write_csv(path1, path2, filename):
    ssn_list = []
    first_name_list = []
    last_name_list = []
    dod_list = []
    dob_list = []
    with open(path1, "r", encoding="utf8") as first_file:
        tsv_reader = csv.reader(first_file, delimiter="\t")
        for row in tsv_reader:
            ssn_list.append(row[0])
            regex = re.compile('[^a-zA-Z]')
            # First parameter is the replacement, second parameter is your input string
            my_str = regex.sub('', row[1])
            my_str = my_str.lower()
            first_name_list.append(my_str)
            my_str = regex.sub('', row[2])
            my_str = my_str.lower()
            last_name_list.append(my_str)
            dod_list.append(row[3])
            dob_list.append(row[4])
    df1 = pd.DataFrame(list(zip(ssn_list, first_name_list, last_name_list, dod_list, dob_list)),
                       columns=['SSN', 'FirstName', 'LastName', 'DoD', 'DoB'])
    filename1 = filename + "ds11_1_5M"
    df1.to_csv(filename1, index=False, header=False)
    ssn_list = []
    first_name_list = []
    last_name_list = []
    dod_list = []
    dob_list = []

    with open(path2, "r", encoding="utf8") as second_file:
        tsv_reader = csv.reader(second_file, delimiter="\t")
        for row in tsv_reader:
            ssn_list.append(row[0])
            regex = re.compile('[^a-zA-Z]')
            # First parameter is the replacement, second parameter is your input string
            my_str = regex.sub('', row[1])
            my_str = my_str.lower()
            first_name_list.append(my_str)
            my_str = regex.sub('', row[2])
            my_str = my_str.lower()
            last_name_list.append(my_str)
            dod_list.append(row[3])
            dob_list.append(row[4])

    df2 = pd.DataFrame(list(zip(ssn_list, first_name_list, last_name_list, dod_list, dob_list)),
                      columns=['SSN', 'FirstName', 'LastName', 'DoD', 'DoB'])
    filename2 = filename + "ds7_1M_original"
    df2.to_csv(filename2, index=False, header=False)


def get_data_and_write_one_csv(path1, path2, filename):
    ssn_list = []
    first_name_list = []
    last_name_list = []
    dod_list = []
    dob_list = []
    with open(path1, "r", encoding="utf8") as first_file:
        tsv_reader = csv.reader(first_file, delimiter="\t")
        for row in tsv_reader:
            ssn_list.append(row[0])
            # regex = re.compile('[^a-zA-Z]')
            # First parameter is the replacement, second parameter is your input string
            # my_str = regex.sub('', row[1])
            # my_str = my_str.lower()
            first_name_list.append(row[1])
            # my_str = regex.sub('', row[2])
            # my_str = my_str.lower()
            last_name_list.append(row[2])
            dod_list.append(row[3])
            dob_list.append(row[4])
    # df1 = pd.DataFrame(list(zip(ssn_list, first_name_list, last_name_list, dod_list, dob_list)),
    #                    columns=['SSN', 'FirstName', 'LastName', 'DoD', 'DoB'])
    # filename1 = filename + "ds11_1_5M"
    # df1.to_csv(filename1, index=False, header=False)
    # ssn_list = []
    # first_name_list = []
    # last_name_list = []
    # dod_list = []
    # dob_list = []

    with open(path2, "r", encoding="utf8") as second_file:
        tsv_reader = csv.reader(second_file, delimiter="\t")
        for row in tsv_reader:
            ssn_list.append(row[0])
            # regex = re.compile('[^a-zA-Z]')
            # First parameter is the replacement, second parameter is your input string
            # my_str = regex.sub('', row[1])
            # my_str = my_str.lower()
            first_name_list.append(row[1])
            # my_str = regex.sub('', row[2])
            # my_str = my_str.lower()
            last_name_list.append(row[2])
            dod_list.append(row[3])
            dob_list.append(row[4])

    df2 = pd.DataFrame(list(zip(ssn_list, first_name_list, last_name_list, dod_list, dob_list)),
                      columns=['SSN', 'FirstName', 'LastName', 'DoD', 'DoB'])
    filename2 = filename + "ds7_1M_original"
    df2.to_csv(filename2, index=False, header=False)


###### Main Func ######
###### This takes tab deliminated <SSN, FirstName, LastName, Dod, Dob> 2 files and merges to one with <SSN, Name, Dod, Dob>
file_path_1 = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/data_processing/ds7.1.1"
file_path_2 = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/data_processing/ds7.1.2"
file_path_3 = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/"
# get_data_and_write_csv(file_path_1, file_path_2, file_path_3)
get_data_and_write_one_csv(file_path_1, file_path_2, file_path_3)