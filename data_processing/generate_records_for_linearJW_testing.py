import math

import numpy as np
import pandas as pd
import csv
import random
from math import comb
import re
import string

random.seed(1)

def get_records(path):
    records = []
    with open(path, "r", encoding="utf8") as file:
        tsv_reader = csv.reader(file, delimiter=",")
        for row in tsv_reader:
            record = []
            regex = re.compile('[^a-zA-Z]')
            my_str = regex.sub('', row[1])
            my_str = my_str.lower()
            record.append(my_str)
            my_str = regex.sub('', row[2])
            my_str = my_str.lower()
            record.append(my_str)
            my_str = regex.sub('', row[3])
            my_str = my_str.lower()
            record.append(my_str)
            records.append(record)
    return records


def get_sampled_records(records, num):
    indices = random.sample(range(len(records)), num)
    sampled_records = [records[i] for i in indices]
    return sampled_records

def add_numerical_mix(records):
    mixed_records = []
    for rec in records:
        length_of_string = len(rec)
        example_list = []  # Create a temporary list
        example_list[:0] = rec  # Get string, turns it into list
        times = round(length_of_string * 0.30)
        for i in range (0, times):
            example_list.insert(random.randrange(0, len(rec)),
                            str(random.randrange(0, 9)))  # Insert a random number into a random position of the list

        list_into_string = ''.join(example_list)  # Get string, turns it into string
        mixed_records.append(list_into_string)
    return mixed_records


def trimStr(records, num):
    trimmedRecs = []
    for rec in records:
        trimmedRecs.append(rec[:num])
    return trimmedRecs


def main():
    path = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/NC_voterData_5M_Source_Annotated.csv"
    num_rec_sampled = 44722
    all_records = get_records(path)
    sampled_records = get_sampled_records(all_records, num_rec_sampled)
    # print(sampled_records)
    fName = []
    for rec in sampled_records:
        fName.append(rec[0])
    fName = add_numerical_mix(fName)
    df_fname = pd.DataFrame(fName)
    df_fname = df_fname.sample(frac=1).reset_index(drop=True)
    print(len(df_fname))
    # print(df_fname)

    fullName = []
    for rec in sampled_records:
        fullName.append(rec[0] + rec[1])
    fullName = add_numerical_mix(fullName)
    df_fullname = pd.DataFrame(fullName)
    df_fullname = df_fullname.sample(frac=1).reset_index(drop=True)
    print(len(df_fullname))
    # print(df_fullname)

    fullName_suburb = []
    for rec in sampled_records:
        fullName_suburb.append(rec[0] + rec[1] + rec[2])
    fullName_suburb = add_numerical_mix(fullName_suburb)
    df_fullname_suburb = pd.DataFrame(fullName_suburb)
    df_fullname_suburb = df_fullname_suburb.sample(frac=1).reset_index(drop=True)
    print(len(df_fullname_suburb))
    # print(df_fullname_suburb)

    string_10 = []
    string_26 = []
    string_36 = []

    for rec in fullName:
        while len(rec) < 36:
            rec = rec + rec
        string_36.append(rec[:36])
        string_26.append(rec[:26])
        string_10.append(rec[:10])
    string_36 = add_numerical_mix(string_36)
    string_26 = add_numerical_mix(string_26)
    string_10 = add_numerical_mix(string_10)

    string_10 = trimStr(string_10, 10)
    string_26 = trimStr(string_26, 26)
    string_36 = trimStr(string_36, 36)

    df_10 = pd.DataFrame(string_10)
    df_10 = df_10.sample(frac=1).reset_index(drop=True)
    print(len(df_10))
    # print(df_10)

    df_26 = pd.DataFrame(string_26)
    df_26 = df_26.sample(frac=1).reset_index(drop=True)
    print(len(df_26))
    # print(df_26)

    df_36 = pd.DataFrame(string_36)
    df_36 = df_36.sample(frac=1).reset_index(drop=True)
    print(len(df_36))
    # print(df_36)

    # filename = "f_name.csv"
    # df_fname.to_csv(filename, index=False, header=False)
    # filename = "fullName.csv"
    # df_fullname.to_csv(filename, index=False, header=False)
    # filename = "fullName_suburb.csv"
    # df_fullname_suburb.to_csv(filename, index=False, header=False)
    # filename = "strSize10.csv"
    # df_10.to_csv(filename, index=False, header=False)
    # filename = "strSize26.csv"
    # df_26.to_csv(filename, index=False, header=False)
    # filename = "strSize36.csv"
    # df_36.to_csv(filename, index=False, header=False)

    filename = "f_name_alpha.csv"
    df_fname.to_csv(filename, index=False, header=False)
    filename = "fullName_alpha.csv"
    df_fullname.to_csv(filename, index=False, header=False)
    filename = "fullName_suburb_alpha.csv"
    df_fullname_suburb.to_csv(filename, index=False, header=False)
    filename = "strSize10_alpha.csv"
    df_10.to_csv(filename, index=False, header=False)
    filename = "strSize26_alpha.csv"
    df_26.to_csv(filename, index=False, header=False)
    filename = "strSize36_alpha.csv"
    df_36.to_csv(filename, index=False, header=False)



if __name__=="__main__":
    main()