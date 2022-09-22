import numpy as np
import pandas as pd
import csv


def get_data(path1, path2):
    with open(path1, "r", encoding="utf8") as first_file:
        i = 0
        tsv_reader = csv.reader(first_file, delimiter="\t")
        for row in tsv_reader:
            print(row[1])
            print(row)
            if i > 10:
                break
            i = i+1


###### Main Func ######

file_path_1 = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/19500671/ds11.1.1"
file_path_2 = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/19500671/ds11.1.2"

get_data(file_path_1, file_path_2)