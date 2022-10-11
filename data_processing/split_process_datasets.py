import numpy as np
import pandas as pd
import csv
import random


def get_data(path1, path2):
    ssn_list = []
    name_list = []
    with open(path1, "r", encoding="utf8") as first_file:
        tsv_reader = csv.reader(first_file, delimiter="\t")
        for row in tsv_reader:
            ssn_list.append(row[0])
            str = row[1] + row[2]
            str = str.replace(" ", "")
            str = str.lower()
            name_list.append(str)

    with open(path2, "r", encoding="utf8") as second_file:
        tsv_reader = csv.reader(second_file, delimiter="\t")
        for row in tsv_reader:
            ssn_list.append(row[0])
            str = row[1] + row[2]
            str = str.replace(" ", "")
            str = str.lower()
            name_list.append(str)

    df = pd.DataFrame(list(zip(ssn_list, name_list)),
                      columns=['SSN', 'Name'])
    return df


###### Main Func ######

file_path_1 = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/19500671/ds11.1.1"
file_path_2 = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/19500671/ds11.1.2"

df = get_data(file_path_1, file_path_2)

unique_ssn_list = df["SSN"]
unique_ssn_list = list(dict.fromkeys(unique_ssn_list))
datasizes = [7500, 10000]

for x in datasizes:
    ssn_list_x = random.sample(unique_ssn_list, x)
    df_x = df.loc[df['SSN'].isin(ssn_list_x)]
    file_name = "ds" + str(x)
    print(len(df_x))
    df_x.to_csv(file_name, index=False, header=False)
