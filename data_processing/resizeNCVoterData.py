import numpy as np
import pandas as pd
import csv
import random


def get_data(path):
    ssn_list = []
    fname_list = []
    lname_list = []
    suburb_list = []
    postcode_list = []
    src_list = []
    with open(path, "r", encoding="utf8") as file:
        tsv_reader = csv.reader(file, delimiter=",")
        for row in tsv_reader:
            ssn_list.append(row[0])
            fname_list.append(row[1])
            lname_list.append(row[2])
            suburb_list.append(row[3])
            postcode_list.append(row[4])
            src_list.append(row[5])

    df = pd.DataFrame(list(zip(ssn_list, fname_list, lname_list, suburb_list, postcode_list, src_list)),
                      columns=['SSN', 'FirstName', 'LastName', 'Suburb', 'PostCode', 'Source'])
    return df


###### Main Func ######

file_path = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/NC_voterData_5M_Source_Annotated.csv"

df = get_data(file_path)

unique_ssn_list = df["SSN"]
print(f"Total Records: {len(unique_ssn_list)}")
unique_ssn_list = list(dict.fromkeys(unique_ssn_list))
print(f"Total Clusters: {len(unique_ssn_list)}")

datasizes = [len(unique_ssn_list)//100, len(unique_ssn_list)//20, len(unique_ssn_list)//10]

for x in datasizes:
    ssn_list_x = random.sample(unique_ssn_list, x)
    df_x = df.loc[df['SSN'].isin(ssn_list_x)]
    file_name = "NC" + str(x)
    print(f"Selected {x} clusters that contains {len(df_x)} records")
    df_x.to_csv(file_name, index=False, header=False)
