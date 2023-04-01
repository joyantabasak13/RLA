import math

import numpy as np
import pandas as pd
import csv
import random
from math import comb
import re
import string


def get_unique_records(path):
    records = []
    with open(path, "r", encoding="utf8") as first_file:
        tsv_reader = csv.reader(first_file, delimiter=",")
        for row in tsv_reader:
            records.append(row)
    return records


def get_original_records(records, num):
    indices = random.sample(range(len(records)), num)
    sampled_records = [records[i] for i in indices]
    return sampled_records


def tag_original_records(records, tag):
    simulated_records = []
    for record in records:
        s_record = []
        for attr in record:
            s_record.append(attr)
        s_record.append(tag)
        simulated_records.append(s_record.copy())
    return simulated_records


def generate_exact_copies(original_records, num_copy, tag):
    simulated_records = []
    for record in original_records:
        for i in range(0, num_copy):
            s_record = []
            for attr in record:
                s_record.append(attr)
            s_record[-1] = tag
            simulated_records.append(s_record.copy())
    return simulated_records


def add_edit_dist_to_date(date, num_dist):
    for i in range(0, num_dist):
        index = random.randint(0, len(date)-1)
        # Do Mutate
        digit = random.choice(string.digits)
        while str(digit) == date[index]:
            digit = random.choice(string.digits)
        # print(f"Date: {date} Length {len(date)}")
        # print(f"Index {index} Replacing {date[index]} by {str(digit)}")
        if index == 0:
            date = str(digit) + date[1:]
            # print(f"Adding to first {date}")
        elif index == len(date) - 1:
            date = date[:index] + str(digit)
            # print(f"Adding to last {date}")
        else:
            date = date[:index] + str(digit) + date[index+1:]
    #         print(f"Adding to middle {date}")
    # print(f"Changed Date: {date}")
    return date


def generate_record_with_change_in_date(record, dist, date_indices, tag):
    sim_rec = record.copy()
    for x in date_indices:
        sim_rec[x] = add_edit_dist_to_date(sim_rec[x], dist)
    sim_rec[-1] = tag
    return sim_rec


def generate_non_blocking_attribute_edited_records(original_records, dist, attr_indices, tag):
    simulated_records = []
    for record in original_records:
        s_record = generate_record_with_change_in_date(record, dist, attr_indices, tag)
        simulated_records.append(s_record)
    return simulated_records


def add_edit_dist_to_str(name, num_dist):
    for i in range(0, num_dist):
        index = random.randint(0, len(name)-1)
        operation = random.randint(0, 2)
        if operation == 0:
            # DO Delete
            if index == len(name) - 1:
                name = name[:-1]
            else:
                name = name[:index] + name[index+1:]
        if operation == 1:
            #DO Addition
            character = random.choice(string.ascii_lowercase)
            if index == len(name) - 1:
                name = name + character
            else:
                name = name[:index] + character + name[index:]
        if operation == 2:
            # Do mutate
            character = random.choice(string.ascii_lowercase)
            while character == name[index]:
                character = random.choice(string.ascii_lowercase)
            if index == 0:
                name = character + name[1:]
            elif index == len(name) - 1:
                name = name[:index] + character
            else:
                name = name[:index] + character + name[index + 1:]
    return name


def generate_record_with_change_in_ch(record, dist, attr_indices, tag):
    sim_rec = record.copy()
    for x in attr_indices:
        sim_rec[x] = add_edit_dist_to_str(sim_rec[x], dist)
    sim_rec[-1] = tag
    return sim_rec


def generate_block_attr_edited_records(original_records, dist, attr_indices, tag):
    simulated_records = []
    for record in original_records:
        s_record = generate_record_with_change_in_ch(record, dist, attr_indices, tag)
        simulated_records.append(s_record)
    return simulated_records


def main():
    # tags
    original_tag = '0'
    exact_copy_tag = '1'
    non_blocking_edit_tag = '2'
    blocking_edit_1_tag = '01'
    incremental_blocking_edit_1_tag = '011'
    blocking_edit_2_tag = '02'
    blocking_edit_1_with_non_blocking_edit_1_tag = '21'
    incremental_blocking_edit_1_with_non_blocking_edit_1_tag = '211'
    blocking_edit_2_with_non_blocking_edit_1_tag = '22'

    # get original records
    target_number_of_original_records = 80000
    original_records = [['0001', 'joyantabasak', '03181996', '03181895']]
    unique_records = get_unique_records("unique_records.csv")
    original_records = get_original_records(unique_records, target_number_of_original_records)
    # print(original_records)

    # add tag to original records
    original_tagged_records = tag_original_records(original_records.copy(), original_tag)

    # get simulated records
    simulated_records = []
    for record in original_tagged_records:
        simulated_records.append(record.copy())

    # generate exact copies
    num_exact_copy = 0
    exact_record_copies = generate_exact_copies(original_tagged_records, num_exact_copy, exact_copy_tag)

    for copy_record in exact_record_copies:
        simulated_records.append(copy_record.copy())

    # generate records with 1 edit distance in blocking attr
    block_attr_dist = 1
    block_attr_indices = [1]
    block_attr_edit_distance_1_records_2 = generate_block_attr_edited_records(original_tagged_records.copy(),
                                                                            block_attr_dist, block_attr_indices,
                                                                            blocking_edit_1_tag)

    block_attr_edit_distance_1_records = generate_block_attr_edited_records(original_tagged_records.copy(),
                                                                            block_attr_dist, block_attr_indices,
                                                                            blocking_edit_1_tag)

    for edited_records in block_attr_edit_distance_1_records:
        simulated_records.append(edited_records.copy())
    for edited_records in block_attr_edit_distance_1_records_2:
        simulated_records.append(edited_records.copy())

    # generate records with edited non-blocking atributes on top of 1 edit distance in blockfield
    non_blocking_dist = 1
    non_blocking_attribute_indices = [2,3]
    non_blocking_attribute_edited_records = generate_non_blocking_attribute_edited_records(block_attr_edit_distance_1_records.copy(), non_blocking_dist, non_blocking_attribute_indices, non_blocking_edit_tag)

    non_blocking_attribute_edited_records_2 = generate_non_blocking_attribute_edited_records(
        block_attr_edit_distance_1_records_2.copy(), non_blocking_dist, non_blocking_attribute_indices,
        non_blocking_edit_tag)

    for edited_records in non_blocking_attribute_edited_records:
        simulated_records.append(edited_records.copy())
    for edited_records in non_blocking_attribute_edited_records_2:
        simulated_records.append(edited_records.copy())


    # print("Printing Simulated Records")
    # for x in simulated_records:
    #     print(x)

    # filename = "simulated_records_" + str(len(df)) + ".csv"
    df = pd.DataFrame(simulated_records)
    df = df.sample(frac=1).reset_index(drop=True)
    print(len(df))
    filename = "simulated_records_UnionFindTest_" + str(target_number_of_original_records) + "_records.csv"
    df.to_csv(filename, index=False, header=False)


if __name__=="__main__":
    main()