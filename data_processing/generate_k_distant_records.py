import math

import numpy as np
import pandas as pd
import csv
import random
from math import comb
import re
import string


def tag_original_records(records, tag):
    simulated_records = []
    for record in records:
        s_record = []
        for attr in record:
            s_record.append(attr)
        s_record.append(tag)
        simulated_records.append(s_record.copy())
    return simulated_records


def add_edit_dist_to_str(name, num_dist):
    for i in range(0, num_dist):
        index = random.randint(0, len(name)-1)
        operation = random.randint(0, 2)
        if operation == 0:
            # DO Delete
            if index == len(name) - 1:
                name = name[:index]
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
            print(name[index])
            name = name[:index-1] + character + name[index+1:]
    return name


def generate_records_with_change_in_ch(record, ch_dist_unit, ch_copy_num, ch_indices):
    records = []
    for i in range(0, ch_copy_num):
        for j in range(0,len(ch_dist_unit))
        for x in ch_indices:
            sim_rec = records[-1]
            sim_rec[x] = add_edit_dist_to_str(sim_rec[x], ch_dist_unit)
            records.append(sim_rec)
    return records


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


def generate_records_with_change_in_date(record, date_dist_unit, date_copy_num, date_indices):
    records = []
    for i in range(0, date_copy_num):
        sim_rec = record.copy()
        for x in date_indices:
            sim_rec[x] = add_edit_dist_to_date(sim_rec[x], date_dist_unit)
        records.append(sim_rec)
    return records


def generate_simulated_records(records, exact_copy_num, ch_dist_unit, ch_copy_num, ch_indices, date_dist_unit, date_copy_num, date_indices):
    sim_records = []
    for i in range(0, len(records)):
        record = records[i]
        # make n exact copies
        for j in range(0, exact_copy_num):
            sim_records.append(record)
        # print("Printing Simulated Records in gen func")
        # for x in sim_records:
        #     print(x)
        # make ch_copy_num * len(ch_dist_unit) mutated copies in each of the ch_indices
        ch_records = generate_records_with_change_in_ch(record.copy(), ch_dist_unit, ch_copy_num, ch_indices)
        # make date_copy_num copies each date_dist_unit digit mutated in date_indices
        date_records = generate_records_with_change_in_date(record.copy(), date_dist_unit, date_copy_num, date_indices)
        for x in ch_records:
            sim_records.append(x)
        for y in date_records:
            sim_records.append(y)
    return sim_records


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
    original_records = [['0001', 'joyantabasak', '03181995', '03181895']]
    # add tag to original records
    original_tagged_records = tag_original_records(original_records.copy(), original_tag)

    # get simulated records
    simulated_records = []
    for record in original_tagged_records:
        simulated_records.append(record.copy())

    # generate exact copies
    num_exact_copy = 1
    exact_record_copies = generate_exact_copies(original_tagged_records, num_exact_copy, exact_copy_tag)

    for copy_record in exact_record_copies:
        simulated_records.append(copy_record.copy())

    # Resume From Here

    
    # generate records with mutated non-blocking atrributes
    non_blocking_dist = 1
    non_blocking_copies = 1
    non_blocking_attribute_indices = [2,3]


    # blockfieldCopy number of record for each of the blockField dist in each of the blockFieldIndices
    # block field dist is
    blockfieldDist = [1, 2]
    blockfieldCopy = 1
    dateFieldDist = 1
    dateFieldCopy = 0
    blockfieldAttrIndices = [1]
    dateFieldDAttrIndices = [2,3]
    simulatedRecords = generate_simulated_records(originalRecords, numExactCopy, blockfieldDist, blockfieldCopy, blockfieldAttrIndices, dateFieldDist, dateFieldCopy, dateFieldDAttrIndices)
    print("Printing Simulated Records")
    for x in simulatedRecords:
        print(x)


if __name__=="__main__":
    main()