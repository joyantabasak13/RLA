import math

import numpy as np
import pandas as pd
import csv
import random
from math import comb
import re


def getTxtData(filePath):
    file = open(filePath, 'r')
    count = 0
    linesC = 0
    smallest = 99999999
    longest = 0

    while True:
        # Get next line from file
        line = file.readline()
        line = line.strip()
        if not line:
            break
        count += len(line)
        linesC += 1

        if len(line) < smallest:
            smallest = len(line)
        if len(line) > longest:
            longest = len(line)

    file.close()
    print(f"Total char: {count} Total lines {linesC}")
    print(f"Smallest string length {smallest}")
    print(f"Longest String length {longest}")
    print(f"Average String length {count/linesC}")


def getCSVData(filePath, numAttributes):
    smallest = [9999999]*numAttributes
    largest = [0]*numAttributes
    total = [0]*numAttributes
    countRow = 0
    with open(filePath, "r", encoding="utf8") as first_file:
        tsv_reader = csv.reader(first_file, delimiter=",")
        for row in tsv_reader:
            countRow += 1
            for attr in range(0,numAttributes):
                curStrLen = len(row[attr])
                total[attr] = total[attr] + curStrLen
                if curStrLen < smallest[attr] and curStrLen!=0:
                    smallest[attr] = curStrLen
                if curStrLen > largest[attr]:
                    largest[attr] = curStrLen
    print(f"Total Rows: {countRow}")
    for attr in range(0,numAttributes):
        print(f"Attribute {attr} ")
        print(f"Smallest: {smallest[attr]}")
        print(f"Largest: {largest[attr]}")
        print(f"Average: {total[attr] / countRow}")


# Read Data
# file_path_1 = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/data_processing/eColi_genes.txt"

# getTxtData(file_path_1)

filePath2 = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/NC_voterData_5M_Source_Annotated.csv"

getCSVData(filePath2, 6)
