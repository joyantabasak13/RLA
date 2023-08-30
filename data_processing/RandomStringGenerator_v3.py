# Generates random string k sizes where k is larger than the alphabet size e.g. 25,50,75,100,125,150

import numpy as np
import pandas as pd
import csv
import random
import string
import re

random.seed(10)


def getRandomStrings(numOfStrings, size, characters):
    strings = []
    for i in range(0, numOfStrings):
        strings.append(''.join(random.choice(string.ascii_lowercase[:characters]) for _ in range(size)))
    print(strings)
    return strings


### main ###

numberOfSequences = 10000
stringSizeArr = [25, 50, 75, 100, 125, 150]
alphabetSizeArr = [5, 10, 15, 20, 26]
for set in range(1,6):
    for stringSize in stringSizeArr:
        for alphabetSize in alphabetSizeArr:
            sequences = getRandomStrings(numberOfSequences, stringSize, alphabetSize)
            print(f"Total sequences: {len(sequences)} of strSize {stringSize} and alphabetSize {alphabetSize}")
            filename = "randomString_" + str(stringSize) + "_alpha_" + str(alphabetSize) + "_set_" + str(set) + ".txt"

            with open(filename, 'w') as f:
                for seq in sequences:
                    f.write("%s\n" %seq)

# print(f"Number of Genes: {len(genes)}")