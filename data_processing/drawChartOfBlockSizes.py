import math

import numpy as np
import pandas as pd
import csv
import random
from math import comb
import re

import numpy as np
import matplotlib.pyplot as plt


counter = 0
path = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/assignment_2/NumCompPerBlock_1M_CTHealth"
blockID = []
blockSizes = []
with open(path, "r", encoding="utf8") as file:
    csv_reader = csv.reader(file, delimiter=",")
    for row in csv_reader:
        blockID.append(int(row[0]))
        blockSizes.append(int(row[1]))

blockID = np.array(blockID)
blockSizes = np.array(blockSizes)


print(f"Data points read: {len(blockID)}")
count = 0
for i in range(0, len(blockSizes)):
    count += blockSizes[i]
print(f"Count {count}")

### Plot all blocks
# fig = plt.figure()
# plt.bar(blockID, blockSizes, color='red')
#
# plt.xlabel("BlockID")
# plt.ylabel("#Comparisions in the Block")
# plt.show()

### Plot large Blocks
norm_blockSizes = []
blockSizes = np.sort(blockSizes)
for i in range(0, len(blockSizes)):
    norm_blockSizes.append(float((blockSizes[i]*100.0) / (count*1.0)))
norm_blockSizes = np.sort(norm_blockSizes)
print(f"Largest % {norm_blockSizes[-1]}")

fig = plt.figure()
plt.bar(blockID[:100], norm_blockSizes[-100:], color='red')

plt.xlabel("BlockID")
plt.ylabel("% Total Comparisions")
plt.show()


