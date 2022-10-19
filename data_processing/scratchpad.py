import numpy as np
import pandas as pd
import csv
import random

arr = [999000, 24995000, 99990000, 624975000, 2499950000]
perc = [1.0, 0.95, 0.90, 0.85, 0.80, 0.75]

for x in arr:
    for y in perc:
        print("FOR: "+ str(x) + " and PERC: "+str(y)+ " The Result is: "+ str(x*y))