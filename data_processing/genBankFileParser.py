import numpy as np
import pandas as pd
import csv
import random
import re


def getGeneLocations(path):
    all_limits = []
    gene_lines = []
    count = 0
    with open(path, "r", encoding="utf8") as file:
        Lines = file.readlines()
        for line in Lines:
            words = re.split('(\W+)', line)
            if words[2] == 'gene' and len(words[1].strip(' ')) == 0:
                gene_lines.append(words)
    # print(f"Total gene Lines: {len(gene_lines)}")
    # should be 4639
    for line in gene_lines:
        # print(line)
        limits = []
        for word in line:
            if word.isnumeric():
                limits.append(int(word))
        # print(limits)
        all_limits.append(limits)
    return all_limits


def getGenome(path):
    genome = ""
    count = 0
    with open(path, "r", encoding="utf8") as file:
        Lines = file.readlines()
        for line in Lines:
            if count > 0 :
                genome = genome + line.strip()
            count = count + 1
    return genome


### main ###

gene_file_path = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/data_processing/eColi.gbff"
genome_file_path = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/RLA/data_processing/eColi.fna"
gene_bounds = getGeneLocations(gene_file_path)
sequence = getGenome(genome_file_path)
print(f"Total genes: {len(gene_bounds)}")
print(f"Sequence size: {len(sequence)}")
genes = []
# for bound in gene_bounds:
#     genes.append(sequence[bound[0]:bound[1]])

sample_size = 900
ptr = 0
while ptr + sample_size < len(sequence):
    genes.append(sequence[ptr:ptr+sample_size])
    ptr = ptr+sample_size

with open('eColi_genes_900.txt', 'w') as f:
    for gene in genes:
        f.write("%s\n" % gene)
print(f"Number of Genes: {len(genes)}")