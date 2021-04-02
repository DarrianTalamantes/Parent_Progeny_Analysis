from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import errno
import os


def siteAnalysis():
    print(f'Hi, What is your name?')
    parent_hmp = importHMP(
        "~/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Filtered_data/parents3000.hmp.txt")
    groupBySite(parent_hmp)


# takes a hmp file and strips it of inimportant data.
def importHMP(file):
    data = pd.read_csv(file, sep='\t', header=0, na_values=" ")
    data.drop(data.columns[[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]], axis=1, inplace=True)
    return data


# Takes importated hmp file and gives all sites a score for each nucleotide, high of 32, low of 0
def groupBySite(hmp):
    counts = hmp[['rs#']]
    hmp = hmp.set_index('rs#')
    # print(counts)
    # print(hmp)
    total_row = len(hmp.index)
    total_col = len(hmp.columns)
    print(total_row)
    print(total_col)
    # This array will be used to sum up the nucleotides in all sites
    array_of_sites = np.zeros((total_row, 4))
    print(array_of_sites)
    i = 0
    j = 0
    x = 0
    hmp_array = hmp.to_numpy()
    print(hmp_array.shape)

    while j < total_row:
        # print("break")
        while i < total_col:
            iupac = hmp_array[j, i]
            array_to_add = iupacNums(iupac)
            while x < 4:
                # print(array_to_add[x])
                array_of_sites[j, x] = array_of_sites[j][x] + array_to_add[x]
                x = x + 1
            x = 0
            # print(iupac)
            i = i + 1
        j = j + 1
        i = 0
    print(array_of_sites)


# This function is a library to make iupac symbols into numerical data
def iupacNums(argument):
    switcher = {
        # array is [A,C,G,T]
        'A': [2, 0, 0, 0],
        'C': [0, 2, 0, 0],
        'G': [0, 0, 2, 0],
        'T': [0, 0, 0, 2],
        'R': [1, 0, 1, 0],
        'Y': [0, 1, 0, 1],
        'S': [0, 1, 1, 0],
        'W': [1, 0, 0, 1],
        'K': [0, 0, 1, 1],
        'M': [1, 1, 0, 0]

    }
    # Get the function from switcher dictionary
    small_array = switcher.get(argument)
    # Execute the function
    # print(type(small_array))
    return small_array


siteAnalysis()
