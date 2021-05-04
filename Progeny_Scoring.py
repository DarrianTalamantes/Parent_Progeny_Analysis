from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import time
import errno
import os


def progenyScoring():
    parent_hmp = importHMP(
        "~/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Filtered_data/parents1500_all.hmp.txt")
    # print(parent_hmp)

    progeny_hmp = importHMP(
        "~/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Filtered_data/FlexSeq_filtered.hmp")
    # print(progeny_hmp)

    parent_array = parent_hmp.to_numpy()
    progeny_array = progeny_hmp.to_numpy()
    get_score_array(parent_array, progeny_array)


def scoring(nuc_parent, nuc_progeny):
    switcher = {
        'A': np.array(["A", "A"]),
        'C': np.array(["C", "C"]),
        'G': np.array(["G", "G"]),
        'T': np.array(["T", "T"]),
        'R': np.array(["A", "G"]),
        'Y': np.array(["C", "T"]),
        'S': np.array(["G", "C"]),
        'W': np.array(["A", "T"]),
        'K': np.array(["G", "T"]),
        'M': np.array(["A", "C"])
    }
    small_array_parent = switcher.get(nuc_parent)
    small_array_progeny = switcher.get(nuc_progeny)
    point = 0
    # print(nuc_parent, nuc_progeny)
    if small_array_parent[0] == small_array_progeny[0] or small_array_parent[0] == small_array_progeny[1]:
        point += 1
    elif small_array_parent[1] == small_array_progeny[1] or small_array_parent[1] == small_array_progeny[0]:
        point += 1

    return point


def get_score_array(parent, progeny):
    ########################################### Stuff here is for testing purposes
    # partial_parents = np.zeros((len(parent), 3), dtype=str)
    # for i in range(3):
    #     for j in range(len(parent)):
    #         partial_parents[j][i] = parent[j][i]
    # partial_progeny = np.zeros((len(progeny), 10), dtype=str)
    # for i in range(10):
    #     for j in range(len(progeny)):
    #         partial_progeny[j][i] = progeny[j][i]
    # parent = partial_parents
    # progeny = partial_progeny
    #########################################
    # This creates an array of zeroes that we will use to add up the scores
    scoreCols = len(parent[0])
    scoreRows = len(progeny[0])
    score_array = np.zeros((scoreRows, scoreCols), dtype=int)
    max_score_array = np.zeros((scoreRows, scoreCols), dtype=int)
    nuc_proj_array = np.zeros((len(progeny), len(progeny[0])), dtype=str)
    nuc_parent_array = np.zeros((len(parent), len(parent[0])), dtype=str)
    print(score_array.shape)
    print(parent.shape)
    print(progeny.shape)
    # print(parent)
    # print(progeny)
    # This will compare the parents and progeny sites to find which are the most similar
    parent_idv = 0
    snp_site = 0
    progeny_idv = 0
    for parent_idv in range(len(parent[0])):  # number of parents 17
        for snp_site in range(len(parent)):  # number of sites 1213
            nuc_parent = parent[snp_site][parent_idv]
            # print(snp_site, parent_idv, "parent site") # For sanity check
            # nuc_parent_array[snp_site][parent_idv] = nuc_parent # For sanity check
            for progeny_idv in range(len(progeny[0])):  # number of progeny 3000
                nuc_progeny = progeny[snp_site][progeny_idv]
                # print(snp_site, progeny_idv, "progeny site") # For sanity check
                # nuc_proj_array[snp_site][progeny_idv] = nuc_progeny
                # If there are any Ns or the parents are heterozygous no points added to max score
                if nuc_progeny != "N" and nuc_parent != "N":
                    # if nuc_parent == "A" or nuc_parent == "C" or nuc_parent == "G" or nuc_parent == "T":
                    #     if nuc_progeny == "A" or nuc_progeny == "C" or nuc_progeny == "G" or nuc_progeny == "T":
                    x = scoring(nuc_parent, nuc_progeny)
                    max_score_array[progeny_idv][parent_idv] += 1
                    score_array[progeny_idv][parent_idv] += x
        # input("Press Enter to continue...") # For sanity check
    print("Score Array")
    print(score_array)
    # Max score array will nto differ within a progeny because there are no Ns in parents
    norm_score_array = np.divide(score_array, max_score_array)
    print("Normal Score Array")
    print(norm_score_array)
    print("Max Score Array")
    print(max_score_array)
    # print(nuc_proj_array) # For sanity check
    # print(nuc_parent_array) # For Sanity check
    # Saving stuff to files to later use in R
    np.savetxt("/home/drt83172/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Filtered_data"
               "/parent_progeny_score_table_normalized.txt", norm_score_array)
    np.savetxt("/home/drt83172/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Filtered_data"
               "/parent_progeny_score_table.txt", score_array)
    np.savetxt("/home/drt83172/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Filtered_data"
               "/max_scores.txt", max_score_array)


def importHMP(file):
    data = pd.read_csv(file, sep='\t', header=0, na_values=" ")
    data.drop(data.columns[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]], axis=1, inplace=True)
    return data


progenyScoring()

# something is making my max score file be the same for every parent per progeny