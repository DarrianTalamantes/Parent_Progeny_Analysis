from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import errno
import os


def progenyScoring():
    parent_hmp = importHMP(
        "~/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Filtered_data/parents1500_all.hmp.txt")
    # print(parent_hmp)

    progeny_hmp = importHMP(
        "~/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Filtered_data/Good_Data.hmp.txt")
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
    # This creates an array of zeroes that we will use to add up the scores
    scoreCols = len(parent[0])
    scoreRows = len(progeny[0])
    score_array = np.zeros((scoreRows, scoreCols), dtype=int)
    max_score_array = np.zeros((scoreRows, scoreCols), dtype=int)
    print(score_array.shape)
    print(parent.shape)
    # print(parent)
    print(progeny.shape)
    # print(progeny)

    # This will compare the parents and progeny sites to find which are the most similar
    parent_idv = 0
    snp_site = 0
    progeny_idv = 0
    while parent_idv < len(parent[0]):  # number of parents 17
        while snp_site < len(parent):  # number of sites 1213
            nuc_parent = parent[snp_site][parent_idv]
            while progeny_idv < len(progeny[0]):  # number of progeny 3000
                nuc_progeny = progeny[snp_site][progeny_idv]
                if nuc_progeny != "N" and nuc_parent != "N":
                    # x = scoring(nuc_parent, nuc_progeny)
                    max_score_array[progeny_idv][parent_idv] += 1
                    score_array[progeny_idv][parent_idv] += 1
                progeny_idv = progeny_idv + 1
            snp_site = snp_site + 1
            progeny_idv = 0
        parent_idv += 1
        snp_site = 0
    print(score_array)
    norm_score_array = np.divide(score_array, max_score_array)
    print(norm_score_array)
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
