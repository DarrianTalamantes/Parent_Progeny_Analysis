#################### Depreciated: Found out I can make numerical genotypes in Tassel.
#################### Main problem with this is that i am changing refrence for each population
import numpy as np
import pandas as pd


def annova_testing():
    progeny_hmp = importHMP(
        "~/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Filtered_data/FlexSeq_filtered_Trans.hmp")
    parent_hmp = importHMP(
        "~/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Filtered_data/parents1500_all.hmp.txt")
    # print(progeny_hmp)
    print(parent_hmp)
    data301 = parent_data_frame_maker(progeny_hmp, str(301))
    print(data301)
    score301 = scoring_data_frames(parent_hmp, data301)
    print(score301)

def importHMP(file):
    data = pd.read_csv(file, sep='\t', header=0, na_values=" ")
    data.drop(data.columns[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]], axis=1, inplace=True)
    return data


def parent_data_frame_maker(progeny_hmp, parentnum):
    new = []
    for s in progeny_hmp.columns:
        if s.startswith(parentnum):
            adder = progeny_hmp[s]
            # new.insert(1, adder)
            new.append(adder)
    df = pd.DataFrame(new)
    return df


def scoring_data_frames(parent_hmp, dataframe):
    #  Parent and progeny data tables are inversed
    parent_array = parent_hmp.to_numpy()
    progeny_array = dataframe.to_numpy()
    score_array = np.zeros((len(progeny_array), len(progeny_array[0])), dtype=int)
    for parent_idv in range(len(parent_array[0])):  # number of parents 17
        for snp_site in range(len(parent_array)):  # number of sites 1213
            nuc_parent = parent_array[snp_site][parent_idv]
            for progeny_idv in range(len(progeny_array)):  # number of progeny 200
                nuc_progeny = dataframe[progeny_idv][snp_site]
                if nuc_progeny == nuc_parent:
                    point = 2
                elif nuc_progeny != nuc_parent:
                    print(nuc_progeny, nuc_parent)
                    point = oneorzero(nuc_progeny, nuc_parent)
                score_array[progeny_idv][snp_site] = point
    return score_array


def oneorzero(nuc_progeny, nuc_parent):
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
        point = 1
    elif small_array_parent[1] == small_array_progeny[1] or small_array_parent[1] == small_array_progeny[0]:
        point = 1
    return point

annova_testing()
