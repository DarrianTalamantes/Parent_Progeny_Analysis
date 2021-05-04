import numpy as np
import pandas as pd


def annova_testing():
    progeny_hmp = importHMP(
        "~/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Filtered_data/FlexSeq_filtered_Trans.hmp")
    parent_hmp = importHMP(
        "~/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/Filtered_data/parents1500_all.hmp.txt")
    # print(progeny_hmp)
    print(parent_hmp)
    new = []
    for s in progeny_hmp.columns:
        if s.startswith('301'):
            adder = progeny_hmp[s]
            # new.insert(1, adder)
            new.append(adder)
    df = pd.DataFrame(new)
    print(df)



    # print(len(new))
    # print(len(new[0]))



def importHMP(file):
    data = pd.read_csv(file, sep='\t', header=0, na_values=" ")
    data.drop(data.columns[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]], axis=1, inplace=True)
    return data

annova_testing()
