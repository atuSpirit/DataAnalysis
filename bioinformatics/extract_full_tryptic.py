#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 22 12:54:32 2018

@author: hao
"""

import pandas as pd

pd.set_option('display.max_rows',100)
pd.set_option('display.max_columns',10)
pd.set_option('display.max_colwidth',20)

"""
Judge whether a peptide is full tryptic or not
"""
def is_full_tryptic(peptide):
    if (peptide[0] == 'K' or peptide[0] == 'R') and \
    (peptide[-3] == 'K' or peptide[-3] == 'R') and \
    (peptide[-1] != 'P'):
        return True
    else:       
        return False

csv_file = "/Users/hao/data/Ramus/2500_rep2.csv"
csv_file = "/Users/hao/data/Ramus/protein-peptides.csv"
uniq_peptides = pd.read_csv(csv_file)
print(len(uniq_peptides))

full_tryptic_peptides = uniq_peptides.loc[map(is_full_tryptic, uniq_peptides['Peptide'])]
print(len(full_tryptic_peptides))
ratio = len(full_tryptic_peptides) / len(uniq_peptides)
print(ratio)