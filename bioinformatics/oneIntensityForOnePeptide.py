#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 12:56:25 2018
Read in extracted peptide_z intensity file
Extract intensity for each pepide
@author: hao
"""


import pandas as pd
import numpy as np
import os
from peptide_set_ratio import plot_ratio, pair_order

def read_peptide_z_intensity(filename):
    peptides = {}
    data = pd.read_csv(filename)
    for index, row in data.iterrows():
            peptide = row[0]
            z = row[1]
            intensity = row[2]
            if peptide in peptides.keys():
                peptides[peptide].append((z, intensity))
            else:
                peptides[peptide] = [(z, intensity)]
                
    return peptides

def choose_higher_intensity(peptides):
    peptides_with_higher_intensity = {}
    for peptide, lst in peptides.items(): 
        highest_intensity = 0
        highest_tuple = (0, 0)
        if len(lst) > 1:
            for item in lst:
                intensity = item[1]
                if intensity > highest_intensity:
                    highest_intensity = intensity
                    highest_tuple = item
            peptides_with_higher_intensity[peptide] = highest_tuple
        else:
            peptides_with_higher_intensity[peptide] = lst[0]
                
    return peptides_with_higher_intensity


def dict_to_dataframe(peptides):
    columns = ['Peptide', 'z', 'Intensity']
    df = pd.DataFrame(columns=columns)
    for peptide, intensity in peptides.items():
        df1 = pd.DataFrame([[peptide, intensity[0], intensity[1]]], index=[peptide],columns=columns)
        df = df.append(df1)
                
    return df
        

csv_path = "/Users/hao/data/ABdata/extracted"
file_name = "constructed_protein_Heavy.peptides.5.csv"
dirs = os.listdir(csv_path)
sample_names = []
sample_peptides = {}

for sample_dir in dirs:        
    if sample_dir.startswith('.') or sample_dir.startswith("ratio") \
            or sample_dir.startswith('fig') or sample_dir.startswith('excluded'):
        continue
    sample_name = sample_dir.split('_')[0]
    sample_names.append(sample_name)
        
    file_path = csv_path + '/' + sample_dir +'/' + file_name
    
    peptides = read_peptide_z_intensity(file_path)
    peptides_dict = choose_higher_intensity(peptides)
    
    sample_peptides[sample_name] = dict_to_dataframe(peptides_dict)  
    
    
for i in range(0, len(sample_names)):
    sample1 = sample_names[i]
    for j in range(i + 1, len(sample_names)):
        sample2 = sample_names[j]
        share2 = pd.merge(sample_peptides[sample1], sample_peptides[sample2], \
                          how='inner', on=['Peptide', 'z'])
        #share2 = share2.loc[share2['z_x'] == share2['z_y']]
        #share2 = share2.drop(columns=['z_x', 'z_y'])
        share2 = share2.set_index('Peptide')
        #base_intensity = find_median_peptide(share2)
        base_intensity = share2.loc['K.NTQPIMDTDGSYFVYSK.L']
        ratios = share2 / base_intensity
        ratios = ratios.rename(columns={'Intensity_x' : sample1, 'Intensity_y':sample2})
        ratios = np.log(ratios.astype(np.float64))
        ratios = ratios.sort_values(sample1)
        title = sample1 + ".vs." + sample2 + ".log(ratio).higher.png"
        
        #plot_ratio(ratios, title, csv_path)
        
        pair_order(ratios, sample1, sample2, csv_path)

peptide_zs = sample_peptides[sample_names[0]]
previous_sample = sample_names[0]
for sample_name in sample_names:
    if sample_name == sample_names[0]:
        continue
    peptide_zs = pd.merge(peptide_zs, sample_peptides[sample_name], \
                          how='inner', on=['Peptide', 'z'])
    print(peptide_zs.columns.values)
    peptide_zs = peptide_zs.rename(columns={'Intensity_x': previous_sample, \
                                            'Intensity':sample_name, 'Intensity_y':sample_name})
    print(peptide_zs.columns.values)
    previous_sample = sample_name

peptide_zs['Peptide_z'] = peptide_zs['Peptide'] + '_' + peptide_zs['z'].map(str)
peptide_zs = peptide_zs.set_index('Peptide')
  
peptide_zs = peptide_zs.drop(columns=['z'])
base_intensity = peptide_zs.loc['K.NTQPIMDTDGSYFVYSK.L']  
peptide_zs = peptide_zs.set_index('Peptide_z')
shared_ratios = peptide_zs / base_intensity.iloc[0:5]
shared_ratios = np.log(shared_ratios.astype(np.float64))
shared_ratios = shared_ratios.sort_values('ab18040')
title = "5 sample.log(ratio)_using highest intensity for each peptide"
plot_ratio(shared_ratios, title, csv_path)
