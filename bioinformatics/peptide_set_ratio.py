#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  7 11:14:32 2018

@author: hao
"""

import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np


def merge_peptides_dict(peptide_zs, current_peptide_zs_df, sample_num, sample_cnt):
    """Merge current_peptide_zs_df into existing peptide_zs as a new column sample_num"""
    occurrence_index = sample_cnt
    for index, row in current_peptide_zs_df.iterrows():
            #print(peptide_z)
            peptide_z = row['Peptide'] + '_' + str(row['z'])
            if peptide_z in peptide_zs.keys():
                peptide_zs[peptide_z][sample_num] = row['Intensity'] 
                peptide_zs[peptide_z][occurrence_index] += 1
            else:
                peptide_zs[peptide_z] = [0] * (sample_cnt + 1)
                peptide_zs[peptide_z][sample_num] = row['Intensity'] 
                peptide_zs[peptide_z][occurrence_index] = 1 
    return peptide_zs

def get_peptides_set(csv_path, file_name):
    """Go through all files in the directory and merge 
       the peptides in one set"""
    #peptide_zs is a dict peptide_z : [list of intensity in each samples]
    peptide_zs = {}
    dirs = os.listdir(csv_path)
    n = 8
    #n = 5
    sample_num = 0
    sample_names = []
    for sample_dir in dirs:        
        if sample_dir.startswith('.') or sample_dir.startswith("ratio") \
            or sample_dir.startswith('fig') or sample_dir.startswith('excluded'):
            continue
       
        #print(sample_dir)   
        sample_name = sample_dir.split('_')[0]
        sample_names.append(sample_name)
        
        csv_file = csv_path + "/" + sample_dir + "/" + file_name
        sample1 = pd.read_csv(csv_file)
        
        peptide_zs = merge_peptides_dict(peptide_zs, sample1, sample_num, n)
        """
        for index, row in sample1.iterrows():
            #print(peptide_z)
            peptide_z = row['Peptide'] + '_' + str(row['z'])
            if peptide_z in peptide_zs.keys():
                peptide_zs[peptide_z][sample_num] = row['Intensity'] 
                peptide_zs[peptide_z][n] += 1
            else:
                peptide_zs[peptide_z] = [0] * (n + 1)
                peptide_zs[peptide_z][sample_num] = row['Intensity'] 
                peptide_zs[peptide_z][n] = 1 
        """
        sample_num += 1
        
    print(len(peptide_zs))
    sample_names.append("#occurrence")
    peptide_zs_df = pd.DataFrame.from_dict(peptide_zs, orient='index', \
                                           columns=sample_names)
    return peptide_zs_df

def pair_order(ratios, label1, label2, csv_path):
    df = ratios[[label1, label2]].dropna()
    print(label1 + " and " + label2 + " shared " + str(len(df)) + " peptides")
    disorder_cnt = 0
    df = df.sort_values(label1)
    for index1 in range(0, len(df)):
        intensity1 = df.ix[index1, label2]
        for index2 in range(index1 + 1, len(df)):
            intensity2 = df.ix[index2, label2]
            if intensity1 > intensity2:
                disorder_cnt += 1
    total_cnt = len(df) * (len(df) - 1) / 2
    print("discorder ratio: " + str(disorder_cnt/total_cnt))
    title = label1 + ".vs." + label2 + ".log(ratio).highest intensity per peptide"
    #plot_ratio(df, title, csv_path)
    
    return disorder_cnt / total_cnt
    
def plot_ratio(ratios, title, csv_path):
    plot = ratios.plot(style = '.')
    fig = plot.get_figure()
    plt.xlabel("Peptide_z")
    plt.ylabel("Log(Sum Intensity / base intensity)")            
    fig_file_name = title + ".png"
    plt.title(title)    
    fig.savefig(csv_path + "/fig/" + fig_file_name)
    
def find_median_peptide(df):
    median_value = df.median()[0]
    for index, row in df.iterrows():
        if row[0] == median_value:
            break
        base_intensity = df.loc[index]
    return base_intensity

def main():
    csv_path = "/Users/hao/data/ABdata/extracted"
    file_name = "constructed_protein_Heavy.peptides.5.csv"
    peptides_zs = get_peptides_set(csv_path, file_name)  
    peptides_zs = peptides_zs.replace(0, np.nan)
    
    sample_num = 8
    #peptides_zs = peptides_zs.drop(columns = ['ab18034'])
    #peptides_zs = peptides_zs.drop(columns = ['ab18034', 'ab18045', 'ab18055'])
    #sample_num = 5
    
    shared_peptides_zs = peptides_zs.dropna()
    base = shared_peptides_zs.iloc[:, range(0, sample_num)]
    
    #normalize the intensity according to the intensity of the peptide with median intensity in sample 1
    base_intensity = find_median_peptide(base)
    
    #normalize the intensity according to the mean intensity of shared peptides
    #base = peptides_zs.loc[peptides_zs['#occurrence'] >= sample_num, :]
    #base_intensity = base.mean()
    
    
    shared_3 = peptides_zs.loc[peptides_zs['#occurrence'] > 2, :]
    ratios = shared_3.iloc[:, range(0, sample_num)] / base_intensity                                  
    #ratios = peptides_zs.iloc[:, range(0, sample_num)] / base_intensity
    #ratios = base / base_intensity
    ratios = np.log(ratios.astype(np.float64))
    
    for i in range(0, len(ratios.columns)):
        label1 = ratios.columns[i]
        for j in range((i + 1), len(ratios.columns)):
            label2 = ratios.columns[j]
            disorder_cnt = pair_order(ratios, label1, label2, csv_path)
            print(label1 + " vs " + label2 + " contains " + \
                  str(disorder_cnt) + " inconsistent order.")
    
    ratios = ratios.sort_values(['ab18040'])
    title = "5 sample.log(ratio)_using highest intensity for each peptide"
    plot_ratio(ratios, title, csv_path)

#shared_peptides = main()