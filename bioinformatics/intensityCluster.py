#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 13:15:42 2018

@author: hao
"""
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.cluster import KMeans
import numpy as np


pd.set_option('display.max_rows',100)
pd.set_option('display.max_columns',10)
pd.set_option('display.max_colwidth',20)

csv_file = "/Users/hao/data/Nuno/nuno.HC.vs.2HC/Nuno2016.HC_PEAKS_16/protein-peptides.csv"
#csv_file = "/Users/hao/data/Nuno/nuno.LC.vs.HCandLC/Nuno2016.LC_PEAKS_10/protein-peptides.csv"
#csv_file = "/Users/hao/data/Nuno/nuno.HC.vs.HC/Nuno2016.HC_PEAKS_22/protein-peptides.csv"
#csv_file = "/Users/hao/data/Nuno/Nuno.2016.HC.abdb/Nuno2016.HC_PEAKS_26/protein-peptides.csv"
proteins = pd.read_csv(csv_file)

#Partition protein-peptides file according to different enzyme
proteinsPepByEnzyme = []
for i in range(1, 9):
    proteinPeps = proteins[['Protein Group', 'Peptide', 'Unique','Intensity Sample ' + str(i)]].dropna().rename(columns = {'Intensity Sample ' + str(i):'Intensity'})    
    proteinPeps = proteinPeps.loc[proteinPeps['Intensity'] >0, :]
    proteinsPepByEnzyme.append(proteinPeps)


#extract uniq
def intensityStatistic(proteinPepTable, i):    
    
    uniq_proteinPep = proteinPepTable.loc[proteinPepTable['Unique'] == 'Y', :]
    print(str(len(proteinPepTable)) + '\t' + str(len(uniq_proteinPep)))
    
    #uniq_proteinPep.boxplot(by = 'Protein Group')
    #proteinPepTable.boxplot(by = 'Prot ein Group')
    
    #savefig("enzyme " + str(i) + ".png")
    #uniq_proteinPep.hist(by = 'Protein Group')
    
    
    
    #clusterAccordingIntensity(uniq_proteinPep['Intensity'].tolist().reshape(-1,1))
    
    """
    hc1 = uniq_proteinPep.loc[uniq_proteinPep['Protein Group'] == 1, 'Intensity']
    hc2 = uniq_proteinPep.loc[uniq_proteinPep['Protein Group'] == 2, 'Intensity']
    plt.hist(hc1, color = 'r')
    plt.hist(hc2, stacked = True, color = 'b')
    plt.show()
    """
    
    #print(hc1.describe())
    #print(hc2.describe())
    #plt.scatter(range(len(hc1)), hc1, color = 'r')
    #plt.scatter(range(len(hc2)), hc2, color = 'b')
    

def clusterAccordingIntensity(intensity):
    kmeans = KMeans(n_clusters=2, random_state=0).fit(intensity)
    print(kmeans.labels_)
    
   
for i in range(8):
    proteinPepTable = proteinsPepByEnzyme[i]    
    intensityStatistic(proteinPepTable, i)
    sortedProteinPepTable = proteinPepTable.sort_values('Intensity')
    sortedProteinPepTable['Intensity'].hist()
    


