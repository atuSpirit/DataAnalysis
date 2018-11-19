#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 16:37:12 2018

@author: hao
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

pd.set_option('display.max_rows',100)
pd.set_option('display.max_columns',10)
pd.set_option('display.max_colwidth',40)

dir = "/Users/hao/data/ABdata/extracted/"
file1 = dir + 'ab18049_Lumos_v20180627_PROTEIN_SEQUENCING_13/light_chain/' + 'peptides.csv'
file2 = dir + 'ab18052_Lumos_v20180627_PROTEIN_SEQUENCING_18/light_chain/' + 'peptides.csv'
file3 = dir + 'ab18054_Lumos_v20180627_PROTEIN_SEQUENCING_16/light_chain/' + 'peptides.csv'

data1 = pd.read_csv(file1)
data2 = pd.read_csv(file2)
data3 = pd.read_csv(file3)

intensity1 = data1['Intensity'].dropna()
intensity2 = data2['Intensity'].dropna()
intensity3 = data3['Intensity'].dropna()

#bins=np.histogram(np.hstack((intensity1,intensity2)), bins = 40)[1]


#plt.hist(intensity1, bins, color = 'r')
#plt.hist(intensity2, bins, color = 'b')

#plt.boxplot(intensity1)
#plt.boxplot(intensity2)

print(intensity1.sum())
print(intensity2.sum())
print(intensity3.sum())
