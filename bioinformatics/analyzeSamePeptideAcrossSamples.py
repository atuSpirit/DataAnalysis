# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 10:21:28 2018
For several peptide intensity csv file, merge them with peptide sequence
@author: Hao Lin
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np



pd.set_option('display.max_rows',100)
pd.set_option('display.max_columns',10)
pd.set_option('display.max_colwidth',40)

sample_num = 3

dir = '/Users/hao/data/ABdata/extracted/'
#dir = '/Users/hao/data/UPSdata/extracted/'

file1 = dir + 'ab18049_Lumos_v20180627_PROTEIN SEQUENCING_13/constructed_protein_Heavy.peptides.5.csv'
#file1 = dir + 'ab18054_Lumos_v20180627_PROTEIN SEQUENCING_16/constructed_protein_Heavy.peptides.5.csv'
#file1 = dir + 'Giai.2016.10fmol_PEAKS_5/peptide_z_intensity.1.csv'
#file1 = dir + 'Arike2013.B1.ups2.ecoli_PEAKS_8/peptide_z_intensity.1.csv'

data1 = pd.read_csv(file1)

file2 = dir + 'ab18052_Lumos_v20180627_PROTEIN SEQUENCING_18/constructed_protein_Heavy.peptides.5.csv'
#file2 = dir + 'ab18054_Lumos_v20180627_PROTEIN SEQUENCING_16/constructed_protein_Light.peptides.5.csv'
#file2 = dir + 'PEAKS_AB_Demo_LC_PEAKS_13/' + 'peptides.csv'
#file2 = dir + 'Arike2013.B1.ups2.ecoli_PEAKS_8/peptide_z_intensity.1.csv'
#file2 = dir + 'Giai.2016.5fmol_PEAKS_5/peptide_z_intensity.1.csv'
#file2 = dir + 'JCox2014.ups1.ecoli_PEAKS_9/peptide_z_intensity.1.csv'
#file2 = dir + 'Arike2013.B2.ups2.ecoli_PEAKS_8/peptide_z_intensity.1.csv'
#file2 = dir + 'JCox2014.ups2.ecoli_PEAKS_9/peptide_z_intensity.1.csv'
data2 = pd.read_csv(file2)

if sample_num == 3:    
    file3 = dir + 'ab18054_Lumos_v20180627_PROTEIN SEQUENCING_16/constructed_protein_Heavy.peptides.5.csv'
    #file3 = dir + 'Arike2013.B2.ups2.ecoli_PEAKS_8/peptide_z_intensity.1.csv'
    #file3 = dir + 'Giai.2016.25fmol_PEAKS_8/peptide_z_intensity.1.csv'
    data3 = pd.read_csv(file3)

"""
tic = [1.0690025e13, 2.8862634e12, 1.43400615E13]
ratio = [tic[0] / min(tic), tic[1] / min(tic), tic[2] / min(tic)]
data1['Intensity'] /= ratio[0]
data2['Intensity'] /= ratio[1]
data3['Intensity'] /= ratio[2]
"""

plt.hist(data1['Intensity'], alpha=0.5, color = 'r')
plt.hist(data2['Intensity'], alpha=0.5, stacked = True, color = 'b')
if sample_num == 3:
    plt.hist(data3['Intensity'], alpha = 0.5, stacked = True, color = 'y')
plt.show()

merged = pd.merge(data1, data2, on = ["Peptide", 'z'], how = "outer").dropna()
if sample_num == 3:
    merged = pd.merge(merged, data3, on = ["Peptide", 'z'], how = "outer").dropna()

#merged = merged.rename(columns = {'Intensity_x' : 'ab18049', 'Intensity_y' : 'ab18052', \
#                         'Intensity' : 'ab18054'})

merged['Peptide'] = merged['Peptide'] + "_" + merged['z'].map(str)
if sample_num == 2:
    merged_intensities = merged[['Peptide', 'Intensity_x', 'Intensity_y']]
else:
    merged_intensities = merged[['Peptide', 'Intensity_x', 'Intensity_y', 'Intensity']]
#merged_intensities = merged[['Peptide', 'ab18049', 'ab18052', 'ab18054']]
merged_intensities = merged_intensities.rename(columns = {'Peptide' : 'Peptide_z'})

sorted_merge = merged_intensities.sort_values(['Intensity_x'])

median_value = sorted_merge.median()[0]
for index, row in sorted_merge.iterrows():
    if row['Intensity_x'] == median_value:
        break

if sample_num == 2:
    base_intensity = sorted_merge.iloc[index, [1, 2]]
else:
    base_intensity = sorted_merge.iloc[index, [1, 2, 3]]

sorted_merge = sorted_merge.replace(0, np.nan)

if sample_num == 2:
    ratios = sorted_merge.iloc[:, [1, 2]] / base_intensity
else:
    ratios = sorted_merge.iloc[:, [1, 2, 3]] / base_intensity
    
ratios = np.log(ratios.astype(np.float64))

ratios['Peptide_z'] = sorted_merge['Peptide_z']
ratios = ratios.sort_values(['Intensity_x'])
ratios = ratios.set_index('Peptide_z')
#ratios = ratios.rename(index = str, columns = {'Intensity_x' : 'Giai.10fmol', 'Intensity_y' : 'JCox'})
#ratios = ratios.rename(index = str, columns = {'Intensity_x' : '10fmol', 'Intensity_y' : '5fmol', 'Intensity' : '25fmol'})
ratios = ratios.rename(index = str, columns = {'Intensity_x' : 'ab18049', 'Intensity_y' : 'ab18052', 'Intensity' : 'ab18054'})

#ratios = ratios.rename(index = str, columns = {'Intensity_x' : 'ab18049', \
#                                               'Intensity_y' : 'ab18052'})
#ratios = ratios.rename(index = str, columns = {'Intensity_x' : 'Arike.ups2.B1', \
#                                               'Intensity_y' : 'Arike.ups2.B2'})
#ratios = ratios.loc[ratios['ab18049'] < 20, :]

#title = "ratio.heavy chain peptide ratios"
#title = "ratio.light.ab18049.vs.ab18052"
#title = "ratio.heavy.ab18049.vs.ab18054"
#title = "Giai 5fmol vs 10fmol vs 25fmol peptides ratio"
#title = "ratio.UPS1.JCox.vs.Giai.10fmol"
#title = "ratio.UPS2.Arike.vs.JCox"
#title = "ratio.UPS2.Arike.B1.vs.B2"
#title = "ratio.light.ab18049.vs.ab18054"
#title = "ratio.light.ab18049.vs.ab18052.vs.ab18054"
title = "ratio.heavy.ab18049.vs.ab18052.vs.ab18054"

fig_file_name = title + ".png"
#fig_file_name = "ratio.Giai.5fmol.vs.10fmol.vs.25fmol.png"

#merged_intensities = merged_intensities.set_index(merged_intensities['Peptide_z']) #,drop = True, inplace = True)

#merged.to_csv("D:/hao/result/UPS1/pxd001819_ramus2016/new_extracted/" + sample_name1 + "vs." + sample_name2 +".csv")

print(len(merged))
if len(merged) > 0:
    #merged.plot.scatter(x = "Intensity_x", y = "Intensity_y")    
    #plot = merged_intensities.plot()    
    plot = ratios.plot(style = '.')
    fig = plot.get_figure()
    plt.xlabel("Peptide_z")
    plt.ylabel("Ratio")
    #plt.title("ab heavy peptides")
    plt.title(title)    
    fig.savefig(dir + fig_file_name)
    #plt.title("heavy chain full tryptic peptides")
    #fig.savefig(dir + "heavy_full_tryptic_peptides.png")