#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 13:48:24 2018

@author: hao
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 14:17:57 2018

@author: hao
"""

import pandas as pd
import os

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



def build_peptide_dict(peptides):
    """Build peptide dictionary"""
    peptides_dict = {}
    """
    for peptide in peptides['Peptide'].tolist():
        peptides_dict[peptide[2:-2]] = peptide
    """
    
    for index, row in peptides.iterrows():
        peptide = row['Peptide']                      
        protein_accession = row['Protein Accession'].split("ups")[0]
        if protein_accession == 'P10636-8':
            protein_accession = 'P10636'
        peptides_dict[peptide[2:-2]] = (peptide, protein_accession)
    
    print("peptide dict: " + str(len(peptides_dict)))
    return peptides_dict

def read_peptide_feature(peptide_feature_file):
    peptide_features = pd.read_csv(peptide_feature_file)
    peptide_features = peptide_features[["z", "Intensity", 'DB seq']].dropna()
    peptide_features = peptide_features.rename(columns = {'DB seq' : 'Peptide'}) 
    
    print("peptide features in peptide_feature_file :" + str(len(peptide_features)))
    return peptide_features
     
    
"""Extract those features appears in unique peptides in protein-peptide table. 
Change the peptide sequence to the one with one AA before and behind
"""
def extract_peptide_z_intensity(peptides_dict, peptide_features):
    peptide_features = peptide_features.loc[map(lambda x : x in peptides_dict.keys(), \
                                              peptide_features['Peptide'])]
    peptide_features['Protein'] = peptide_features['z']
    for i in peptide_features.index: 
        (peptide, protein) = peptides_dict.get(peptide_features.loc[i, 'Peptide'])
        peptide_features.loc[i, 'Peptide'] = peptide
        peptide_features.loc[i, 'Protein'] = protein
          
    print("peptide feature number: " + str(len(peptide_features)))
    
    #Get the highest intensity for each charge 
    peptide_z_intensity = {}
    for i in peptide_features.index:
        key = peptide_features.loc[i, 'Peptide'] + ',' + str(peptide_features.loc[i, 'z'])
        if key in peptide_z_intensity.keys():
            if peptide_features.loc[i, 'Intensity'] > peptide_z_intensity[key]['Intensity']:
                peptide_z_intensity[key] = peptide_features.loc[i, ['Intensity', 'Protein']]
        else:
            peptide_z_intensity[key] = peptide_features.loc[i, ['Intensity', 'Protein']]
            
    print("peptide charge num: " + str(len(peptide_z_intensity)))
    return peptide_z_intensity

def export_result(out_csv_name, peptide_z_intensity):    
    with open(out_csv_name, 'w') as f:
        f.write("Peptide,z,Intensity,Protein\n")
        
        for key in peptide_z_intensity.keys():
            intensity = peptide_z_intensity.get(key)['Intensity']
            protein = peptide_z_intensity.get(key)['Protein']
            f.write(key + ',' + str(intensity) +',' + protein + '\n')
            

def extract_result(protein_peptide_file, peptide_feature_file_prefix, dir, sample):
    df = pd.read_csv(protein_peptide_file)

    data = df.loc[:, ['Protein Accession','Peptide', 'Unique', 'Intensity Sample 1', \
                      'Intensity Sample 2', 'Intensity Sample 3']]
    uniq_data = data.loc[data['Unique'] == 'Y']
    
    """Show the distribution of the intensity from three replicates"""
    #print(uniq_data.describe())
    #uniq_data.plot.box()
 
    rep_num = 3
    peptide_number = 0
    for i in range(1, rep_num + 1):        
        """export the peptides with intensity greater than zero"""
        rep_label = "Intensity Sample " + str(i)        
        rep = uniq_data.loc[uniq_data[rep_label] > 0, \
                             ['Protein Accession','Peptide', rep_label]].dropna()
                
        """rename the Intensity column"""
        rep_peptides = rep.rename(columns = {rep_label : "Intensity"})
        print("peptides number in " + rep_label + " is " + str(len(rep_peptides)))
        rep_peptides_dict = build_peptide_dict(rep_peptides)
        
        rep_peptide_feature_file = peptide_feature_file_prefix + str(i) + ".csv"
        rep_peptide_features = read_peptide_feature(rep_peptide_feature_file)
        
        rep_peptide_z_intensity = extract_peptide_z_intensity(rep_peptides_dict, rep_peptide_features)
        
        out_dir = "/Users/hao/data/UPSdata/extracted/" + sample
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
    
        out_csv_file = out_dir  + "/peptide_z_intensity." + str(i) + ".csv"
        export_result(out_csv_file, rep_peptide_z_intensity)
        
        peptide_number += len(rep_peptide_z_intensity)
        
        
    return peptide_number
 
#dir = "/Users/hao/data/UPSdata/JCox2014/"
#samples = ['JCox2014.ups1.ecoli_PEAKS_9', 'JCox2014.ups2.ecoli_PEAKS_9']
#dir = "/Users/hao/data/UPSdata/Giai.2016.ups1/" 
#samples = ['Giai.2016.10fmol_PEAKS_5', 'Giai.2016.25fmol_PEAKS_8', 'Giai.2016.5fmol_PEAKS_5']
dir = "/Users/hao/data/UPSdata/Arike2013/"
samples = ['Arike2013.B1.ups2.ecoli_PEAKS_8', 'Arike2013.B2.ups2.ecoli_PEAKS_8']

for sample in samples:
    protein_peptide_file =  dir + sample + "/protein-peptides.csv"
    peptide_feature_file_prefix = dir + sample + "/peptide features_"
    
    peptide_number = extract_result(protein_peptide_file, peptide_feature_file_prefix, dir, sample)
                                    
    print(peptide_number)
    
def read_AB_protein_peptides(main_protein, tryptic, protein_peptide_file):    
    protein_peptides = pd.read_csv(protein_peptide_file)  

    peptides = protein_peptides.loc[protein_peptides['Protein Accession'] == main_protein, \
                                          ['Protein Accession','Peptide', 'Unique', 'Area']]
    
    uniq_peptides = peptides.loc[peptides['Unique'] == 'Y'].dropna()
    full_tryptic_peptides = uniq_peptides.loc[map(is_full_tryptic, uniq_peptides['Peptide'])]
    
    if tryptic == False:
        peptides = uniq_peptides
    else:
        peptides = full_tryptic_peptides
    
    return peptides

"""
#main_protein = "constructed_protein_Light"
main_protein = "constructed_protein_Heavy"
            
for dir in ["ab18049_Lumos_v20180627_PROTEIN SEQUENCING_13", "ab18052_Lumos_v20180627_PROTEIN SEQUENCING_18", "ab18054_Lumos_v20180627_PROTEIN SEQUENCING_16"]:
    for tryptic in [False, True]:
        print(dir + " " + main_protein + " full-tryptic: " + str(tryptic))
        extract_abResult(tryptic, main_protein, dir)
         
        
    #peptide_feature.to_csv(out_csv_name, index = False)
"""