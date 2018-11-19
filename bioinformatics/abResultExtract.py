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

    
def extract_abResult(tryptic, main_protein, path, sample_name):
    protein_peptide_file = path + sample_name + "/protein-peptides.csv"
    protein_peptides = pd.read_csv(protein_peptide_file) 

    peptides_in_main_protein = protein_peptides.loc[protein_peptides['Protein Accession'] == main_protein, \
                                          ['Protein Accession','Peptide', 'Unique', 'Area']]
    uniq_peptides = peptides_in_main_protein.loc[peptides_in_main_protein['Unique'] == 'Y'].dropna()
    
    
    if tryptic == False:
        peptides = uniq_peptides
    else:
        full_tryptic_peptides = uniq_peptides.loc[map(is_full_tryptic, uniq_peptides['Peptide'])]
        peptides = full_tryptic_peptides
    
    """Build peptide dictionary"""
    peptides_dict = {}
    for peptide in peptides['Peptide'].tolist():
        peptides_dict[peptide[2:-2]] = peptide
    
    #print("peptide dict: " + str(len(peptides_dict)))
    
    for enzyme_num in range(1, 6):
        """Read in feature files according to enzyme_num"""
        peptide_feature_file = path + sample_name + "/peptide features_" + str(enzyme_num) + ".csv"
        peptide_feature = pd.read_csv(peptide_feature_file)
        peptide_feature = peptide_feature[["z", " Area Intensity", 'DB seq']].dropna()
        
        """Extract those features appears in unique peptides of selected protein in 
        protein-peptide table. Change the peptide sequence to the one with one AA before and behind
        """
        peptide_feature = peptide_feature.loc[map(lambda x : x in peptides_dict.keys(), \
                                                  peptide_feature['DB seq'])]
        for i in peptide_feature.index:           
            peptide_feature.loc[i, 'DB seq'] = peptides_dict.get(peptide_feature.loc[i, 'DB seq'])
        peptide_feature = peptide_feature.rename(columns = {' Area Intensity' : 'Intensity', \
                                                            'DB seq' : 'Peptide'}) 
        #print("peptide feature number: " + str(len(peptide_feature)))
        
        #Get the highest intensity for each charge 
        peptide_charge = {}
        for i in peptide_feature.index:
            key = peptide_feature.loc[i, 'Peptide'] + ',' + str(peptide_feature.loc[i, 'z'])
            if key in peptide_charge.keys():
                if peptide_feature.loc[i, 'Intensity'] > peptide_charge[key]:
                    peptide_charge[key] = peptide_feature.loc[i, 'Intensity']
            else:
                peptide_charge[key] = peptide_feature.loc[i, 'Intensity']
                
        print("peptide charge num: " + str(len(peptide_charge)))
        
        
        """Export the results to csv file"""
        out_dir = path + "extracted/" + sample_dir
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
                
        if tryptic == True:
            out_csv_name = out_dir + "/" + main_protein + ".full_tryptic." + str(enzyme_num) + ".csv"
        else:
            #out_csv_name = "/Users/hao/data/ABdata/extracted/" + out_dir + "/peptides." + str(enzyme_num) + ".csv"
            out_csv_name = out_dir + "/" + main_protein + ".peptides." + str(enzyme_num) + ".csv"
        
        with open(out_csv_name, 'w') as f:
            f.write("Peptide,z,Intensity\n")
            
            for key, intensity in peptide_charge.items():
                f.write(key + ',' + str(intensity) + '\n')
    


proteins = ["constructed_protein_Heavy", "constructed_protein_Light"]
path = '/Users/hao/data/ABdata/samples/'  
for sample_dir in os.listdir(path):
    if sample_dir.startswith('.') or sample_dir == "extracted":
        continue
    print("Dealing with " + sample_dir)    
    for protein in proteins:
#        for tryptic in [False, True]:                
        tryptic = False
        extract_abResult(tryptic, protein, path, sample_dir)
     

        
    #peptide_feature.to_csv(out_csv_name, index = False)
