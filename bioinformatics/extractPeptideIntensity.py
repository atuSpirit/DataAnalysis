# -*- coding: utf-8 -*-
"""
Created on Wed Aug  8 10:38:14 2018
Read in peaks protein_peptide.csv and export peptide sequence and intensity
@author: Hao Lin
"""
import pandas as pd

pd.set_option('display.max_rows',100)
pd.set_option('display.max_columns',10)
pd.set_option('display.max_colwidth',20)


def extract_result(sample_name):
    csv_file_name = "D:/hao/result/UPS1/pxd001819_ramus2016/UPS1." + sample_name + "mol.yeast_PEAKS_5/protein-peptides.csv"    
    #csv_file_name = "D:/hao/result/UPS1/ups1.Giai2016/PXD002370_Giai2016." + sample_name + "fmol_PEAKS_5/protein-peptides.csv"    
    print(csv_file_name)

    df = pd.read_csv(csv_file_name)

    data = df.loc[:, ['Protein Accession','Peptide', 'Unique', 'Intensity Sample 1', 'Intensity Sample 2', 'Intensity Sample 3']]
    uniq_data = data.loc[data['Unique'] == 'Y']
    
    """Show the distribution of the intensity from three replicates"""
    #print(uniq_data.describe())
    #uniq_data.plot.box()
    
    """dump the results to three files respectively"""
    extract_labels = ['Intensity Sample 1', 'Intensity Sample 2', 'Intensity Sample 3']
    i = 1
    peptide_number = 0
    for label in extract_labels:
        """export the peptides with intensity greater than zero"""
        rep_label = "Intensity Sample " + str(i)        
        rep = uniq_data.loc[uniq_data[rep_label] > 0, ['Protein Accession','Peptide', rep_label]].dropna()
        
        """rename the Intensity column"""
        rep = rep.rename(columns = {rep_label : "Intensity"})
        print(rep.head(2))
        
        peptide_number += rep.size / 3   # 3 is the column numbers
        print(rep.size / 3)
        
        if (rep.size > 0):            
            #out_csv_name = "D:/hao/result/UPS1/pxd001819_ramus2016/new_extracted/" + sample_name + "_rep" + str(i) + ".csv" 
            out_csv_name = "D:/hao/result/UPS1/ups1.Giai2016/extracted/" + sample_name + "_rep" + str(i) + ".csv" 
            rep.to_csv(out_csv_name)
            
        i += 1
           
        
        
    return peptide_number
 

samples = ["250", "500", "2500", "5000", "12500", "25000","50000"]
#samples = ["25", "10", "12.5", "5"]
total_peptides_number = 0
for sample_name in samples:
    total_peptides_number += extract_result(sample_name)
print("Total peptide number is :" + str(total_peptides_number))