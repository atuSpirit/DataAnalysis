#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 11:05:14 2018
Compute the theoretic peptide number given protein database, enzyme and number
of allowed miscleavage
@author: Hao Lin
"""

def importDatabase(databasePath):
    proteins = []
    
    with open(databasePath) as f_obj:
        for line in f_obj:
            if not line.startswith('>'):                
                proteins.append(line.strip())    
                
    return proteins


def cutByPattern(protein, patternList, excludePattern):
    """ Cut protein by characters in patternList if they are not followed 
    by excludePattern. For trypsin, patternList is [K, R], excludePattern is P.
    For LysC, patternList is K, excludePattern is ''. """
    peptideList = []
    peptide = ''
    for index in range(len(protein)):
        flag = False
        peptide += protein[index]
        
        for pattern in patternList:
            if (protein[index] == pattern):
                if index < (len(protein) - 1):
                    if protein[index + 1] != excludePattern:
                        flag = True
                        break
                else:
                    flag = True
                    break
        
        if (flag and len(peptide) > 0):
            peptideList.append(peptide)
            peptide = '' 
    if len(peptide) > 0:                
        peptideList.append(peptide)
       
    return peptideList

def mergeMiscleavage(peptideList, miscleave_num):
    """Allowing miscleavage will allow peptides not be cut. So concantenate i 
    peptides i = 1, 2, ..., miscleave_num peptides will generate a new
    peptide list allowing miscleavege"""
    miscleaved_peptides = []
    for peptide_index in range(len(peptideList)):
        newPeptide = peptideList[peptide_index]
        
        miscleaved_peptides.append(newPeptide)           
        i = 1
        while i <= miscleave_num:
            if peptide_index + i < len(peptideList):
                newPeptide += peptideList[peptide_index + i]
                miscleaved_peptides.append(newPeptide)
            i += 1
           
    return miscleaved_peptides

def peptideLenFilter(peptideList, min_len):
    filteredPeptides = []
    for peptide in peptideList:
        if len(peptide) >= min_len:   
            filteredPeptides.append(peptide)
    return filteredPeptides
        
def uniquePeptideFilter(peptideList):
    filteredPeptides = []
    uniquePeptides = {}
    for peptide in peptideList:
        if peptide in uniquePeptides.keys():
            uniquePeptides[peptide] += 1  
            #print(peptide + ": ")
            #print(uniquePeptides[peptide])
        else:
            uniquePeptides[peptide] = 1
    
    for peptide in uniquePeptides.keys():
        if uniquePeptides[peptide] == 1:
            filteredPeptides.append(peptide)
            
            
    return filteredPeptides
        

def statisticOfUniqPeptideBeCut(fastaFile, enzymeName, miscleave_num, peptide_min_len):
    proteins = importDatabase(fastaFile)
    peptides = []
    for protein in proteins:
        peptides += cutByPattern(protein, cutPatterns[enzymeName], excludePatterns[enzymeName])
    
    miscleaved_peptides = mergeMiscleavage(peptides, miscleave_num)
    filteredPeptides = peptideLenFilter(miscleaved_peptides, peptide_min_len)
    uniquePeptides = uniquePeptideFilter(filteredPeptides)
    print('\tenzymeName\tmiscleave_num\tpeptide_min_len')
    print('\t' + enzymeName + '\t' + str(miscleave_num) + '\t' + str(peptide_min_len))
    print('uniquePeptide num: ' + str(len(uniquePeptides)))
    

    
"""
#cut by LysC first
peptides = []
for protein in ups:
   peptides += cutByPattern(protein, ['K'], '')   

#cut by trypsin following
newPeptides = []
for peptide in peptides:
    newPeptides += cutByPattern(peptide, ['K', 'R'], 'P')
#print(len(newPeptides))
"""


upsFile = "/Users/hao/data/database/ups1-ups2-sequences.fasta"
cutPatterns = {'trypsin': ['K','R'], 'LysC': ['K']}
excludePatterns = {'trypsin' : 'P', 'LysC' : '' }

for enzymeName in ['trypsin', 'LysC']:
    for peptide_min_len in [1, 7]:
        for miscleave_num in [0, 1, 3]:
            statisticOfUniqPeptideBeCut(upsFile, enzymeName, miscleave_num, peptide_min_len)



    



