import pandas as pd 
import matplotlib.pyplot as plt

pd.set_option('display.max_rows',100)
pd.set_option('display.max_columns',10)
pd.set_option('display.max_colwidth',20)

def read_protein_peptide_csv(protein_peptide_csv):
	protein_peptides = pd.read_csv(protein_peptide_csv)
	
	protein_peptides = protein_peptides.loc[:, ['Protein ID', 'Peptide', 'Unique', 'Intensity Sample 2']]
	unique_protein_peptides = protein_peptides.loc[protein_peptides['Unique'] == 'Y', :].dropna()
	
	return unique_protein_peptides

def intensity_scatter_by_protein_ID(unique_protein_peptides):
	plt.scatter(unique_protein_peptides['Protein ID'], unique_protein_peptides['Intensity Sample 2'], s=2)
	plt.show()

protein_peptide_csv = '/Users/hao/data/Nuno.spider/nuno.vs.HCLC/Nuno2016.LC_SPIDER_13/protein-peptides.csv'
unique_protein_peptides = read_protein_peptide_csv(protein_peptide_csv)
print(len(unique_protein_peptides))
intensity_scatter_by_protein_ID(unique_protein_peptides)