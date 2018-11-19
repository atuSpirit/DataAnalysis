"""Transform peaks db result to the de novo format could be accepted by TBnovo"""
import pandas as pd
import re
from generate_dn_for_peaks_ab import fake_conf_score 

def extract_db(psm_file):
	"""Extract peptide and intensity from psm file and assign conf score to 99"""
	db = pd.read_csv(psm_file)
	db_result = db.loc[:, ['Peptide', 'Mass', 'Length', 'ppm', 'm/z', 'RT', 'Fraction', 'Scan']]

	return db_result

def map_fraction_to_enzyme(frac_num):
	switcher = {
		1 : 'PEPSIN', 
		2 : 'PROK',
		3 : 'TRYPSIN',
		4 : 'CHYMOTRYPSIN',
		5 : 'TRYPSIN',
		11 : 'PEPSIN',
		12 : 'PROK',
		13 : 'TRYPSIN' 
	}
	return switcher.get(frac_num, "Invalid enzyme")

def remove_modification(peptide):
	peptide = re.sub('\(\S+\)', '', peptide)
	return peptide


def db_to_tbnovo(db_result):
	"""Transform peaks db result to TBnovo bottom_up_seq"""
	tbnovo_columns = ['file_name', 'enzyme', 'Scan', 'Peptide', 'ppm', 'ALC', 'm/z', 'z', 'RT', 'Mass', 'ppm2', 'N', 'conf', 'Peptide2']

	tbnovo = pd.DataFrame(columns=tbnovo_columns)
	tbnovo['file_name'] = db_result['Fraction'].apply(map_fraction_to_enzyme)
	tbnovo['enzyme'] = tbnovo['file_name']
	tbnovo['Scan'] = db_result['Scan'].apply(lambda x: x.split(':')[1])
	tbnovo['Peptide'] = db_result['Peptide'].apply(remove_modification)
	tbnovo['ppm'] = db_result['ppm'].apply(abs)
	tbnovo['ALC'] = 99	#For all DB result, set the ALC to be 99
	tbnovo['m/z'] = db_result['m/z']
	tbnovo['z'] = db_result['Mass'] / (db_result['m/z'] - 1)
	tbnovo['z'] = tbnovo['z'].apply(round)
	tbnovo['z'] = tbnovo['z'].apply(int)
	tbnovo['RT'] = db_result['RT']
	tbnovo['Mass'] = db_result['Mass']
	tbnovo['ppm2'] = tbnovo['ppm']
	tbnovo['N'] = 'N'
	tbnovo['conf'] = tbnovo['Peptide'].apply(fake_conf_score)
	tbnovo['Peptide2'] = tbnovo['Peptide']

	print(tbnovo.head(10))
	return tbnovo
	
def extract_dn(dn_file):
	dn = pd.read_csv(dn_file)
	dn_result = dn.loc[:, ['Fraction', 'Scan', 'Peptide', 'ALC (%)', 'm/z', 'z', 'RT', 'Mass', 'ppm', 'local confidence (%)']]
	return dn_result

def dn_to_tbnovo(dn_result):
	"""Transform de novo result of peaks to TBnovo bottom_up_seq format"""
	tbnovo_columns = ['file_name', 'enzyme', 'Scan', 'Peptide', 'ppm', 'ALC', 'm/z', 'z', 'RT', 'Mass', 'ppm2', 'N', 'conf', 'Peptide2']

	tbnovo = pd.DataFrame(columns=tbnovo_columns)
	tbnovo['file_name'] = dn_result['Fraction'].apply(map_fraction_to_enzyme)
	tbnovo['enzyme'] = tbnovo['file_name']
	tbnovo['Scan'] = dn_result['Scan'].apply(lambda x: x.split(':')[1])
	tbnovo['Peptide'] = dn_result['Peptide'].apply(remove_modification)
	tbnovo['ppm'] = dn_result['ppm'].apply(abs)
	tbnovo['ALC'] = dn_result['ALC (%)']
	tbnovo['m/z'] = dn_result['m/z']
	tbnovo['z'] = dn_result['z']

	tbnovo['RT'] = dn_result['RT']
	tbnovo['Mass'] = dn_result['Mass']
	tbnovo['ppm2'] = tbnovo['ppm']
	tbnovo['N'] = 'N'
	tbnovo['conf'] = dn_result['local confidence (%)']
	tbnovo['Peptide2'] = tbnovo['Peptide']

	print(tbnovo.head(10))
	print(tbnovo.loc[[0,1], ['ALC', 'm/z', 'z', 'RT']])
	return tbnovo



psm_file = "/Users/hao/Data/topDown/mab.bottom_up/peaks.ab.vs.trim/DB search psm.csv"
output_file = psm_file.strip(".csv") + ".forTBnovo.csv"
db_result = extract_db(psm_file)
db_tbnovo = db_to_tbnovo(db_result)
db_tbnovo.to_csv(output_file, header=False,index=False)

"""
#dn_file = "/Users/hao/Data/topDown/mab.bottom_up/hcd.peaks.db/de novo only peptides.csv"
dn_file = "/Users/hao/Data/topDown/mab.bottom_up/peaks.ab.vs.trim/de novo only peptides.csv"
output_file = dn_file.strip(".csv") + ".forTBnovo.csv"
dn_result = extract_dn(dn_file)
dn_tbnovo = dn_to_tbnovo(dn_result)
dn_tbnovo.to_csv(output_file, header = False, index = False)
"""