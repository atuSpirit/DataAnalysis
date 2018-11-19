import pandas as pd 
import re

pd.set_option('display.max_rows',100)
pd.set_option('display.max_columns',10)
pd.set_option('display.max_colwidth',20)

def extract_dn(denovo_file):
	"""Extract peptide, conf score and intensity from de novo only file"""
	dn = pd.read_csv(denovo_file)
	dn_result = dn.loc[:, ['Scan', 'Peptide', 'local confidence (%)', 'Intensity']]
	dn_result.columns = ['scan', 'peptide', 'conf' , 'area']

	#Keep the first scan number
	dn_result['scan'] = dn_result['scan'].apply(lambda x : x.split(' ')[0])
	dn_result = dn_result.fillna(100)

	return dn_result

def fake_conf_score(db_seq):
	"""Set the conf score for each AA in peptide to 99"""
	conf = ""
	for a in db_seq:
		if a.isupper():
			conf += '99 '
		elif a.islower():
			print(db_seq)

	return conf.strip()

def modify_AA(spider_peptide):
	modified_peptide = []
	length = len(spider_peptide)
	i = 0
	flag = 0
	while i < length:
		a = spider_peptide[i]
		#If a sequence is in form of A(sub T), A is replaced by T
		if a == '(':
			tmp = spider_peptide[(i + 1):(i+4)]
			if  tmp == 'sub':
				a = spider_peptide[i + 5]
				modified_peptide.pop()
				i += 6   #skip the ')'
			elif tmp == 'ins':
				a = ''
				i += 4
			elif tmp == 'del':
				a = ''
				i += 6
		modified_peptide.append(a)
		i += 1

	return ''.join(modified_peptide)


def extract_db(psm_file, is_spider):
	"""Extract peptide and intensity from psm file and assign conf score to 99"""
	db = pd.read_csv(psm_file)
	db_result = db.loc[:, ['Scan', 'Peptide', 'Intensity']]
	db_result.columns = ['scan', 'peptide', 'area']

	if (is_spider):
		db_result['peptide'] = db_result['peptide'].apply(modify_AA)

	db_result['conf'] = db_result['peptide'].apply(fake_conf_score)
	db_result = db_result.fillna(100)
	
	db_result = db_result[['scan', 'peptide', 'conf', 'area']]
	
	return db_result



"""
peaks_result_dir = '/Users/hao/data/Nuno.spider/nuno.vs.HCLC/Nuno2016.HC_SPIDER_25'
peaks_result_dir = '/Users/hao/data/topDown/mab.bottom_up'
is_spider = True
dn_result = extract_dn(peaks_result_dir + '/de novo only peptides.csv')
dn_result.to_csv(peaks_result_dir + '/peptide_0.csv', index=False)
db_result = extract_db(peaks_result_dir + '/DB search psm.csv', is_spider)
db_result.to_csv(peaks_result_dir + '/peptide_db.csv', index=False)

peptides = pd.concat([db_result, dn_result])
peptides.to_csv(peaks_result_dir + '/peptide_combine.csv', index=False)
"""
