"""Exclude the proteins which is similar to known antibody. 
   Extract the protein ids which has high similarity reported 
   in blastp result of a given antibody sequence vs antibody database. 
   Exclude these ids to generate a new antibody database
""" 
import re

def extract_protein_id_from_blastp_result(blastp_output):
	seq_ids = {}
	with open(blastp_output) as fobj:
		for line in fobj:
			match_obj = re.match(r'^  sp\|(\d+)\|', line, re.M|re.U);
			if match_obj:
				seq_id = match_obj.group(1)
				seq_ids[seq_id] = 1

	return seq_ids

def exclude_seq_from_antibody_database(database, seq_ids):
	new_database_seqs = [] 
	with open(database) as dobj:
		lines = dobj.readlines()

	seq_num = len(lines)
	i = 0
	skip = False
	while i < seq_num:
		line = lines[i];
		match_obj = re.match(r'^>sp\|(\d+)\|', line, re.M|re.I)
		if match_obj:
			seq_id = match_obj.group(1)
			if seq_id in seq_ids.keys():
				skip = True
			else:
				skip = False
		if skip == False:
			new_database_seqs.append(line)
	
		i += 1
	return new_database_seqs

def seqs_to_file(sequences, output_file):
	with open(output_file, 'w') as out_obj:		
		for line in sequences:
			out_obj.write(line)


blastp_file = "D:/Hao/data/database/antibody_trim/true_seq.vs.mAB_database.blastp.out"
database_file = "D:/Hao/data/database/antibody_trim/mAB_database.ann"
output_file = database_file + '.trim'
seq_ids = extract_protein_id_from_blastp_result(blastp_file)
new_database_seqs = exclude_seq_from_antibody_database(database_file, seq_ids)
seqs_to_file(new_database_seqs, output_file)




