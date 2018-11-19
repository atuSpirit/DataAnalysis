import re

"""Add SEQ=EMPTY to each spectrum in MGF file to be used as input of DeepNovo"""
mgf_file_path = "/Users/hao/data/topDown/mab.top_down/1.mgf"
modified_file = mgf_file_path + "for_DeepNovo.mgf"
with open(mgf_file_path) as fb:
	lines = fb.readlines()

seq_string = "SEQ=EMPTY\n"
with open(modified_file, 'w') as fout:
	for line in lines:
		match = re.search("RTINSECONDS", line)
		if match:
			print(line)
			fout.write(line)
			fout.write(seq_string)
		else:
			fout.write(line)



