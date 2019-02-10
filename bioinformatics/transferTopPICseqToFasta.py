#Covert the proteoform column of the result of ToPIC software to fasta file, 
#so that the fasta file can be pairwise aligned with template protein sequence
#One example ToPIC sequence is G.TQTYI((C)[Carbamidomethylation]NVNHKPSNTKVDKRVEPKS(C)[Carbamidomethylation]D)[-4.98935]KT.H
#The output should be 
#>1
#TQTYICNVNHKPSNTKVDKRVEPKSCDKT
import re

def transfer_each_seq(proteo_seq):
	proteo_seq = re.sub(r'^\S?\.', '', proteo_seq);
	proteo_seq = re.sub(r'\.\S?$', '', proteo_seq);
	proteo_seq = re.sub(r'\[.*?\]', '', proteo_seq)
	proteo_seq = re.sub(r'\(', '', proteo_seq)
	proteo_seq = re.sub(r'\)', '', proteo_seq)
	return proteo_seq

file = "d:/Hao/data/for_analysis/topPIC/hcd1.txt"
fasta_seq = ""
with open(file) as fobj:	
	line = fobj.readline()	#header
	line = fobj.readline()
	cnt = 1
	while line:				
		new_seq = transfer_each_seq(line.strip())
		fasta_seq += ">" + str(cnt) + "\n" + new_seq + "\n";
		cnt += 1
		line = fobj.readline()

print(fasta_seq)


