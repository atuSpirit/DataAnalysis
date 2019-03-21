#Given a set of psms and template number k, find k proteins whose psm number are maximum.
#Print the list of top k template id and the number of psm mapped to it.  If the id is empty,
#it's because of psms with no significant protein
import re

def read_psm_file(psm_file):
	"""For each psm in psm_file, extract its corresponding protein id list"""
	with open(psm_file) as fobj:
		title_line = fobj.readline()	#The title line
		(psm_index, pid_index) = parseTitleLine(title_line)

		line = fobj.readline()
		psm_pids = {}
		while line:			
			(psm, pid_list) = get_pid_list(line, psm_index, pid_index)
			#If psm has pid, then store it 
			if len(pid_list) > 0:
				psm_pids[psm] = pid_list
			line = fobj.readline()

	return psm_pids

def get_pid_psms_map(psm_pids):
	pid_psms = {}

	for psm, pids in psm_pids.items():
		for pid in pids:
			if pid in pid_psms.keys():
				pid_psms[pid].append(psm)
			else:
				pid_psms[pid] = [psm]

	return pid_psms


def parseTitleLine(title_line):
	"""Parse the title line, return the index of scan and the index of Accession."""
	title_line.rstrip()
	field_names = title_line.split(",")
	size = len(field_names)

	psm_index = -1
	pid_index = -1
	for i in range(size):
		if field_names[i] == "Scan":
			psm_index = i
		if field_names[i] == "Accession":
			pid_index = i

	return (psm_index, pid_index)

def get_pid_list(line, psm_index, pid_index):
	"""Extract psm scan and its corresponding protein id list from current line""" 
	line.rstrip()
	fields = line.split(',')
	psm = fields[psm_index]	
	pid_string = fields[pid_index]
	pid_list = pid_string.split(':')

	return (psm, pid_list)

def get_pid_with_maximum_psms(pid_psms):
	"""Find the protein id who has the maximum psm number in the pid_psms dictionary"""
	max_psm_num = 0
	pid_with_max_psms = 0
	for pid,psms in pid_psms.items():		
		if len(psms) > max_psm_num:			
			max_psm_num = len(psms)
			pid_with_max_psms = pid
	return pid_with_max_psms

def remove_pid_with_maximum_psms(pid_psms, pid_with_max_psms):
	"""Remove the pid entry with maximum psms from pid_psms.
	   Remove the psms shared with those of pid """	
	psm_list = pid_psms[pid_with_max_psms]
	del pid_psms[pid_with_max_psms]

	psm_set = set(psm_list)

	for pid, psms in pid_psms.items():
		#print("before: " + str(len(psms)))
		new_psms = []
		for psm in psms:			
			if psm not in psm_set:
				new_psms.append(psm)
				
		#print("after: " + str(len(new_psms)))
		pid_psms[pid] = new_psms
		
	return pid_psms

def pick_top_k_pids(pid_psms, k):
	pid_list = []
	"""Extract k proteins with maximum non overlapped psms"""
	for i in range(k):
		pid_with_max_psms = get_pid_with_maximum_psms(pid_psms)
		print(pid_with_max_psms + " " + str(len(pid_psms[pid_with_max_psms])))
		if pid_with_max_psms != "":
			pid_list.append(pid_with_max_psms)

		pid_psms = remove_pid_with_maximum_psms(pid_psms, pid_with_max_psms)

	return pid_list

def export_template_seqs(template_pid_list, database_file, template_file):
	template_pid_map = {pid : 1 for pid in template_pid_list}	
	template_seq = ""
	export = False
	with open(database_file) as fobj:
		for line in fobj:
			if line.startswith('>'):
				export = False
				match_obj = re.match(r'^>ab(\S+)\n', line, re.M|re.I)
				if match_obj:
					template_pid = match_obj.group(1)					
					if template_pid in template_pid_map.keys():
						export = True
			if export == True:
				template_seq += line

	with open(template_file, 'w') as fobj:
		fobj.write(template_seq)



def test():
	print(len(psm_pids))
	print(len(psm_pids["F8:4669"]))
	print(len(pid_psms))
	print(len(pid_psms["|Q6GI34|SSPA_STAAR"]))
	print(pid_psms["|Q6GI34|SSPA_STAAR"])
	pid = get_pid_with_maximum_psms(pid_psms)
	print(pid)
	print(len(pid_psms[pid]))

psm_file = "C:/Users/hlin/PeaksExports/Nuno2016_LC_SPIDER_12/DB search psm.csv"
psm_file = "C:/Users/hlin/PeaksExports/PolyClonal_ab19001_SPIDER_7/DB search psm.csv"
psm_file = "C:/Users/hlin/PeaksExports/Hieu_WlgG1_heavy_SPIDER_7/DB search psm.csv"
psm_file = "D:/Hao/result/Waters_mAB_SPIDER_9/DB search psm.csv"
psm_file = "D:/Hao/result/ab19001.5enzymes_SPIDER_66/DB search psm.csv"
psm_pids = read_psm_file(psm_file)
pid_psms = get_pid_psms_map(psm_pids)
top_k = 8
template_pid_list = pick_top_k_pids(pid_psms, top_k + 1)	# plus 1 due to empty pid

database_file = "D:/Hao/data/database/antibodies_with_contaminants.fasta"
database_file = "D:/Hao/result/ab19001.5enzymes_SPIDER_66/candidate_template3.fasta"
template_file = "D:/Hao/data/Polyclonal/ab19001Result/ab19001.template_top" + str(top_k) + ".fasta"
template_file = "D:/Hao/data/Hieu_data/WlgG1_heavy.template_top" + str(top_k) + ".fasta";
template_file = "D:/Hao/result/Waters_mAB_SPIDER_9/Waters_mAB.template_top" + str(top_k) + ".fasta"
template_file = "D:/Hao/result/ab19001.5enzymes_SPIDER_66/ab19001.5enzymes.template_top" + str(top_k) + ".fasta"
export_template_seqs(template_pid_list, database_file, template_file)



