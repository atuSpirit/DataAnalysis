"""Extract ios with given mass from the msalign file exported by ms-deconv software. """
import re
import matplotlib.pyplot as plt
from os import listdir 

def load_msalign_file(msalign_file):
	with open(msalign_file) as fobj:
		lines = fobj.readlines()
	return lines


def filter_ios_by_mass(ion_lines, pep_mass, mass_diff_tolerance):
	ions_of_given_mass = {}
	peak_num_of_ions = {}
	n = len(ion_lines)

	i = 0
	mass_equal = False
	while (i < n):
		line = ion_lines[i]
		match_obj = re.match(r'^BEGIN', line, re.M|re.I)
		if match_obj:
			mass_equal = False
			ion = line

			#ID line
			i += 1
			id = ion_lines[i]
			#print(title)
			ion += id

			#SCANS line
			i += 1
			scans = ion_lines[i]
			scans_match = re.match(r'^SCANS=(\d+)', scans, re.M|re.I)
			if scans_match:
				scan = scans_match[1]
				ion += scans
				#print(scan)
				i += 1
				activation = ion_lines[i]
				ion += activation

				i += 1
				mz = ion_lines[i]
				ion += mz

				i += 1
				z = ion_lines[i]
				ion += z

				i += 1
				precursor_mass_line = ion_lines[i]
				#precursor_mass_match = re.match(r'PRECURSOR_MASS=(\d+)\.(\d+)', precursor_mass_line, re.M|re.I)
				precursor_mass_match = re.match(r'PRECURSOR_MASS=(\S+)', precursor_mass_line, re.M|re.I)
				if precursor_mass_match:
					#mass = precursor_mass_match[1] + '.' + precursor_mass_match[2]
					mass = precursor_mass_match[1]
					mass = float(mass)

					if abs(mass - pep_mass) < mass_diff_tolerance:
						mass_equal = True
						ion += precursor_mass_line
					else:
						ion = ""

					#Load peak list of each ion
					peak_list = {}
					i += 1
					peak_num = 0
					while i < n:
						line = ion_lines[i]

						end_ion = re.match(r'^END', line, re.M|re.I)
						if end_ion:					
							if mass_equal:
								if bool(peak_list)== False:
									break
								#If the peak list is not empty, sort them according to mass ascendingly
								sorted_peak_list = sorted(peak_list.items(), key=lambda kv:kv[0])
								sorted_peak_list = list(map(lambda x : x[1], sorted_peak_list))
								ion += ''.join(sorted_peak_list)
								ion += line
								ions_of_given_mass[scan] = ion
								peak_num_of_ions[scan] = peak_num
							break
						else:
							if mass_equal:
								#ion += line
								peak_num += 1
								peak_line_match = re.match(r'(\d+)\.(\d+)\S', line, re.M|re.I)
								if peak_line_match:
									peak_mass = peak_line_match[1] + '.' + peak_line_match[2]
									peak_mass = float(peak_mass)
									peak_list[peak_mass] = line
							i += 1
		i += 1

	print("Number of ions with mass equals to " + str(pep_mass) + " is " + str(len(ions_of_given_mass)))

	return (ions_of_given_mass, peak_num_of_ions)

def peaks_list_num_distribution(peak_num_of_ions):
	"""Draw the distribution of the peak number of each ion"""
	min = 10000
	max = 0
	peak_num_list = peak_num_of_ions.values()

	plt.title("The distribution of the peak number of each ion")
	plt.hist(peak_num_list)
	plt.show()
	
def export_ions(ions_of_given_mass, peak_num_of_ions, output_file, top_num):
	"""Export the ions according to the number of peaks of ion decendingly."""
	sorted_by_peak_num = sorted(peak_num_of_ions.items(), key=lambda kv:kv[1], reverse=True)

	i = 0
	with open(output_file, 'w') as fobj:
		for t in sorted_by_peak_num:
			scan = t[0]
			fobj.write(ions_of_given_mass[scan])
			fobj.write("\n")
			i += 1
			if i == top_num:
				break


def extract_msalign_file(pep_mass, msalign_file, mass_diff_tolerance, top_num):
	output_file = ".".join(msalign_file.split('.')[0:-1]) + ".pep_mass" + str(pep_mass) + ".top" + str(top_num) + ".msalign"
	
	ion_lines = load_msalign_file(msalign_file)
	print(len(ion_lines))
	(ions_of_given_mass, peak_num_of_ions) = filter_ios_by_mass(ion_lines, pep_mass, mass_diff_tolerance)

	#peaks_list_num_distribution(peak_num_of_ions)
	print("out:")
	print(output_file)
	export_ions(ions_of_given_mass, peak_num_of_ions, output_file, top_num)



msalign_dir = "D:/Hao/data/for_analysis/TopPICmsalign/valuablePeaks_twoAA"
dirs = listdir(msalign_dir)

pep_mass = 23556.8
mass_diff_tolerance = 0.5
top_num = 100

for file in dirs:
	msalign_file = msalign_dir + "/" + file
	print(msalign_file)
	extract_msalign_file(pep_mass, msalign_file, mass_diff_tolerance, top_num)

#msalign_file = "D:/Hao/data/for_analysis/TopPICmsalign/FAB_MAB_TCEP_HCD_110626213244_ms2_valublePeak.msalign"
#extract_msalign_file(pep_mass, msalign_file, mass_diff_tolerance, top_num)





					





