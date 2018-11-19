"""Extract ios with given mass from the mgf file exported by ms-deconv software. """
import re
import matplotlib.pyplot as plt

def load_mgf_file(mgf_file):
	with open(mgf_file) as fobj:
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

			i += 1
			title = ion_lines[i]
			#print(title)
			ion += title
			title_match = re.match(r'^TITLE=Scan_(\d+)', title, re.M|re.I)
			if title_match:
				scan = title_match[1]
				#print(scan)
				i += 1
				ms_level = ion_lines[i]
				ion += ms_level

				i += 1
				pep_mass_line = ion_lines[i]
				pep_mass_match = re.match(r'PEPMASS=(\d+)\.(\d+)', pep_mass_line, re.M|re.I)
				if pep_mass_match:
					mass = pep_mass_match[1] + '.' + pep_mass_match[2]
					mass = float(mass)

					if abs(mass - pep_mass) < mass_diff_tolerance:
						mass_equal = True
						ion += pep_mass_line
					else:
						ion = ""

					i += 1
					charge_line = ion_lines[i]
					ion += charge_line

					i += 1
					peak_num = 0
					while i < n:
						line = ion_lines[i]
						if mass_equal:
							ion += line
							peak_num += 1
						end_ion = re.match(r'^END', line, re.M|re.I)
						if end_ion:
							if mass_equal:
								ions_of_given_mass[scan] = ion
								peak_num_of_ions[scan] = peak_num
							break
						else:
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
	
def export_ions(ions_of_given_mass, peak_num_of_ions, output_file):
	"""Export the ions according to the number of peaks of ion decendingly."""
	sorted_by_peak_num = sorted(peak_num_of_ions.items(), key=lambda kv:kv[1], reverse=True)

	with open(output_file, 'w') as fobj:
		for t in sorted_by_peak_num:
			scan = t[0]
			fobj.write(ions_of_given_mass[scan])
			fobj.write("\n")

		

mgf_file = "/Users/hao/data/tmp/FAB_MAB_TCEP_HCD_110626213244_msdeconv.mgf"
pep_mass = 23556
output_file = mgf_file.split('.')[0] + ".pep_mass" + str(pep_mass) + ".mgf"
mass_diff_tolerance = 1
ion_lines = load_mgf_file(mgf_file)
print(len(ion_lines))
(ions_of_given_mass, peak_num_of_ions) = filter_ios_by_mass(ion_lines, pep_mass, mass_diff_tolerance)

#peaks_list_num_distribution(peak_num_of_ions)
export_ions(ions_of_given_mass, peak_num_of_ions, output_file)





					





