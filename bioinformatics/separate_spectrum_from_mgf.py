#Extract each spectrum in a mgf file and form a file with only one spectrum

spectrum_dir = "d:/Hao/data/herceptin.top_down/mgf_dir"
mgf_file = "d:/Hao/data/herceptin.top_down/herceptin_filter5940-7830s_ms2.addCharge.msalign.mgf"

with open(mgf_file) as fobj:
	for line in fobj:
		if line.startswith("BEGIN IONS"): 
			spectrum = line
		else:
			spectrum += line
			if line.startswith("TITLE"):
				#TITLE=Scan_2713				
				fields = line.split("=")
				scan = fields[1].strip('\n')				
			if line.startswith("END IONS"):
				out_file = spectrum_dir + "/" + scan + ".mgf";
				with open(out_file, "w") as out:
					out.write(spectrum)



		

