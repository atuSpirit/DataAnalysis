#Plot the ppm distribution histogram given the true peaks file and theoretic peaks of reference sequence
import pandas as pd
import matplotlib.pyplot as plt

true_peaks_file = "D:/Hao/data/for_analysis/ppmDistribution/real_peaks.all.txt"
theoretic_peaks_file = "D:/Hao/data/for_analysis/ppmDistribution/theoretic_peaks.txt"

true_peaks = pd.read_csv(true_peaks_file, sep='\t')
theoretic_peaks = pd.read_csv(theoretic_peaks_file, sep=' ');
tolerance = 1

len1 = len(true_peaks)

len2 = len(theoretic_peaks)
ppm_values = []

j = 0
for i in range(0, len1):
	mass = true_peaks.loc[i, 'mass']
	#print("mass" + str(mass))
	theo_mass = theoretic_peaks.iloc[j, 1]
	while j < len2 and mass > (theo_mass - tolerance):    	
		theo_mass = theoretic_peaks.iloc[j, 1]
		#print(theo_mass)
		if abs(mass - theo_mass) < tolerance:			
			ppm = (mass-theo_mass) * 1e6 / mass
			ppm_values.append(ppm)	
			if (ppm > 10):	
				print("mass: " + str(mass))
				print("matched: " + str(theo_mass))
			break				

		j += 1

print(len(ppm_values))
plt.title("Histogram of ppm distribution of matched peaks.all")
plt.xlabel("ppm")
plt.ylabel("frequency")
plt.hist(ppm_values)


plt.title("mass vs ppm of matched peaks.all")
plt.xlabel("mass")
plt.ylabel("ppm")
mass_list = true_peaks['mass'].tolist()[1:-1]	#The first is 0, the last of true_peaks are precurMass
plt.scatter(mass_list, ppm_values)

plt.show()



