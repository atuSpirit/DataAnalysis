import pandas as pd
import matplotlib.pyplot as plt

file = "D:/Hao/data/mab/top-down/TopPIC.1.2.2.vs.antibody_with_true/mab.topdown.lengthDistri.csv"
data = pd.read_csv(file)

length = data['length'];
"""
binwidth = 10

plt.hist(length, bins = range(min(length), max(length) + binwidth, binwidth))
plt.title("Length distribution of db sequence of top down spectra")
plt.xlabel("length")
plt.ylabel("freq")
plt.show()
"""

ratio = data['matched_peaks'] / data['length'] 
plt.scatter(length, ratio)
plt.title("Ratio distribution of number of matched peaks / length")
plt.xlabel("Length")
plt.ylabel("Ratio of matched peaks number")
plt.show()
