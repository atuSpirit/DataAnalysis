import pandas as pd
import matplotlib.pyplot as plt
import pygal 

peaks_dataframe = pd.read_csv("D:/Hao/data/for_analysis/tmp.txt", sep = '\t');

peak_bar_list = []
for index, row in peaks_dataframe.iterrows():
	height = row[1]
	#print(height)
	mass = row[0]
	bar = (height, mass, mass + 0.01)
	peak_bar_list.append(bar)

hist = pygal.Histogram()
hist.title = "HCD first spectrum"

hist.add('peaks', peak_bar_list)


hist.render_in_browser()
