import pandas as pd
import matplotlib.pyplot as plt
import pygal 

def load_peaks(peak_file):
	peaks = pd.read_csv(peak_file, sep = '\t')
	print(peaks.head(10))
	return peaks

def draw_peaks(real_peak_list, noise_peak_list, fig_title, fig_file):
	hist = pygal.Histogram()
	hist.title = fig_title
	print(fig_title)
	
	hist.add('noise peaks', noise_peak_list)
	
	hist.add('real_peaks', real_peak_list)
	hist.render_to_file(fig_file)
	#hist.render_in_browser()

def transform_peaks_to_bar_list(peaks_dataframe, height_index):
	"""Transfer peaks in df into a list of (height, mass, mass + 0.01)
	   height_index indicate which column to draw, the peak num(index 1) 
	   or intensity(index = 4)"""
	peak_bar_list = []
	for index, row in peaks_dataframe.iterrows():
		height = row[height_index]
		#print(height)
		mass = row[2]
		bar = (height, mass, mass + 0.01)
		peak_bar_list.append(bar)

	return peak_bar_list



dir = "D:/Hao/data/for_analysis/result"
real_peaks_file = dir + "/real_peaks.valuablePeak.twoAA.all.txt"
noise_peaks_file = dir + "/noise_peaks.valuablePeak.twoAA.all.txt"

real_peaks = load_peaks(real_peaks_file)
noise_peaks = load_peaks(noise_peaks_file)

#Draw intensities of real peaks vs noise peaks
height_index = 4
fig_title = "Intensity of Real peaks vs noise peaks of filtered peaks using one or two AA tag"
fig_file = dir + "/fig/intensity.real.vs.noise_peaks.filteredByTwoAAtag.svg"


real_peak_bar_list = transform_peaks_to_bar_list(real_peaks, height_index)
noise_peak_bar_list = transform_peaks_to_bar_list(noise_peaks, height_index)
draw_peaks(real_peak_bar_list, noise_peak_bar_list, fig_title, fig_file)


#Draw supporting number of peaks
height_index = 1
fig_title = "Supporting number of Real peaks vs noise peaks using one or two AA tag"
fig_file = dir + "/fig/supporting_num.1.real.vs.noise_peaks..filteredByTwoAAtag.svg"

real_peak_bar_list = transform_peaks_to_bar_list(real_peaks, height_index)
noise_peak_bar_list = transform_peaks_to_bar_list(noise_peaks, height_index)
draw_peaks(real_peak_bar_list, noise_peak_bar_list, fig_title, fig_file)



