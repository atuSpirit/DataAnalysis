import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
intensity_file = "D:/Hao/result/ab19001.5enzymes_SPIDER_46/intensity.csv"

intensity = pd.read_csv(intensity_file, names=['pos', 'AA', 'inten'])
#plt.scatter(intensity['pos'], intensity['inten'])
#plt.show()
sns.lmplot(x = "pos", y = "inten", data = intensity, fit_reg = False, hue = "AA", legend=False)
plt.legend(loc='lower right')
plt.show()

