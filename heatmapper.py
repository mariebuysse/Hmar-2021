# pseudogenization level representation
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import pandas as pd
import sys
# requires table with number of COGs for category for each genome
data=sys.argv[1]
tab = pd.read_csv(data, index_col=0,sep="\t").T
tab_minmax=(tab/tab.max())
print(tab_minmax)
colour = sns.color_palette("rocket_r", as_cmap=True)
sns.heatmap(tab_minmax, robust=True, xticklabels=1,cmap=colour,linewidths=.5)
plt.show()
