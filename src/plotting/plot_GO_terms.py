from argparse import ArgumentParser
import matplotlib
from matplotlib import pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42 # Editable text when exporting pdfs
import numpy as np
import pandas as pd
import scipy as sp
import seaborn as sns

# Create a dataframe
df = pd.read_csv("results/communication/GO_enrichment_filtered.csv") 
# Reorder it following the values:
ordered_df = df.sort_values(by='FDR')
my_range=range(1,len(df.index)+1)
 
# The vertival plot is made using the hline function
# I load the seaborn library only to benefit the nice looking feature
plt.figure(figsize=(4,4))
#plt.hlines(y=my_range, xmin=0, xmax=ordered_df['FDR'], color='black')
plt.plot(-np.log10(ordered_df['FDR']), my_range, "o", c = 'k')
 
# Add titles and axis names
plt.yticks(my_range, ordered_df["GO molecular function complete"].values)
#plt.title("Gene enrichment", loc='left')
plt.xlabel('-log10(FDR)',fontsize = 14)
#plt.ylabel('GO term')
plt.savefig("figures/communication/GO_terms.pdf")