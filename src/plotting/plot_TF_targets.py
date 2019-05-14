from argparse import ArgumentParser
import matplotlib
from matplotlib import pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42 # Editable text when exporting pdfs
import numpy as np
import pandas as pd
import scipy as sp
import seaborn as sns

def normalize_arcsinh(count_matrix,cofactor = 1000):
    total_counts = np.squeeze(np.sum(count_matrix,axis = 1))
    normalized_counts = np.arcsinh(count_matrix / total_counts[:,np.newaxis] * cofactor)
    return(normalized_counts)


TF = "Stat3"
df = pd.read_csv("results/Dorothea/TF_target_expression/"+TF+"_target_expression.csv",index_col = 0) 
df = np.abs(df)
df = df.loc[df.any(axis=1),:]
sign = np.sign(df)
#df = df.loc[:,df.any()]


normalized_counts = normalize_arcsinh( np.array( df ) )
normalized_counts = normalized_counts*sign

# Get sample colors
sample_names = [ ("_").join(x.split("_")[0]) for x in df.index.values]
samples = np.unique(sample_names)
color_map = plt.get_cmap("Pastel1")
sample_color_key = dict(zip( samples, color_map.colors[0:len(samples)] ) )
sample_color_legend = pd.Series(sample_names,index = normalized_counts.index.values).map(sample_color_key)



cg = sns.clustermap( normalized_counts.transpose(),
               center = 0, cmap = "coolwarm",
               col_colors= sample_color_legend,
               yticklabels = 1)
cg.ax_row_dendrogram.set_visible(False)
cg.ax_col_dendrogram.set_visible(False)
cg.ax_heatmap.set_xticklabels([])
plt.savefig("figures/Dorothea/"+TF+"_target_expression.pdf")