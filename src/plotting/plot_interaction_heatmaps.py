from argparse import ArgumentParser
import matplotlib
from matplotlib import pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42 # Editable text when exporting pdfs
import numpy as np
import pandas as pd
import scipy as sp
import seaborn as sns
from sklearn.cluster import bicluster

#%%
def make_heatmap( interaction_scores, predicted_cell_types, sample_color_legend, ct_color_legend,
                 cell_types = None, pct_cells = 0.2, score_cutoff = 0.5 ):
    
    if cell_types is None:
        cell_types = np.unique(predicted_cell_types)
    
    ct_idx = np.array([ct in cell_types for ct in predicted_cell_types.predicted_cell_type])
    num_cells = np.sum(ct_idx)
    ct_scores = interaction_scores.loc[ct_idx,:]

    # Filter by percentage of cells above threshold
    ct_scores = ct_scores.loc[:,ct_scores.any()] # drop scores that don't occur in cell type
    score_idx = np.sum(ct_scores > score_cutoff, axis = 0 ) > (pct_cells*num_cells)
    ct_scores = ct_scores.loc[:,score_idx]
    
    df = pd.DataFrame(ct_scores).transpose()
    cg = sns.clustermap( df, cmap = "Greys",
                         col_colors= [ ct_color_legend.loc[ct_idx], sample_color_legend.loc[ct_idx]],
                         yticklabels = 1)
    cg.ax_row_dendrogram.set_visible(False)
    cg.ax_col_dendrogram.set_visible(False)
    cg.ax_heatmap.set_xticklabels([])
    #plt.savefig("figures/communication/heatmaps/" + ct + ".pdf")
    
#%%
interaction_scores = pd.read_csv("results/communication/interaction_scores.csv",index_col = 0)
interaction_scores = interaction_scores .loc[:,interaction_scores .any()]

predicted_cell_types = pd.read_csv("results/classification/processed/predicted_cell_types.csv",index_col = 0 )
#ct_colors = predicted_cell_types["predicted_cell_type"].map(color_key)
max_val = np.max(interaction_scores,axis = 0)
min_val = np.min(interaction_scores,axis = 0)   
normalized_scores = (interaction_scores - min_val) / (max_val - min_val)
#%% Color legends
sample_names = np.array([("_").join(x.split("_")[0:3]) for x in predicted_cell_types.index.values] )
samples = np.unique(sample_names)
color_map = plt.get_cmap("Pastel1")
sample_color_key = dict(zip( samples, color_map.colors[0:len(samples)] ) )
sample_color_legend = pd.Series(sample_names,index = predicted_cell_types.index.values).map(sample_color_key)

cell_types = np.unique( predicted_cell_type )
color_map = plt.get_cmap("tab20")
ct_color_key = dict(zip( cell_types, color_map.colors[0:len(cell_types)] ) )
ct_color_key["Unassigned"] = (0.8,0.8,0.8)
ct_color_legend = predicted_cell_types["predicted_cell_type"].map(ct_color_key)

#%% 
make_heatmap( normalized_scores, predicted_cell_types, sample_color_legend, ct_color_legend,
                 cell_types = ["Bcell","Tcell"], pct_cells = 0.05, score_cutoff = 0.33 ) # Lymphocyte
plt.savefig("figures/communication/lymphocyte_interactions.pdf")

make_heatmap( normalized_scores, predicted_cell_types, sample_color_legend, ct_color_legend,
                 cell_types = ["Monocyte","Macrophage"], pct_cells = 0.125, score_cutoff = 0.33 ) # Myeloid
plt.savefig("figures/communication/myeloid_interactions.pdf")

make_heatmap( normalized_scores, predicted_cell_types, sample_color_legend, ct_color_legend,
                 cell_types = ["Colonocyte","Goblet","Stem","Tumor"], pct_cells = 0.03, score_cutoff = 0.3 ) # Epithelial
plt.savefig("figures/communication/epithelial_interactions.pdf")

make_heatmap( normalized_scores, predicted_cell_types, sample_color_legend, ct_color_legend,
                 cell_types = ["Myofibroblast","CAF"], pct_cells = 0.25, score_cutoff = 0.33 ) # stromal
plt.savefig("figures/communication/stromal_interactions.pdf")