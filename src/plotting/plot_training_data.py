from argparse import ArgumentParser
import matplotlib
from matplotlib import pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42 # Editable text when exporting pdfs
import numpy as np
import pandas as pd
import scipy as sp
import seaborn as sns




training_data = np.array( sp.sparse.load_npz( "results/classification/sparse_training_counts.npz").todense())
training_labels = pd.read_csv( "results/classification/training_labels.csv", index_col = 0) 

umap_coordinates = pd.read_csv("results/UMAP/processed_UMAP_coordinates.csv",index_col = 0 )


df = pd.merge( umap_coordinates, training_labels, how = "inner", left_index = True, right_index = True )

cell_types = np.unique(df.training_label.values )

plt.figure(figsize=(9,9))
for ct in cell_types:
    idx = df.training_label.values  == ct
    plt.scatter( df.UMAP1.values[idx], df.UMAP2.values[idx], s = 2)
plt.legend(cell_types)




plt.figure(figsize=(9,9))
plt.scatter( df.UMAP1.values, df.UMAP2.values, s = 2, 
             c = np.max( training_probabilities, axis = 1 ))
plt.colorbar()
