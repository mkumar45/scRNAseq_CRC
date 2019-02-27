from argparse import ArgumentParser
import matplotlib
from matplotlib import pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42 # Editable text when exporting pdfs
import numpy as np
import pandas as pd
import scipy as sp
import seaborn as sns

scRNA_expression = np.array(sp.sparse.load_npz("data/combined/processed_counts.npz").todense())
predicted_cell_type = pd.read_csv("results/classification/predicted_cell_types.csv",index_col=0)
prediction_probabilities = pd.read_csv("results/classification/prediction_probabilities.csv",index_col=0)

high_confidence_idx =  np.max( prediction_probabilities, axis = 1 ) > 0.5


confident_counts = scRNA_expression[ high_confidence_idx, :]
confident_ct_labels = predicted_cell_type[high_confidence_idx]
group_label = [ "_".join(x.split("_")[0:-1]) for x in confident_ct_labels.index.values]

sp.sparse.save_npz( "results/classification/processed/processed_counts.npz", sp.sparse.coo_matrix( confident_counts) )
confident_ct_labels.to_csv("results/classification/processed/predicted_cell_types.csv")
pd.DataFrame( group_label, index = confident_ct_labels.index.values).to_csv(
              "results/classification/processed/sample_labels.csv")