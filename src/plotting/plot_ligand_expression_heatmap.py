import pandas as pd
import numpy as np
import scipy as sp
import matplotlib
from matplotlib import pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42 # Editable text when exporting pdfs
import seaborn as sns

from src import cell_cell_communication as ccc

scRNA_expression = np.array( sp.sparse.load_npz( "results/classification/processed/processed_counts.npz" ).todense() )
gene_symbols = np.squeeze(np.array(pd.read_csv("data/mouse_gene_symbols.csv",header=None)))
predicted_cell_type = np.squeeze(np.array(pd.read_csv("results/classification/processed/predicted_cell_types.csv",index_col = 0)))
sample_labels = np.squeeze(np.array(pd.read_csv("results/classification/processed/sample_labels.csv",index_col=0)))
interaction_scores = pd.read_csv("results/communication/interaction_scores.csv",index_col = 0)


sorted_tumor_diff_interactions = pd.read_csv("results/communication/tumor_diff_interactions.csv")


ligand_names = [ x.split("_")[0] for x in sorted_tumor_diff_interactions.loc[:,"name"]]
receptor_names = [ x.split("_")[1] for x in sorted_tumor_diff_interactions.loc[:,"name"]]
unique_ligand = np.unique(np.array(ligand_names))


model = np.array([x.split("_")[0] for x in sample_labels])
tissue = np.array([x.split("_")[1] for x in sample_labels])
idx = np.logical_and( model == "AOM", tissue == "tumor")


cell_type_ligand_expresssion = ccc.cell_type_specific_expression(scRNA_expression[idx,:], gene_symbols, unique_ligand, predicted_cell_type[idx])
average_ligand_expression = np.array([np.mean(x,axis = 1) for x in cell_type_ligand_expresssion])
df = pd.DataFrame( sp.stats.zscore(average_ligand_expression), columns = unique_ligand, index = np.unique(predicted_cell_type))

plt.figure(figsize=(10,8))
hMap = sns.heatmap( df, center = 0, cmap = "coolwarm",square=True,
                   xticklabels=True,yticklabels = True)
#hMap.set_xticklabels(labels = unique_ligand,rotation=45)
plt.savefig("figures/communication/tumor_diff_interactions_ligand_expression.pdf")



unique_ligands = np.unique([col.split("_")[0] for col in interaction_scores.columns])
unique_receptors = np.unique([col.split("_")[1] for col in interaction_scores.columns])