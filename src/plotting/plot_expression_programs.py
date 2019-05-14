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
    
    
def plot_DPC( umap_coordinates, cluster_id ):
    fig = plt.figure(figsize=(4,4))
    plt.scatter(umap_coordinates.UMAP1, umap_coordinates.UMAP2, s = 2, c = cluster_id.Cluster, cmap="tab20")
    plt.xlabel("UMAP1",fontsize=16)
    plt.ylabel("UMAP2",fontsize=16)
    plt.xticks([])
    plt.yticks([])
    # Overlay cluster label at center of each cluster
    for n in np.unique(cluster_id):
        member_idx = cluster_id.Cluster == n
        x_center = np.mean( umap_coordinates.UMAP1[member_idx] )
        y_center = np.mean( umap_coordinates.UMAP2[member_idx])            
        plt.text(x_center,y_center, str(n), fontsize = 18 )
    plt.title("DPC Clustering")
    
def umap_gene_expression(count_matrix,gene_symbols,gene,umap_coordiates):
    plt.figure(figsize=(3,3))
    plt.scatter( umap_coordiates.UMAP1, umap_coordiates.UMAP2, s = 2, 
               c = np.squeeze( count_matrix[:, np.squeeze(gene_symbols == gene)] ),
               cmap = "inferno",
               norm=matplotlib.colors.LogNorm() )
    plt.xlabel("UMAP1",fontsize=16)
    plt.ylabel("UMAP2",fontsize=16)
    plt.xticks([])
    plt.yticks([])
    plt.colorbar()
    #cbar = plt.colorbar(ticks = [0.001,5,10,14.91])
    #cbar.ax.set_yticklabels(['0','5','10','15'],fontsize=14)  # vertically oriented colorbar
    plt.title(gene)
    
def umap_signature_expression(count_matrix, gene_symbols, signature_genes, umap_coordinates):
    
    signature_idx = [sym in signature_genes for sym in gene_symbols]
    signature_expression = count_matrix[:,signature_idx]
    max_val = np.max(signature_expression,axis = 0)
    min_val = np.min(signature_expression,axis = 0)
    
    normalized_signature = (signature_expression - min_val) / (max_val - min_val)
    signature_score = np.sum(normalized_signature,axis = 1) / normalized_signature.shape[1]
    plt.figure(figsize=(3,3))
    plt.scatter( umap_coordinates.UMAP1, umap_coordinates.UMAP2, s = 1, 
                c = signature_score, cmap = "inferno" )
    plt.xlabel("UMAP1",fontsize=16)
    plt.ylabel("UMAP2",fontsize=16)
    plt.xticks([])
    plt.yticks([])
    cbar = plt.colorbar(ticks = [0,np.max(signature_score)])
    cbar.ax.set_yticklabels(['0',"{:.2f}".format(np.max(signature_score))],fontsize=14)
    return(signature_score)
#%% Import
UMAP_coordinates = pd.read_csv( "results/UMAP/processed_UMAP_coordinates.csv", index_col = 0)
count_matrix = np.array(sp.sparse.load_npz("data/combined/processed_counts.npz").todense())
sequencing_metrics = pd.read_csv( "data/combined/processed_metrics.csv",index_col = 0)
gene_symbols = np.squeeze(np.array(pd.read_csv("data/mouse_gene_symbols.csv",header=None)))
DPC = pd.read_csv("results/DPC/processed_DPC_coordinates.csv",index_col = 0)
normalized_counts = normalize_arcsinh( count_matrix )
expression_programs = pd.read_csv("data/gene_signatures.csv")

#%% Immune cells 
umap_gene_expression(count_matrix,gene_symbols,"Ptprc",UMAP_coordinates)
plt.savefig("figures/expression/CD45.pdf")

umap_gene_expression(count_matrix,gene_symbols,"Epcam",UMAP_coordinates)
plt.savefig("figures/expression/Epcam.pdf")

umap_gene_expression(count_matrix,gene_symbols,"Thy1",UMAP_coordinates)
plt.savefig("figures/expression/Thy1.pdf")

#%% Plot literature expression programs

for program in expression_programs.columns:
    signature_genes = [x.title().strip() for x in expression_programs[program].dropna()]
    umap_signature_expression(normalized_counts, gene_symbols, signature_genes, UMAP_coordinates)
    plt.title(program)
    plt.savefig("figures/expression/" + program + ".pdf")

#%% subset epithelial cells

epithelial_idx = [c in [14,6,13,5,4,8,3,5] for c in DPC.Cluster ]
epithelial_UMAP = UMAP_coordinates.loc[epithelial_idx,:]
epithelial_counts = count_matrix[epithelial_idx,:]
epithelial_metrics = sequencing_metrics.loc[epithelial_idx,:]
epithelial_DPC = DPC.loc[epithelial_idx,:]
normalized_epithelial = normalize_arcsinh( epithelial_counts )

#%% Plot literature expression programs only using epithelial 
expression_programs = pd.read_csv("data/gene_signatures.csv")

for program in expression_programs.columns:
    signature_genes = [x.title().strip() for x in expression_programs[program].dropna()] #title case "converts" human to mouse by changing capitalization. Should fix
    score = umap_signature_expression(normalized_epithelial, gene_symbols, signature_genes, epithelial_UMAP)
    plt.title(program)
    plt.savefig("figures/expression/" + program + "_epithelial.pdf")
