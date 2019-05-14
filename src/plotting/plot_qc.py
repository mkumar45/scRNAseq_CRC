import os
import pandas as pd
import scipy as sp
from scipy import sparse
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42 # Editable text when exporting pdfs
import os


if not os.path.exists("figures/qc/"):
    os.makedirs("figures/qc/")



UMAP_coordinates = pd.read_csv( "results/UMAP/processed_UMAP_coordinates.csv", index_col = 0)
count_matrix = np.array(sp.sparse.load_npz("data/combined/processed_counts.npz").todense())
sequencing_metrics = pd.read_csv( "data/combined/processed_metrics.csv",index_col = 0)
gene_symbols = np.squeeze(np.array(pd.read_csv("data/mouse_gene_symbols.csv",header=None)))
num_genes_detected = np.squeeze(np.array(np.sum(count_matrix>0,axis=1)))
total_counts = np.squeeze(np.sum(count_matrix,axis = 1))
mt_mouse_idx = pd.Series( gene_symbols ).str.contains("mt-").values # index of mitochondrial genes
mt_counts = np.sum( count_matrix[:,mt_mouse_idx],axis=1 )
mt_pct = mt_counts/total_counts*100


plt.figure(figsize=(4,4))
plt.scatter( UMAP_coordinates.UMAP1.values, UMAP_coordinates.UMAP2.values, s = 2, 
            c = num_genes_detected,
            norm=matplotlib.colors.LogNorm(),cmap = "inferno")
plt.xticks([])
plt.xlabel("UMAP1",fontsize=16)
plt.yticks([])
plt.ylabel("UMAP2",fontsize=16)
cbar = plt.colorbar(ticks = [501,1000,2000,5000,8000])
cbar.ax.set_yticklabels(['500','1000','2000','5000','8000'],fontsize=14)  # vertically oriented colorbar
plt.savefig("figures/qc/num_genes_detected.pdf")


plt.figure(figsize=(4,4))
plt.scatter( UMAP_coordinates.UMAP1.values, UMAP_coordinates.UMAP2.values, s = 2, 
            c = mt_pct, cmap = "inferno")
plt.xticks([])
plt.xlabel("UMAP1",fontsize=16)
plt.yticks([])
plt.ylabel("UMAP2",fontsize=16)
cbar = plt.colorbar(ticks = [0.001,5,10,14.91])
cbar.ax.set_yticklabels(['0','5','10','15'],fontsize=14)  # vertically oriented colorbar
plt.savefig("figures/qc/mt_pct.pdf")


sample_type = np.array(["_".join(x.split("_")[0:3]) for x in sequencing_metrics.index.values])
sequencing_metrics["sample"] = sample_type
sequencing_metrics.groupby("sample").mean().to_csv("results/qc/sample_metrics.csv")



plt.figure(figsize=(4,4))
plt.scatter( num_genes_detected, sequencing_metrics.UMIFM, s = 2, c = 'k' )
plt.xlabel("Number of genes detected",fontsize=16)
plt.ylabel("Number of unique molecules",fontsize=16)
cbar = plt.colorbar(ticks = [0.001,5,10,14.91])
'''

# Import  data
df = pd.read_csv("results/qc_info.csv")
tsne_coordinates = np.array(pd.read_csv("data/filtered/filtered_tsne_coordinates.csv",header=None))
sample_type = sample_type = np.array(pd.read_csv("data/filtered/filtered_labels.csv",header=None))
# Full t-SNE plot
plt.figure(num=None,figsize=(8,8) )
for type in np.unique(sample_type):
    idx = np.squeeze(sample_type == type)
    plt.scatter( tsne_coordinates[idx,0],tsne_coordinates[idx,1],s=2)
plt.legend(np.unique(sample_type))
plt.savefig("figures/tsne_filtered.pdf")


# Plot PCA of mitochondrial gene expression
plt.scatter( df.PC1, df.PC2, s = 2, c = df.average )
plt.colorbar()
plt.xlabel("Principal Component 1")
plt.ylabel("Principal Component 2")
plt.title("PCA of mitochondrial gene expression")
#plt.axis("square")
plt.savefig("figures/qc/mitochondrial_pca.pdf")

#heatmap = sns.clustermap( mito_expression[:,idx] )
#heatmap.savefig("figures/mitochondrial_heatmap.pdf")


# Find CDF of mitochondrial gene expression
thresholds = np.linspace(0, np.max(df.average), 100 )
num_cells = list()
for thresh in thresholds:
    num_cells.append( np.sum( df.average < thresh) )
  
# Plot distribution of mitochondrial gene expression
plt.subplot(211)
plt.hist( df.average, bins = 20,color = 'k' )
plt.hold()
plt.axvline(x = 3,linewidth=2,linestyle = 'dashed',color='r')
plt.xlabel("Average mitochondrial gene expression")
plt.ylabel("Number of cells")
plt.title("Histogram of mitochondrial gene expression")
plt.subplot(212)
plt.plot( thresholds, num_cells,'k')
plt.hold()
plt.axvline(x = 3,linewidth=2,linestyle = 'dashed',color='r')
plt.xlabel("Average mitochondrial gene expression")
plt.ylabel("Number of cells")
plt.title("Cumulative distribution of mitochondrial gene expression")
plt.savefig("figures/qc/mitochondrial_distribution.pdf")


# tSNE plot of mitochondrial expression
plt.figure(num=0,figsize=(4,4))
plt.scatter(tsne_coordinates[:,0],tsne_coordinates[:,1], s =2 , c = df.average )
plt.axis("square")
plt.colorbar()
plt.savefig("figures/qc/tsne_mitochondrial.pdf")

# Find CDF of number of genes detected
thresholds = np.linspace(0, np.max(df["Number of genes detected"]), 100 )
num_cells = list()
for thresh in thresholds:
    num_cells.append( np.sum( df["Number of genes detected"]< thresh) )
  
# Plot distribution of mitochondrial gene expression
plt.subplot(211)
plt.hist( df["Number of genes detected"], bins = 25, color = 'k' )
plt.hold()
plt.axvline(x = 750,linewidth=2,linestyle = 'dashed',color='r')
plt.xlabel("Number of genes detected")
plt.ylabel("Number of cells")
plt.title("Histogram of number of genes detected")
plt.subplot(212)
plt.plot( thresholds, num_cells,'k')
plt.hold()
plt.axvline(x = 750,linewidth=2,linestyle = 'dashed',color='r')
plt.xlabel("Number of genes detected")
plt.ylabel("Number of cells")
plt.title("Cumulative distribution of number of genes detected")
plt.savefig("figures/qc/detected_genes_distribution.pdf")



# tSNE plot of number of detected genes
plt.scatter(tsne_coordinates[:,0],tsne_coordinates[:,1], s =2 , c = df["Number of genes detected"] < 750 )
plt.axis("square")
plt.savefig("figures/qc/tsne_detected_genes.pdf")

# Plot number of genes detected vs. average mitochondrial expression
plt.scatter(df["Number of genes detected"],df.average, c = 'k', s =2 )
plt.xlabel("Number of genes detected")
plt.ylabel("Average mitochondrial gene expression")
plt.axvline(x = 750,linewidth=2,linestyle = 'dashed',color='r')
plt.axhline(y = 3,linewidth=2,linestyle = 'dashed',color='r')

plt.savefig("figures/qc/detected_v_mitochondrial.pdf")
'''