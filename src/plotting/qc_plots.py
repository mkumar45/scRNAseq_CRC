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



# Import  data
df = pd.read_csv("results/qc_info.csv")
tsne_coordinates = np.array(pd.read_csv("results/combined_sparse_tsne.csv",header=None))


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
