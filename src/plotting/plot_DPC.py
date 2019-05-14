from argparse import ArgumentParser
import matplotlib
from matplotlib import pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42 # Editable text when exporting pdfs
import numpy as np
import pandas as pd
import scipy as sp
import seaborn as sns




def plot_umap( umap_coordinates, cluster_id ):
    fig = plt.figure(figsize=(5,4))
    plt.scatter(umap_coordinates.UMAP1, umap_coordinates.UMAP2, s = 10, c = cluster_id.Cluster, cmap="tab20")
    plt.colorbar()
    # Overlay cluster label at center of each cluster
    for n in np.unique(cluster_id):
        member_idx = cluster_id.Cluster == n
        x_center = np.mean( umap_coordinates.UMAP1[member_idx] )
        y_center = np.mean( umap_coordinates.UMAP2[member_idx])            
        plt.text(x_center,y_center, str(n), fontsize = 18 )
    plt.title("UMAP Clustering")
    plt.xticks([])
    plt.yticks([])
    plt.xlabel("UMAP1",fontsize=16)
    plt.ylabel("UMAP2",fontsize=16)
    
    

plot_umap( umap_coordinates, DPC )
plt.savefig("figures/DPC_clustering.pdf")