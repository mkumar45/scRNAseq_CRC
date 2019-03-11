from argparse import ArgumentParser
import matplotlib
from matplotlib import pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42 # Editable text when exporting pdfs
import numpy as np
import os
import pandas as pd
from pydpc import Cluster
import scipy as sp
import seaborn as sns
from sklearn.decomposition import PCA
import umap

#%%
def normalize_arcsinh(count_matrix,cofactor = 1000):
    """
    Description
    -----------   
    Normalize count matrix using inverse hyperbolic sine transform. 
    Parameters
    ----------    
    count_matrix : numpy array containing counts. row = cell, columsn = genes
    cofactor : cofactor for arcsinh transform. 
    Returns
    -------
    normalized_counts : numpy array normalized count matrix
    """
    total_counts = np.squeeze(np.sum(count_matrix,axis = 1))
    normalized_counts = np.arcsinh(count_matrix / total_counts[:,np.newaxis] * cofactor)
    return(normalized_counts)

def runPCA(count_matrix, n_comp):
    """
    Description
    -----------   
    run principal componenet analysis (PCA) on count matrix  
    Parameters
    ----------    
    count_matrix : numpy array containing counts. row = cell, columns = genes
    n_comp : number of componenets (see  sklearn.decomposition pca documentation)
    Returns
    -------
    pca_components : numpy array of principal componenets
    """
    _pca = PCA(n_components = n_comp)
    pca_components = _pca.fit(count_matrix).transform(count_matrix)
    return(pca_components)
    
def runUMAP( pca_components, n_neighbors = 10,min_dist = 0.5):
    """
    Description
    -----------   
    run uniform manifold approximation and projection (UMAP)
    Parameters
    ----------    
    pca_components : numpy array of principal componenets ( could be any matrix)
    n_neighbors : number of neighbors (see umap documentation)
    min_dist : minimum distance (see umap documentation)
    Returns
    -------
    UMAP_coordinates : numpy array UMAP projection
    """
    UMAP_coordinates = umap.UMAP(n_neighbors = n_neighbors, min_dist = min_dist, metric = 'correlation').fit_transform(pca_components)
    return(UMAP_coordinates)

def runDPC(X,x_cutoff,y_cutoff):
    """
    Description
    -----------   
    run density peak-based clustering (DPC)
    Parameters
    ----------    
    X : numpy array of to cluster( could be any matrix) row = cell, columns = umap coordinates
    x_cutoff : (see pydpc documentation)
    y_cutoff: (see pydpc documentation)
    Returns
    -------
    UMAP_coordinates : numpy array UMAP projection
    """
    DPC = Cluster(X.astype('float64'))
    DPC.assign(x_cutoff,y_cutoff)
    return(DPC)
    

#%%
def _argparse():
    """
    Description
    -----------   
    Parse arguments when running from command line
    Parameters
    ----------    
    input_file : path to .csv file to process
    gene_symbols : path to .csv file containing gene symbols
    Returns
    -------
    """
    parser = ArgumentParser('Import data')
    parser.add_argument("input_file", help = "input files to process",type=str,default="")
    parser.add_argument("gene_symbols", help = "path to gene symbols (columns of input file)",type=str,default="")

    return parser

def main(args):
    """
    Description
    -----------   
    Reads input count file and converts to sparse format for faster I/O downstream
    Files should be preprocessed toremove droplets with low counts (not cells). 
    Runs UMAP and DPC and filters sequencing metrics file for subsequent processing
    
    Parameters
    ----------    
    See _argparse
    Returns
    -------
    saves output files
        1) sparse matrix ( in data/qc/*_sparse.npz )
        2) sequencing metrics ( in data/qc/*_metrics.csv)
        3) UMAP coordinates ( in results/qc*_UMAP.csv)
        4) DPC cluster membership ( in results/qc*_DPC.csv)
        5) saves mouse gene symbols corresponding to columns of sparse matrix
        (data/mouse_gene_symbols.csv)
    """
    parser = _argparse()
    argp = parser.parse_args(args[1:])
    # Import data
    input_file = argp.input_file
    sample_name = os.path.basename(input_file).split(".")[0]    
    print("Procssing " + sample_name)
    sequencing_metrics = pd.read_csv("data/sequencing_metrics/" + sample_name + "_metrics.tsv", sep = "\t",
                                     index_col = 0 )
    cell_barcodes = pd.read_csv("data/sequencing_metrics/" + sample_name + "_index.csv",
                                index_col = 0, header = None )
    cell_barcodes = cell_barcodes.iloc[:,0].values
    gene_symbols = np.squeeze( np.array( pd.read_csv(argp.gene_symbols, header=None)) )
    
    # Drop reporter genes from count matrix and gene symbols
    df = pd.read_csv(input_file,header=None)
    df.columns = gene_symbols
    df = df.drop(["KL000001.1","KL000002.1"],axis = "columns")
    
    idx = np.argwhere( gene_symbols == "KL000001.1" )
    gene_symbols = np.delete(gene_symbols,idx)
    idx = np.argwhere( gene_symbols == "KL000002.1" )
    gene_symbols = np.delete(gene_symbols,idx)
    
    #%% Second pass to only import cells/mouse genes
    metrics_df = sequencing_metrics.loc[cell_barcodes] # get metrics for corresponding cells
    barcodes = [sample_name + "_" + barcode for barcode in cell_barcodes] # prepend barcode with sample name
    lib_size = np.sum(df,axis=1)
    

    # Normalize
    cofactor = 1000
    normalized_counts = np.array(np.arcsinh(df.div(lib_size,axis = 'rows')*cofactor))
    print("Running PCA")    
    pca_coordinates = runPCA(normalized_counts, n_comp = 100)
    print("Running UMAP")
    umap_coordinates = runUMAP( pca_coordinates, n_neighbors = 10,min_dist = 0.5)
    print("Running DPC")
    DPC = runDPC(umap_coordinates, x_cutoff = 0,y_cutoff= 0.5 )
    
    #%% Write data
    sparse_mat = sp.sparse.coo_matrix( df.values )
    sp.sparse.save_npz( "data/sparse/" + sample_name + ".npz", sparse_mat )
    pd.DataFrame(gene_symbols).to_csv("data/mouse_gene_symbols.csv",index=False, header = False )
    pd.DataFrame( metrics_df.values, columns=metrics_df.columns, index=barcodes).to_csv("data/sparse/" + sample_name + "_metrics.csv")
    # save umap coordinates and cluster ids
    pd.DataFrame( umap_coordinates, columns= ["UMAP1","UMAP2"], index = barcodes).to_csv("results/UMAP/" + sample_name + "_UMAP.csv",index=True)
    pd.DataFrame( DPC.membership, columns= ["Cluster"], index = barcodes).to_csv("results/DPC/" + sample_name + "_DPC.csv",index=True)
    
    
#%%
if __name__ == '__main__':
    from sys import argv
    exit(main(argv))