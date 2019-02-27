from argparse import ArgumentParser
import matplotlib
from matplotlib import pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42 # Editable text when exporting pdfs
import numpy as np
import os
import pandas as pd
from progress.bar import Bar
from pydpc import Cluster
import re
import scipy as sp
from sklearn.decomposition import PCA
import sys
import umap


#%%
def find_inflection(lib_size, plot_name = None ):
    """
    Description
    -----------   
    Permissive filtering of droplets to remove debris and only use droplets
    containing cells. 
    Parameters
    ----------    
    lib_size : numpy array containing library size (sum of counts) for each
    barcode
    plot_name : name of outplot if desired. Otherwise no plot is saved. 
    Returns
    -------
    index of inflection point to use as threshold for removing barcodes. 
    Assumes barcodes are sorting by library size from least to greatest. 
    """
    num_barcodes = len(lib_size)
    sorted_lib_size = np.sort( lib_size )
    print(num_barcodes)
    cumsum = np.cumsum( sorted_lib_size )        
    x_vals = np.arange(0,num_barcodes)
    # Secant line
    secant_coef=cumsum[num_barcodes-1]/num_barcodes
    secant_line=secant_coef*x_vals
    # Distance b/w secant and cumulative sum
    secant_dist=abs(cumsum-secant_line)
    max_dist=np.where(np.max(secant_dist)==secant_dist)[0][0] # inflection point
    
    if plot_name:
        plt.figure(figsize=(8,8))
        plt.plot(np.array(cumsum), label="Library size cumulative sum", linewidth=2)
        plt.plot(np.array(secant_line), linestyle='--',label="Secant Line", linewidth=2)
        plt.plot(np.array(secant_dist), label="CS-SL Distance", linewidth=2)
        
        plt.plot(max_dist, np.max(secant_dist),".",markersize = 24, c='red',label="Inflection Point")
        plt.axvline(x=max_dist, ymax = np.max(secant_dist)/np.max(cumsum),c='red',linestyle='--',linewidth=2)
        plt.legend()
        plt.xlabel("Barcode Rank")
        plt.ylabel("Counts")
        #plt.savefig(plot_name)        
        
    print("Inflection point at {}".format(max_dist))
    return( max_dist )

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
    input_file : path to .tsv file to process
    Returns
    -------
    """
    parser = ArgumentParser('Import data')
    parser.add_argument("input_file", help = "input files to process",type=str,default="")
    return parser

#%%
def main(args):
    """
    Description
    -----------   
    Reads input count file and removes droplets with low counts (not cells).
    Removes human genes from count matrix. 
    Runs UMAP and DPC on filtered droplets and saves data in sparse format for
    faster I/O downstream
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
    input_file = argp.input_file
    #input_file = "data/counts/APC_adjacent_S3.tsv"
    sample_name = os.path.basename(input_file).split(".")[0]    
    print("Procssing " + sample_name)
    sequencing_metrics = pd.read_csv("data/sequencing_metrics/" + sample_name + "_metrics.tsv", sep = "\t",
                                     index_col = 0 )
    
    #%% First pass through the data to count library size (sum of counts) for each barcode
    # One line at a time to avoid having to read the entire file at once
    lib_size = list()
    barcodes = list()
    
    total_file_size = os.path.getsize( input_file )
    progress_bar = Bar("First pass read:",max=100)
    read_bytes = 0 

    with open( input_file ) as f:
        header = f.readline() # save gene names for later filtering
        read_bytes += sys.getsizeof(header)
        split_header = header.split("\t")
        gene_symbols = split_header[1:]

        for line in f:
            read_bytes += sys.getsizeof(line)
            split_line = line.split("\t")
            barcodes.append( split_line[0].strip() )
            counts = np.array(split_line[1:],dtype=int)
            lib_size.append( np.sum(  counts ) )
            
            if (read_bytes > total_file_size/100):
                read_bytes = 0
                progress_bar.next()
    progress_bar.next()
    progress_bar.finish()
    f.close()
    
    #%% 
    # Find barcodes with high library cells (cells)
    inflection_idx = find_inflection(lib_size, plot_name = "figures/qc/" + sample_name + "_inflection.pdf")    
    sort_idx = np.argsort( lib_size )
    barcodes_to_keep = np.array(barcodes)[sort_idx[inflection_idx:]]
        
    # Use only mouse genes
    regex = re.compile("(?=.*[a-z])") # regular expression to find any lowercase
    mouse_idx = [ bool( re.match( regex, symbol ) ) for symbol in gene_symbols]
    mouse_symbols = np.array(gene_symbols)[mouse_idx]
    
    #%% Second pass to only import cells/mouse genes
    mouse_counts = np.empty((len(barcodes_to_keep),np.sum(mouse_idx)))
    barcodes = list()
    n = 0 
    
    progress_bar = Bar("Second pass read",max = len(barcodes_to_keep))
    with open( input_file ) as f:
        next(f) # skip header
        for line in f:
            split_line = line.split("\t")
            barcode = split_line[0].strip() 
            if barcode in barcodes_to_keep:
                barcodes.append( split_line[0].strip() ) # first entry is barcode
                counts = np.array(split_line[1:],dtype=int) # everything else is counts 
                mouse_counts[n,:] = (counts[mouse_idx])
                n = n + 1
                progress_bar.next()
                if n == len(barcodes_to_keep):
                    break
    progress_bar.finish()
    f.close()
    
    
    metrics_df = sequencing_metrics.loc[barcodes] # use only corresponding sequencing metric
    barcodes = [sample_name + "_" + barcode for barcode in barcodes] # prepend barcode with sample name
    count_df = pd.DataFrame(mouse_counts, index = barcodes, columns = mouse_symbols )
    lib_size = np.sum(count_df,axis=1)
    

    # Normalize
    cofactor = 1000
    normalized_counts = np.array(np.arcsinh(count_df.div(lib_size,axis = 'rows')*cofactor))
    print("Running PCA")    
    pca_coordinates = runPCA(normalized_counts, n_comp = 100)
    print("Running UMAP")
    umap_coordinates = runUMAP( pca_coordinates, n_neighbors = 10,min_dist = 0.5)
    print("Running DPC")
    DPC = runDPC(umap_coordinates, x_cutoff = 0,y_cutoff= 0.5 )
    
    #%% Write data
    sparse_mat = sp.sparse.coo_matrix( count_df.values )
    sp.sparse.save_npz( "data/qc/" + sample_name + "_sparse", sparse_mat )
    pd.DataFrame( mouse_symbols ).to_csv("data/mouse_gene_symbols.csv",header=None,index=False)
    pd.DataFrame( metrics_df.values, columns=metrics_df.columns, index=barcodes).to_csv("data/qc/" + sample_name + "_metrics.csv")
    # save umap coordinates and cluster ids
    pd.DataFrame( umap_coordinates, columns= ["UMAP1","UMAP2"], index = barcodes).to_csv("results/qc/" + sample_name + "_UMAP.csv",index=True)
    pd.DataFrame( DPC.membership, columns= ["Cluster"], index = barcodes).to_csv("results/qc/" + sample_name + "_DPC.csv",index=True)
    
    
#%%
if __name__ == '__main__':
    from sys import argv
    exit(main(argv))
    