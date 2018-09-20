# -*- coding: utf-8 -*-
"""
Created on Wed Jun 27 12:04:27 2018

@author: Manu

"""

from argparse import ArgumentParser, FileType
import os
import glob
import numpy as np
import pandas as pd
import scipy as sp
from scipy import sparse

def _argparse():
    parser = ArgumentParser('Combine tumor expression data')
    parser.add_argument("processed_directory", help = "directory of processed files to combine",type=str,default="")
    return parser

def main(args):
    """
    Description
    -----------   
    Parameters
    ----------    
    Returns
    -------
    """
    parser = _argparse()
    argp = parser.parse_args(args[1:])
    
    sparse_files = glob.glob(argp.processed_directory + "*.npz")
    gene_files = glob.glob(argp.processed_directory + "*gene_names.csv")
    metric_files = glob.glob(argp.processed_directory + "*processed_metrics.csv")
        
    num_samples = len(sparse_files)
    gene_names = []
    sequencing_metrics = []
    
    for f in range(num_samples):
        gene_names.append(  pd.read_csv(gene_files[f],header = None,index_col = 0) )
        sequencing_metrics.append( pd.read_csv( metric_files[f], index_col=0))
        
    combined_metrics = pd.concat( sequencing_metrics )
    merged_genes = pd.concat(gene_names, axis = 1,join = "inner") # find interaction of gene names
    csc_matrices = []
    sample_names = []
    # Subset each matrix
    for f in range(num_samples):
        file_name = os.path.splitext(os.path.basename(sparse_files[f]))[0]
        idx = [gene in merged_genes.index.values for gene in gene_names[f].index.values]
        sparse_mat =  sp.sparse.load_npz( sparse_files[f] )
        csc_mat = sparse_mat.tocsc()
        csc_matrices.append( csc_mat[:,np.squeeze(np.nonzero(idx))] )
        sample_names.append([file_name] * csc_mat.shape[0])
    # Combine sparse matrices
    combined_sparse = sp.sparse.vstack( csc_matrices, format = "coo")
    sp.sparse.save_npz("data/combined/combined_sparse.npz",combined_sparse)
    
    merged_genes.to_csv("data/combined/combined_gene_names.csv", columns = [], header = False )
    combined_metrics.to_csv("data/combined/combined_metrics.csv")
    
    
    sample_labels = [ item for sublist in sample_names for item in sublist]
    pd.DataFrame( sample_labels,index=combined_metrics.index.values ).to_csv("data/combined/combined_labels.csv",header=False)

if __name__ == '__main__':
    from sys import argv
    exit(main(argv))