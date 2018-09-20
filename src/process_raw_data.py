from argparse import ArgumentParser, FileType
from sys import stderr, stdin, stdout
import numpy as np
import scipy as sp
from scipy import sparse
import pandas as pd

def _argparse():
    """
    Description
    -----------   
    Parameters
    ----------    
    Returns
    -------
    """
    parser = ArgumentParser('Combine tumor expression data')
    parser.add_argument("input_file", help = "input files to process",type=str,default="")
    parser.add_argument("metric_file", help = "name of the file containing sequencing metrics", type=str,default="")
    parser.add_argument("sparse_file", help = "name of output file for sparse matrix",type=str,default="")
    parser.add_argument("gene_file", help = "name of output file for gene names",type=str,default="")
    parser.add_argument("metric_out", help = "name of output file for metrics",type=str,default="")

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
    # Import expression data
    df = pd.read_csv( argp.input_file, header=0, index_col=0 )
    # Import sequencing metrics
    metrics = pd.read_csv( argp.metric_file, delimiter="\t",index_col = 0 )
    # Merge based on cell barcodes
    concat_df = pd.concat([df, metrics], axis = 1, join = "inner" )

    
    if "tSNE1" in df.columns:
        df = df.drop( labels = "tSNE1",axis = "columns") # pre-processed data has t-SNE columns
        
    if "tSNE2" in df.columns:
        df = df.drop( labels = "tSNE2",axis = "columns") # pre-processed data has t-SNE columns
        
    if "ID" in df.columns:
        df = df.drop( labels = "ID",axis = "columns") # pre-processed data has t-SNE columns
        
    # Reseparate sequencing and expression data
    df = concat_df.drop( metrics.columns.values, axis = "columns" )
    metrics = concat_df.drop( df.columns.values, axis = "columns")
        
    
    # Convert expression data to sparse format (for faster file I/O)
    coo_mat = sp.sparse.coo_matrix( df )
    sp.sparse.save_npz( argp.sparse_file ,coo_mat )
    # Write gene names
    pd.DataFrame(df.columns.values ).to_csv( argp.gene_file,index= False, header = False )
    # Write sequencing metrics
    metrics.to_csv( argp.metric_out )

if __name__ == '__main__':
    from sys import argv
    exit(main(argv))
    