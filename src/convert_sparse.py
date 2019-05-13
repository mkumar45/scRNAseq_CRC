from argparse import ArgumentParser
import numpy as np
import pandas as pd
import scipy as sp
from scipy import sparse

def _argparse():
    """
    Description
    -----------   
    Parse arguments when running from command line
    Parameters
    ----------    
    input_file : path to .npz file to process
    Returns
    -------
    """
    parser = ArgumentParser('Import sparse matrix and write to dense format')
    parser.add_argument("input_file", help = "sparse_matrix (npz) to save as dense",type=str,default="")
    parser.add_argument("output_file", help = "filename for output",type=str,default="")

    return parser

#%%
def main(args):
    """
    Description
    -----------   
    Reads input sparse file and save in dense format as .csv
    Parameters
    ----------    
    See _argparse
    Returns
    -------
    saves output files
        1) dense matrix ( in data/qc/*_sparse.npz )
    
    """
    parser = _argparse()
    argp = parser.parse_args(args[1:])
    sparse_matrix = np.array(sp.sparse.load_npz( argp.input_file ).todense())
    
    pd.DataFrame(sparse_matrix).to_csv(argp.output_file,index = False, header = False)
    
    
    
#%%
if __name__ == '__main__':
    from sys import argv
    exit(main(argv))
    