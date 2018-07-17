from argparse import ArgumentParser, FileType
from sys import stderr, stdin, stdout
import numpy as np
import scipy as sp
import pandas as pd
from bhtsne import bhtsne


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
    parser.add_argument("output_file", help = "name of output file",type=str,default="")
    parser.add_argument('-n', '--num_genes', type=int, default=None)


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
    sparse_mat = sp.sparse.load_npz( argp.input_file )
    full_mat = sparse_mat.todense()


    if argp.num_genes is None: 
        result = bhtsne.run_bh_tsne( full_mat )
    else:
        std_dev = np.squeeze( np.array( np.std( full_mat, axis = 0 ) ) )
        sort_idx = np.squeeze( np.array( np.argsort(std_dev) ) ) 
        expression_subset = full_mat[:, sort_idx[-argp.num_genes:]]
        
        result = bhtsne.run_bh_tsne( expression_subset )
        
        
        
        
    np.savetxt(argp.output_file,result,delimiter = ",",fmt="%10.5f")


if __name__ == '__main__':
    from sys import argv
    exit(main(argv))
    