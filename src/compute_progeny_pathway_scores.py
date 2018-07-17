import numpy as np
import scipy as sp
from scipy import sparse
import pandas as pd
from argparse import ArgumentParser, FileType

import matplotlib.pyplot as plt
import seaborn as sns

def compute_pathway_score(scRNA_expression, gene_symbols,  PROGENy_coefficients ):
                        
        normalized_expression = sp.stats.zscore( scRNA_expression )
        nan_idx = np.logical_not(np.all(np.isnan( normalized_expression ), axis = 0 ))
        
        normalized_expression = normalized_expression[:,nan_idx]
        gene_symbols = gene_symbols[nan_idx]
        
        # Pathway genes only 
        idx = [x in PROGENy_coefficients.index.values for x in gene_symbols.index.values]
        expression_subset = normalized_expression[:,idx]

        
        coeff_subset = PROGENy_coefficients.loc[ gene_symbols.index.values[idx] ] # only use coefficients found in expression
        coeff_subset = coeff_subset.reindex( gene_symbols.index.values[idx] ) # Reorder for multiplication

        pathway_scores = np.matmul( expression_subset, coeff_subset)
        df = pd.DataFrame(pathway_scores,columns=coeff_subset.columns)
        return(df)
        
def plot_score_heatmap(pathway_scores, cell_type_labels ):
    
        sns.heatmap( pathway_scores )
        plt.figure()
        sns.heatmap( sp.stats.zscore(pathway_scores) )
        
        cell_types = np.unique(cell_type_labels)
        num_cell_types = len(cell_types)

        pathway_average_scores = np.empty((num_cell_types,pathway_scores.shape[1],))
        
        for ct in range(num_cell_types):
            ct_idx = cell_type_labels.values == cell_types[ct]
            pathway_average_scores[ct,:] = np.mean( pathway_scores.loc[np.squeeze(ct_idx),:], axis=0)
                    
        
        plt.figure()
        sns.heatmap( pd.DataFrame(pathway_average_scores, columns = pathway_scores.columns, index = cell_types))
        plt.figure()
        sns.heatmap(pd.DataFrame(sp.stats.zscore(pathway_average_scores), columns = pathway_scores.columns, index = cell_types))
        plt.show()
def _argparse():
    """
    Description
    -----------   
    Parameters
    ----------    
    Returns
    -------
    """
    parser = ArgumentParser('Classify cell-types using scRNA expression and predefined markers')
    parser.add_argument("input_file", help = "sparse expression matrix",
                        type=str,default="")
    parser.add_argument("gene_names", help = "gene names corresponding to rows of expression matrix",
                        type=str,default="")
    parser.add_argument("cell_type_labels", help = "directory to save output files",
                        type=str,default="")
    parser.add_argument("PROGENy",
                        help = "file containing coefficients for PROGENy",
                        type=str,default="")
    parser.add_argument("output_dir", help = "directory to save output files",
                        type=str,default="")
    parser.add_argument("-g","--group",
                        help = "Number of partitions to use for k-fold cross validation"
                        ,type=str,default=None)

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
    # Parse inputs
    parser = _argparse()
    argp = parser.parse_args(args[1:])
    # Read data
    scRNA_expression = sp.sparse.load_npz( argp.input_file )
    gene_names = pd.read_csv(argp.gene_names,index_col = 0,header=None)
    progeny_coeff = pd.read_csv( argp.PROGENy, index_col = 0 )
    cell_type_labels = pd.read_csv(argp.cell_type_labels,header=None)
    #scRNA_expression = sp.sparse.load_npz( "data/combined_sparse.npz" )
    #gene_names = pd.read_csv("data/combined_gene_names.csv",index_col = 0,header=None)
    #RL_list = import_RL( "data/mouse_receptor_ligand.csv" )
    #cell_type_labels = pd.read_csv("results/predicted_cell_types.csv",header=None)
    
    
    PROGENy_coefficients = pd.read_csv( "data/PROGENy_mouse_model_v2.csv", index_col = 0 )
    cell_type_labels = pd.read_csv("results/predicted_cell_types.csv",header=None)
    
    scRNA_expression = np.array( sp.sparse.load_npz("data/combined_sparse.npz").todense() )
    gene_symbols = pd.read_csv("data/combined_gene_names.csv" ,index_col = 0,header=None)

    sample_labels = pd.read_csv("data/combined_labels.csv",header = None)
    sample_type = [ x.split("_")[0] for x in sample_labels.iloc[:,0] ]
    
    for sample in np.unique(sample_type):
        print(sample)
        idx = [ typ == sample for typ in sample_type]
        pathway_scores = compute_pathway_score(scRNA_expression[idx,:], gene_symbols, PROGENy_coefficients )
        plot_score_heatmap(pathway_scores, cell_type_labels[idx] )
        pathway_scores.to_csv( "results/"+sample+"_PROGENy_scores.csv")

    
    pathway_scores = compute_pathway_score(scRNA_expression, gene_symbols, PROGENy_coefficients )
    pathway_scores.to_csv( "results/PROGENy_scores.csv")

    
    #idx = [ typ == label for typ in sample_type]   
    
    
    
    #df.to_csv( "results/PROGENy_scores.csv")
    
if __name__ == '__main__':
    from sys import argv
    exit(main(argv))