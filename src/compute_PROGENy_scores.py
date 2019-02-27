import numpy as np
import scipy as sp
from scipy import sparse, stats
import pandas as pd
from argparse import ArgumentParser, FileType

#%% Normalize count matrix using arcsinh transformation. Assumes each row is a sample/cell and columns correspond to genes
def normalize_arcsinh(count_matrix,cofactor = 1000):
    total_counts = np.squeeze(np.sum(count_matrix,axis = 1))
    normalized_counts = np.arcsinh(count_matrix / total_counts[:,np.newaxis] * cofactor)
    return(normalized_counts)
    
    
    
def compute_pathway_score(scRNA_expression, gene_symbols,  PROGENy_coefficients ):
                        
        normalized_expression = sp.stats.zscore( scRNA_expression )
        nan_idx = np.logical_not(np.all(np.isnan( normalized_expression ), axis = 0 ))
        
        normalized_expression = normalized_expression[:,nan_idx]
        gene_symbols = gene_symbols[nan_idx]
        
        # Pathway genes only 
        idx = [x in PROGENy_coefficients.index.values for x in gene_symbols]
        expression_subset = normalized_expression[:,idx]

        
        coeff_subset = PROGENy_coefficients.loc[ gene_symbols[idx] ] # only use coefficients found in expression
        coeff_subset = coeff_subset.reindex( gene_symbols[idx] ) # Reorder for multiplication

        pathway_scores = np.matmul( expression_subset, coeff_subset)
        df = pd.DataFrame(pathway_scores,columns=coeff_subset.columns)
        return(df)
        
        
def cell_type_specific_expression( scRNA_expression, gene_symbols, return_genes, cell_type_labels ):
    """
    Description
    -----------
    Return cell-type specific expression for each gene in list of genes     
    Parameters
    ----------       
    scRNA_expression : sparse matrix
        expression matrix with rows corresponding to cells and columns to genes
    gene_symbols : 
        corresponding gene_symbols for each column of scRNA_expression
    return_genes :
        genes to return cell-types expression
    cell_type_labels:
        cell-type labels for each row of scRNA_expression matrix
    Returns
    -------
    """
    
    
    cell_types = np.unique( cell_type_labels )
    num_cell_types = len(cell_types)
    num_genes = len( return_genes )
    
    gene_idx = np.empty((num_genes,1))
    unique_symbols = np.unique( return_genes )
    for sym in unique_symbols:
        
        gene_idx[ return_genes == sym ] = np.nonzero( gene_symbols.index.values == sym )[0]
        
    ct_expression = [[] for n in range(num_cell_types)]

    for ct in range(num_cell_types):
        ct_idx = cell_type_labels.values == cell_types[ct]
        
        ct_expression[ct] = scRNA_expression[ np.squeeze(ct_idx), gene_idx.astype(int)] # cast to int for indexing
    return(ct_expression)
        
        
def compute_randomized_scores(scRNA_expression, gene_symbols,  PROGENy_coefficients, cell_type_labels):
    
    normalized_expression = sp.stats.zscore( scRNA_expression )
    idx = [x in PROGENy_coefficients.index.values for x in gene_symbols.index.values]
    expression_subset = normalized_expression[:,idx]
    
    cell_type_expression = cell_type_specific_expression( scRNA_expression, gene_symbols, PROGENy_coefficients.index.values,cell_type_labels)
    
    
    np.random.permutation()
    
        
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
    parser.add_argument("cell_type_labels", help = "file containing cell type labels",
                        type=str,default="")
    parser.add_argument("PROGENy",
                        help = "file containing coefficients for PROGENy",
                        type=str,default="")
    parser.add_argument("output_dir", help = "directory to save output files",
                        type=str,default="results/PROGENy/PROGENy_scores.csv")
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
    scRNA_expression = np.array( sp.sparse.load_npz( argp.input_file ).todense() )
    normalized_expression = normalize_arcsinh( scRNA_expression )
    gene_symbols = np.squeeze( np.array( pd.read_csv(argp.gene_names, header=None)) )
    cell_type_labels = pd.read_csv(argp.cell_type_labels,header=None)
    PROGENy_coefficients = pd.read_csv( argp.PROGENy, index_col = 0 )
    
    '''
    for sample in np.unique(sample_type):
        print(sample)
        idx = [ typ == sample for typ in sample_type]
        pathway_scores = compute_pathway_score(scRNA_expression[idx,:], gene_symbols, PROGENy_coefficients )
        plot_score_heatmap(pathway_scores, cell_type_labels[idx] )
        pathway_scores.to_csv( "results/"+sample+"_PROGENy_scores.csv")
    '''
    
    pathway_scores = compute_pathway_score(normalized_expression, gene_symbols, PROGENy_coefficients )
    pathway_scores = pathway_scores.set_index( cell_type_labels.index.values )
    pathway_scores.to_csv( argp.output_dir )
    
if __name__ == '__main__':
    from sys import argv
    exit(main(argv))