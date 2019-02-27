import numpy as np
import scipy as sp
import pandas as pd

from argparse import ArgumentParser
from sklearn import mixture
from sklearn.model_selection import KFold

#%% Normalize count matrix using arcsinh transformation. Assumes each row is a sample/cell and columns correspond to genes
def normalize_arcsinh(count_matrix,cofactor = 1000):
    total_counts = np.squeeze(np.sum(count_matrix,axis = 1))
    normalized_counts = np.arcsinh(count_matrix / total_counts[:,np.newaxis] * cofactor)
    return(normalized_counts)


#%% Create subset of data with only cell-type marker expression data
def create_training_data( scRNA_expression, gene_symbols, marker_data, GMMs ):
    """
    Description
    -----------   
    Parameters
    ----------
    scRNA_expression : numpy array
        count matrix with each column is a gene and each row is a single cell
    gene_symbols :  numpy array
        gene symbols corresonding to columsn of scRNA_expressoin
    marker_data : pandas dataframe    
        Dataframe containing user specified marker genes defining cell-types.
        Each row is a cell-type and each column is a marker gene. Entries contain
        AND,OR,NOT
    GMMs: list
        List containing Gaussian mixture models for each marker gene. Fit using the 
        fit_GMM function
    
    Returns
    -------
    training_data: training expression data in COO sparse format
        subset of input data to use for training classifier
    training_labels : numpy array
        assigned cell-type based on marker expression to use for training classifier
    training_idx : numpy array
        boolean array indicating which cells from full expression are used for training
    """
    
    # Subset expression data to only marker expressoin
    marker_idx = [gene in marker_data.columns.values for gene in gene_symbols]
    marker_expression = scRNA_expression[:,marker_idx]
    
    num_cell_types = marker_data.shape[0]
    num_cells = scRNA_expression.shape[0]
    is_cell_type = np.zeros((num_cell_types,num_cells),dtype=bool)
    
    # Assign cluster 
    cluster_idx = np.array([ GMMs[n].predict( marker_expression[:,n].reshape(-1,1) ) for n in range(marker_expression.shape[1]) ])
    max_idx = np.array([ np.argmax(GMMs[x].means_)  for x in range(len(GMMs))])
    min_idx = np.array([ np.argmin(GMMs[x].means_)  for x in range(len(GMMs))])
    
    # Loop through each cell type and determine if cell matches marker profile
    for ct in range(num_cell_types):
        ct_markers = marker_data.iloc[ct, :]

        AND_markers = ct_markers.index.values[ct_markers == "AND"]
        OR_markers = ct_markers.index.values[ct_markers == "OR"]
        NOT_markers = ct_markers.index.values[ct_markers == "NOT"]
        
        AND_idx = [ x in AND_markers for x in marker_data.columns.values]
        OR_idx = [ x in OR_markers for x in marker_data.columns.values]
        NOT_idx = [ x in NOT_markers for x in marker_data.columns.values]
        
        AND_flag = np.all(cluster_idx[AND_idx,] == max_idx[AND_idx,np.newaxis],axis=0) if np.any(AND_idx) else np.ones((num_cells,),dtype=bool)
        OR_flag = np.any(cluster_idx[OR_idx,] == max_idx[OR_idx,np.newaxis],axis=0) if np.any(OR_idx) else np.ones((num_cells,),dtype=bool)
        NOT_flag = np.all(cluster_idx[NOT_idx,] == min_idx[NOT_idx,np.newaxis],axis=0 ) if np.any(NOT_idx) else np.ones((num_cells,),dtype=bool)
        

        is_cell_type[ct,:] = np.all( np.stack( [AND_flag,OR_flag,NOT_flag] ),axis=0 )                       

    # Only use cells assigned to a single cell type
    training_idx = np.sum( is_cell_type, axis=0) == 1
    training_data = scRNA_expression[training_idx,:]
    # Find index of unique cell type and add cell type name to column data
    cell_type_idx = np.where( np.transpose(is_cell_type[:,training_idx]) )
    training_labels = marker_data.index.values[cell_type_idx[1]]
    
    return(training_data,training_labels,training_idx)
       
#%% Compute expression thresholds for each marker using otsu's method
#   Assumes two groups (high expressing and low expressing)
def fit_GMM( marker_expression, num_components = 5, k_fold = 5 ):
    """
    Description
    -----------   
    Fit gaussian mixture models to gene expression data using KFold cross validation
    and select the best model based on bayesian information criteria (BIC). 
    Parameters
    ----------    
    marker_expression : numpy array
        count matrix with each column is a marker gene and each row is a single cell
    num_components : int 
        specifies maximimum number of mixture componenets to test
    k_fold : int
        specify number of folds for cross-validation
    Returns
    -------
    best_GMM : list
        returns a list containing one Gaussian mixture model per marker gene (row of
                                                                              input data) 
    GMM_dict : dictionary
        contains Gaussian mixture model parameters (means, covariances, and weights)
    BIC_dict : dictionary
         contains BIC values from Gaussian mixture modeling fitting
         (means and standard error of the means)    
    """
    num_markers = marker_expression.shape[1]
    n_components_range = range(1, num_components+1)
    CV = KFold( n_splits = k_fold )

    BIC = np.empty((num_markers,len(n_components_range),k_fold))
    # Fit models 
    for num in range(num_markers):
        for num_comp in n_components_range:
            fold_bic = []
            X =  marker_expression[:,num].reshape(-1,1)
            for train,test in CV.split( X ):
                # Fit a Gaussian mixture with EM
                GMM = mixture.GaussianMixture(n_components=num_comp)
                GMM.fit( X[train] )
                fold_bic.append(GMM.bic(X[test]))
                
            BIC[num,num_comp-1,:] = np.array( fold_bic )
            
    # Select model with lowest BIC within 1SE of minimum
    mean_BIC = np.mean( BIC, axis = 2 )
    sem_BIC = sp.stats.sem( BIC, axis = 2)
    
    min_idx = np.argmin(mean_BIC,axis = 1 )
    threshold = mean_BIC.min(axis = 1) + sem_BIC[np.arange( sem_BIC.shape[0] ), min_idx ]
    num_components = np.argmax( mean_BIC < threshold[:,np.newaxis],axis=1 ) + 1 # +1 due to zero-indexing
    
    # Fit model using all cells
    best_GMM = []
    for num in range(num_markers):
        X =  marker_expression[:,num].reshape(-1,1)
        GMM = mixture.GaussianMixture(n_components=num_components[num])
        GMM.fit( X )
        best_GMM.append(GMM)
        
    GMM_dict = { "means" : [x.means_ for x in best_GMM],
                 "covariances" : [x.covariances_ for x in best_GMM],
                 "weights" : [x.weights_ for x in best_GMM]}
    BIC_dict = { "mean" : mean_BIC,
                 "sem" : sem_BIC}
    return(best_GMM, GMM_dict, BIC_dict )
    



def _argparse():
    """
    Description
    -----------
    Parse command line arguments
    Parameters
    ----------    
    See individual arguments
    Returns
    -------
    """
    parser = ArgumentParser('Classify cell-types using scRNA expression and predefined markers')
    parser.add_argument("input_file", help = "preprocessed input file to classify",
                        type=str,default="")
    parser.add_argument("gene_names", help = "file containing gene names",
                        type=str,default="")
    parser.add_argument("marker_file",
                        help = "file containing marker genes for classification",
                        type=str,default="")
    parser.add_argument("-k","--k_fold",
                        help = "Number of partitions to use for k-fold cross validation"
                        ,type=int,default=5)
    parser.add_argument("-m","--n_components",
                        help = "Number of mixture componenets to fit for Gaussian Mixture models"
                        ,type=int,default=5)
    parser.add_argument("-n","--num_genes",
                        help = "Only use the n highest genes by standard deviation"
                        ,type=int ,default=1000)
    parser.add_argument("-pca","--pca_variance",
                        help = "% of PCA variance to keep when defining training inputs"
                        ,type=float  ,default=0.90)

    return parser
        

#%% Compute expression thresholds for each marker using otsu's method
#   Assumes two groups (high expressing and low expressing)
def main(args):
    """
    Description
    -----------   
    Create a training-data set of gene expression and high-confidence cell type
    labels based on the manually defined cell-type markers 
    
    Parameters
    ----------
    see argparse_              
    
    Returns
    -------
    saves output files
        1)  sparse matrix containg counts for cell in training data 
            ( results/classification/sparse_training_counts.npz )
        2)  pandas dataframe containing cell type label for each cell in the 
            training data set
            ( results/classification/training_lables.csv )
    """
    # Parse inputs
    parser = _argparse()
    argp = parser.parse_args(args[1:])
    # Read data
    scRNA_expression = sp.sparse.load_npz( argp.input_file )
    scRNA_expression = np.array( scRNA_expression.todense() )
    normalized_expression = normalize_arcsinh(scRNA_expression)
    
    sequencing_metrics = pd.read_csv("data/combined/processed_metrics.csv",index_col = 0)
    gene_symbols = np.squeeze( np.array( pd.read_csv(argp.gene_names, header=None)) )
    marker_df = pd.read_csv( argp.marker_file, index_col = 0 )
    
    # Use only markers are present in expression data
    idx = [ elem in gene_symbols for elem in marker_df.columns.values ]
    marker_subset = marker_df.loc[:,idx] 
    marker_subset = marker_subset.dropna(axis = 0, how = "all") # drop cell type if no markers remain
    marker_subset = marker_subset.reindex_axis(sorted(marker_subset.columns), axis=1)

    # Gaussian Mixture models 
    idx = [ gene in marker_subset.columns.values for gene in gene_symbols ]
    marker_expression = normalized_expression[:,idx]    
    GMMs, GMM_dict, BIC_dict = fit_GMM(marker_expression, 5, 5)
    training_data, training_labels, training_idx = create_training_data(scRNA_expression, gene_symbols, marker_subset, GMMs )
    
    # Save outputs    
    sp.sparse.save_npz("results/classification/sparse_training_counts.npz",
                       sp.sparse.coo_matrix( scRNA_expression[training_idx,:]) )
    pd.DataFrame( training_labels,
                  index = sequencing_metrics.index.values[training_idx],
                  columns = ["training_label"]).to_csv(
                    "results/classification/training_labels.csv")

if __name__ == '__main__':
    from sys import argv
    exit(main(argv))
    