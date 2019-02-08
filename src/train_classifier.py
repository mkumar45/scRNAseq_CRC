from argparse import ArgumentParser
import matplotlib
from matplotlib import pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42 # Editable text when exporting pdfs
import numpy as np
import pandas as pd
import scipy as sp
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import confusion_matrix

#%% Normalize count matrix using arcsinh transformation. Assumes each row is a sample/cell and columns correspond to genes
def normalize_arcsinh(count_matrix,cofactor = 1000):
    total_counts = np.squeeze(np.sum(count_matrix,axis = 1))
    normalized_counts = np.arcsinh(count_matrix / total_counts[:,np.newaxis] * cofactor)
    return(normalized_counts)

#%% Train a classifer using training_object generated by create_training_data
def train_LR( training_data, training_labels, num_genes = 1024, pct_var = 0.90, save_results = True ):
    """
    Description
    -----------
    Trains logistic regression classification model.
        1) Select subset of features based on standard deviation
        2) Run PCA using only selected features
        3) Use PCA components as inputs into logistic regression classifier
    Calls sklearn logistic regression, which has cross-validation built-in

    Parameters
    ----------  
    training_data : numpy array
        contains count data
    training_labels : numpy array
        cell type labels corresponding to each row of training data
    num_genes : int
        number of genes to select
    pct_var : int
        percentage of variance to explain by PCA. Determines number of components
        used as inputs to logistic regression
    save_results : boolean
        specify whether to save cross validation results
    See individual arguments
    
    Returns
    -------
    a dictionary containing handles to functions for each step in the classification
    procedure (i.e. normalization, feature selection, pca transformation, and prediction)
    """
    
    normalized_expression = normalize_arcsinh(training_data)
    # Select subset of genes to use based on standard deviation
    std_dev = np.std( normalized_expression, axis = 0)
    sorted_idx = np.argsort(std_dev)
    # Train on PCA componenets rather than gene expression
    pca = PCA(n_components = pct_var)
    training_pca_coordinates = pca.fit_transform( normalized_expression[:,sorted_idx[-num_genes:]] )
    # Fit logisitc regression model and cross validate
    LR = LogisticRegressionCV( penalty = "l2", Cs = np.logspace(-7,1,9), dual = False, class_weight = "balanced", multi_class = "ovr", solver = "sag")
    LR.fit( training_pca_coordinates, training_labels )
    
    # Save cross validation results
    if save_results:
        directory = "results/classification/training_results/"
        parameters = str(num_genes) + "_" + str(int(pct_var*100))
        for ct in np.unique( training_labels ):
            pd.DataFrame( LR.scores_[ct], columns=[f'{x:.20f}' for x in LR.Cs_]).to_csv(
                         directory + ct + "_" + parameters + ".csv")
        confusion_mat = pd.DataFrame(confusion_matrix( training_labels, LR.predict( training_pca_coordinates )),
                                     index = np.unique(training_labels), columns = np.unique(training_labels))
        confusion_mat.to_csv(directory + "confusion_matrix" + "_" + parameters + ".csv")
        
    def feature_selection(X):
        return(X[:,sorted_idx[-num_genes:]])
    # Pack functions for classification into dictionary
    return {'normalization': normalize_arcsinh,
            'feature_selection': feature_selection,
            'pca_function': pca.transform, 
            'prediction_function': LR.predict,
            'prediction_probability': LR.predict_proba}
    #    return(LR)
#%%
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
    parser = ArgumentParser('Train classifier using given parameters')
    parser.add_argument("training_data", help = "sparse matrix (.npz containing training data)",
                        type=str,default="")
    parser.add_argument("training_labels", help = ".csv file containing training labels (cell types)",
                        type=str,default="")
    #parser.add_argument("gene_names", help = "file containing gene names",
    #                    type=str,default="")
    parser.add_argument("-k","--k_fold",
                        help = "Number of partitions to use for k-fold cross validation"
                        ,type=int,default=5)
    parser.add_argument("-n","--num_genes",
                        help = "Only use the n highest genes by standard deviation"
                        ,type=int ,default=1024)
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
    Predict cell-types of each cell in the single cell expression data using
    pre-defined marker genes
    
    Parameters
    ----------
    see _argparse
    
    Returns
    -------
    No return value. 
    train_LR function will save cross-validation results if specified. 
    For each cell type in the training labels, train_LR will write accuracy 
    scores over a range of regularization parameters. Will also save a confusion matrix
    """
    # Parse inputs
    parser = _argparse()
    argp = parser.parse_args(args[1:])
    # Read data
    sparse_training = sp.sparse.load_npz( argp.training_data )
    training_counts = np.array( sparse_training.todense() )
    training_labels = np.squeeze(np.array(pd.read_csv("results/classification/training_labels.csv",index_col = 0 )))
    
    train_LR( training_counts, training_labels, num_genes = argp.num_genes, pct_var = argp.pca_variance/100 )


if __name__ == '__main__':
    from sys import argv
    exit(main(argv))
    