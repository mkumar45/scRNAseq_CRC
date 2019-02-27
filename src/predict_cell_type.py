from argparse import ArgumentParser
import matplotlib
from matplotlib import pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42 # Editable text when exporting pdfs
import numpy as np
import pandas as pd
import scipy as sp
import seaborn as sns
import train_classifier #import train_LR

from sklearn.metrics import confusion_matrix


def normalize_arcsinh(count_matrix,cofactor = 1000):
    total_counts = np.squeeze(np.sum(count_matrix,axis = 1))
    normalized_counts = np.arcsinh(count_matrix / total_counts[:,np.newaxis] * cofactor)
    return(normalized_counts)
    
#%% Use trained classifier to predict cell-types
def predict_cell_types(X, trained_classifier):
    """
    Description
    -----------   
    Predict cell-types using trained classifier
    Parameters
    ----------    
    X : numpy array
        count data to use as input for prediction
    trained_classifier : dict
        output of train_LR fuction. Dictionary containing handles to functions 
        used for classification
    Returns
    -------
    predicted_cell_types : numpy array
        predicted cell type label
    prediction_probabilities : numpy array
        probability each cell belongs to each class
    """
    # Unpack functions for prediction
    normalize = trained_classifier["normalization"]
    select_features = trained_classifier["feature_selection"]
    pca_transform = trained_classifier["pca_function"]
    predict_function = trained_classifier["prediction_function"]
    predict_prob = trained_classifier["prediction_probability"]
    # Apply functions
    predictors = pca_transform( select_features(normalize(X) ) )
    predicted_cell_types = predict_function( predictors )
    prediction_probabilities = predict_prob( predictors )
    
    return predicted_cell_types, prediction_probabilities
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
    parser.add_argument("test_data", help = "sparse matrix (.npz containing data to classify)",
                        type=str,default="")
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

def main(args):
    """
    Description
    -----------   
    Predict cell-type of each cell in the single cell expression data using
    pre-defined marker genes
    
    Parameters
    ----------
    see argparse_
    Returns
    -------
    Saves two outputs
        1) predicted_cell_types ( in results/classification/)
        2) prediction_probabilities ( in results/classification/)
    """
    # Parse inputs
    parser = _argparse()
    argp = parser.parse_args(args[1:])
    # Read training data
    sparse_training = sp.sparse.load_npz( argp.training_data )
    training_data = np.array( sparse_training.todense() )
    training_labels = np.array(pd.read_csv("results/classification/training_labels.csv",index_col = 0 ))

    # Train classifier        
    classifier = train_classifier.train_LR( training_data, training_labels, num_genes = argp.num_genes, pct_var = argp.pca_variance/100, save_results = False  )

    training_predictions, training_probabilities = predict_cell_types( training_data, classifier )
    confusion_mat = pd.DataFrame(confusion_matrix( training_labels, training_predictions ),
                             index = np.unique(training_labels), columns = np.unique(training_labels))
    # Make predictions for all cells
    sparse_scRNA_expression = sp.sparse.load_npz( argp.test_data )
    scRNA_expression = np.array( sparse_scRNA_expression.todense() )
    predicted_cell_types, prediction_probabilities = predict_cell_types( scRNA_expression, classifier ) 
    
    sequencing_metrics = pd.read_csv("data/combined/processed_metrics.csv",index_col = 0)
    # Save outputs
    pd.DataFrame(prediction_probabilities, 
                 index = sequencing_metrics.index.values,
                 columns = np.unique(predicted_cell_types) ).to_csv("results/classification/prediction_probabilities.csv")
    pd.DataFrame(predicted_cell_types, 
                 index = sequencing_metrics.index.values,
                 columns = ["predicted_cell_type"] ).to_csv("results/classification/predicted_cell_types.csv")

if __name__ == '__main__':
    from sys import argv
    exit(main(argv))
    