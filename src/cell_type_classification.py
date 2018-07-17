import numpy as np
import scipy as sp
import pandas as pd
import os 

from argparse import ArgumentParser, FileType
from sys import stderr, stdin, stdout
from sklearn import mixture
from sklearn.model_selection import KFold, GridSearchCV
from sklearn.decomposition import PCA
from sklearn.feature_selection import SelectKBest
from sklearn import svm 
from sklearn.pipeline import Pipeline

#%% Use trained classifier to predict cell-types
def predict_cell_types(self, trained_classifier_dict):
    """
    Description
    -----------   
    Parameters
    ----------    
    Returns
    -------
    """
    # Unpack functions for prediction
    pca_transform = trained_classifier_dict["pca_function"]
    n_comp = trained_classifier_dict["pca_num_comp"]
    predict_function = trained_classifier_dict["prediction_function"]
    predict_prob = trained_classifier_dict["prediction_probability"]

    predictors = pca_transform( np.transpose(self.data))[:,:n_comp]
    predicted_cell_types = predict_function( predictors )
    prediction_probabilities = predict_prob( predictors )
    
    
    # Add column data to object
    self.column_data["predicted_cell_type"] = predicted_cell_types
    self.column_data["prediction_probability"] = [x[0] for x in np.vsplit(prediction_probabilities, prediction_probabilities.shape[0])]
    #np.vsplit(prediction_probabilities, prediction_probabilities.shape[0])
    
#%% Train a classifer using training_object generated by create_training_data
def train_classifier(training_data, training_labels ):
    """
    Description
    -----------   
    Train a classifer using training_object generated by create_training_data
    Perform PCA decomposition for dimensionality reducation
    Parameters
    ----------    
    training_data: training expression data in COO sparse format
        subset of input data to use for training classifier
    training_labels : numpy array
        assigned cell-type based on marker expression to use for training classifier
    
    Returns
    -------
    estimator : result of using grid search
    """

    #def std_dev_func (X,y):
    #    return(np.std(X,axis=0))
    #std_dev_filter = SelectKBest( std_dev_func, num_genes )
    #training_subset = training_data.iloc[sort_idx[0:num_genes],:]
    
    training_data = training_data.todense()
    
    # Perform PCA
    pca = PCA()
    svm_classifier = svm.SVC(probability = True)#
    pipe = Pipeline(steps=[('pca', pca), ('svm', svm_classifier)])

    n_components = 2**np.linspace(4,8,5)
    Cs = np.logspace(-2, 2, 5)

    estimator = GridSearchCV(pipe,
                         dict(pca__n_components=n_components.astype(int),
                              svm__C=Cs), cv = 5,verbose=10 )
    estimator.fit( training_data, training_labels )
    return(estimator)
    
    
    #predictors = pca.fit_transform( training_subset.transpose() )
    #n_comp = np.nonzero( np.cumsum(pca.explained_variance_ratio_) > pca_var_to_keep)[0][0]
    #print(n_comp)
    #print("\n")
    #predictors = predictors[:,:n_comp]
    # Fit Decision Tree classifier using pca scores as predictors
    #decision_tree.fit( predictors, training_labels)
    #print( cross_val_score( decision_tree, predictors, training_labels, cv=10))
    # Pack functions for classification into dictionary
    #return {'pca_function': pca.transform, 
    #        'pca_num_comp': n_comp,
    #        'prediction_function': decision_tree.predict,
    #        'prediction_probability': decision_tree.predict_proba}
    
#%% Create subset of data with only cell-type marker expression data
def create_training_data( scRNA_expression, gene_names, marker_data, GMMs ):
    """
    Description
    -----------   
    Parameters
    ----------
    scRNA_expression : sparse matrix in coo format
        Dataframe containing input scRNA expression data. Each row is a unique gene
        and each column is a single cell
    gene_names :  
    
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
    marker_idx = [gene in marker_data.columns.values for gene in gene_names.index.values]
    marker_expression = scRNA_expression.tocsc()[:,np.squeeze(np.nonzero(marker_idx))]
    marker_expression = np.array(marker_expression.todense())
    
    num_cell_types = marker_data.shape[0]
    num_cells = scRNA_expression.shape[0]
    is_cell_type = np.zeros((num_cell_types,num_cells),dtype=bool)
    
    # Assign cluster 
    cluster_idx = np.array([ GMMs[n].predict( marker_expression[:,n].reshape(-1,1) ) for n in range(marker_expression.shape[1]) ])
    
    max_idx = np.array([ np.argmax(GMMs[x].means_)  for x in range(len(GMMs))])
    min_idx = np.array([ np.argmin(GMMs[x].means_)  for x in range(len(GMMs))])
    
    # Loop through 
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
    training_data = scRNA_expression.tocsc()[training_idx,:]
    # Find index of unique cell type and add cell type name to column data
    cell_type_idx = np.where( np.transpose(is_cell_type[:,training_idx]) )
    training_labels = marker_data.index.values[cell_type_idx[1]]
    
    return(training_data.tocoo(),training_labels,training_idx)
       
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
    marker_expression : sparse matrix in csc format
        Row for each marker gene and columns represent single cells
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
    Parameters
    ----------    
    Returns
    -------
    """
    parser = ArgumentParser('Classify cell-types using scRNA expression and predefined markers')
    parser.add_argument("input_file", help = "preprocessed input file to classify",
                        type=str,default="")
    parser.add_argument("output_dir", help = "directory to save output files",
                        type=str,default="")
    parser.add_argument("marker_file",
                        help = "file containing marker genes for classification",
                        type=str,default="")
    parser.add_argument("-k","--k_fold",
                        help = "Number of partitions to use for k-fold cross validation"
                        ,type=int,default=5)
    parser.add_argument("-n","--n_components",
                        help = "Number of mixture componenets to fit for Gaussian Mixture models"
                        ,type=int,default=5)
    parser.add_argument("-n","--num_genes",
                        help = "Only use the n highest genes by standard deviation"
                        ,type=int ,default=500)
    parser.add_argument("-pca","--pca_variance",
                        help = "% of PCA variance to keep when defining training inputs"
                        ,type=float  ,default=0.90)

    return parser
        

#%% Compute expression thresholds for each marker using otsu's method
#   Assumes two groups (high expressing and low expressing)
def classify_cell_types(args):
    """
    Description
    -----------   
    Predict cell-types of each cell in the single cell expression data using
    pre-defined marker genes
    
    Parameters
    ----------
    scRNA_expression : pandas dataframes
    genes in rows, single cells in columns.
    
    gene_markers : pandas dataframe
    
    
    Returns
    -------

    """
    # Parse inputs
    parser = _argparse()
    argp = parser.parse_args(args[1:])
    # Read data
    scRNA_expression = sp.sparse.load_npz( argp.input_file )
    gene_names = pd.read_csv("data/combined_gene_names.csv",index_col = 0,header=None)
    scRNA_expression = sp.sparse.load_npz( "data/combined_sparse.npz" )
    marker_df = pd.read_csv( argp.marker_file, index_col = 0 )
    
    # Use only markers are present in expression data
    idx = [ elem in gene_names.index.values for elem in marker_df.columns.values ]
    marker_subset = marker_df.loc[:,idx] 
    marker_subset = marker_subset.dropna(axis = 0, how = "all") # drop cell type if no markers remain
    marker_subset = marker_subset.reindex_axis(sorted(marker_subset.columns), axis=1)

    # Gaussian Mixture models 
    idx = [ gene in marker_subset.columns.values for gene in gene_names.index.values]
    csc_mat = scRNA_expression.tocsc()
    marker_expression = csc_mat[:,np.squeeze(np.nonzero(idx))].todense()
    marker_expression = np.array(marker_expression)
    
    #fit_GMM(marker_expression, argp.num_components, argp.k_fold)
    GMMs, GMM_dict, BIC_dict = fit_GMM(marker_expression, 5, 5)
    
    
    training_data, training_labels, training_idx = create_training_data(scRNA_expression, gene_names, marker_subset, GMMs )
    
    trained_classifier_dict = train_classifier( training_data, training_labels )   
    estimator = predict_cell_types(trained_classifier_dict)
    
    
    predicted_cell_type = np.array(estimator.best_estimator_.predict( scRNA_expression.todense() ))
    prediction_probabilities = np.array(estimator.best_estimator_.predict_proba( scRNA_expression.todense() ))
    
    np.savetxt(argp.output_dir + "predicted_cell_types.csv", predicted_cell_type, delimiter = ",",fmt="%s")
    np.savetxt(argp.output_dir + "prediction_probabilities.csv", prediction_probabilities, delimiter = ",",fmt = "%10.5f", header = ",".join(map(str, np.unique(training_labels))))
    
    return(estimator.best_estimator_)
    
    #df = pd.DataFrame(columns=['col'])
    #for i in range(10):
        #df.loc[i,'col'] = np.zeros(3)



if __name__ == '__main__':
    from sys import argv
    exit(classify_cell_types(argv))
    