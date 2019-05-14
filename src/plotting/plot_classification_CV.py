from argparse import ArgumentParser
import glob
import matplotlib
from matplotlib import pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42 # Editable text when exporting pdfs
import numpy as np
import os
import pandas as pd
import scipy as sp
import seaborn as sns

cv_files =  glob.glob("results/classification/training_results/*.csv")
filenames = [ os.path.basename(file).split(".")[0] for file in cv_files]
ct = [ "_".join(x.split("_")[0:-2]) for x in filenames]
num_genes = [ int(x.split("_")[-2]) for x in filenames]
pct_var = [ int(x.split("_")[-1]) for x in filenames]

pca_parameters = np.unique(pct_var)
num_gene_parameters = np.unique(num_genes)
cell_types = np.unique(ct)


mean_cell_type_accuracy = {}
sem_cell_type_accuracy = {}
for ct in cell_types:
    mean_cell_type_accuracy[ct] = pd.DataFrame( [], index = num_gene_parameters, columns= pca_parameters )
    sem_cell_type_accuracy[ct] = pd.DataFrame( [], index = num_gene_parameters, columns= pca_parameters )    
    
    
#%%
for file in cv_files:
    
    filename = os.path.basename(file).split(".")[0]
    ct = "_".join(filename.split("_")[0:-2])
    num_genes = int(filename.split("_")[-2])
    pct_var = int(filename.split("_")[-1]) 

    accuracy = pd.read_csv( file, index_col = 0)    
    
    mean_accuracy = np.mean( accuracy, axis = 0 )
    sem_accuracy = np.std( accuracy, axis = 0 ) / np.sqrt( accuracy.shape[0] ) 
    regularization_parameter = np.array([ float(x) for x in accuracy.columns.values])    
    
    mean_cell_type_accuracy[ct].loc[num_genes,pct_var] = np.max( mean_accuracy) 
    sem_cell_type_accuracy[ct].loc[num_genes,pct_var] = sem_accuracy[np.argmax(mean_accuracy)]
    
    #%% 
    plt.figure(figsize=(4,4))
    plt.errorbar( x = regularization_parameter,
               y =  mean_accuracy,
               yerr = sem_accuracy,
               c = 'k', linewidth = 2 )
    plt.xscale('log')
    plt.xlabel('Regularization parameter',fontsize=16)
    plt.ylabel('Accuracy',fontsize=16)
    plt.yticks(np.linspace(0,1,6))
#%%    
for ct in cell_types:
    mean_df = mean_cell_type_accuracy[ct]
    sem_df = sem_cell_type_accuracy[ct]
    pca_variances = mean_df.columns.values
    num_genes= mean_df.index.values
    
    plt.figure(figsize=(5,5))
    color_map = plt.get_cmap("Greys")
    color_key = color_key = dict(zip( pca_variances, color_map(np.linspace(0.25,1,len(pca_variances))) ) )
    for pct in pca_variances:
        plt.errorbar(num_genes, mean_df.loc[:,pct].values,
                     yerr = sem_df.loc[:,pct].values,
                 c = color_key[pct], linewidth = 2)
    plt.xscale('log')
    plt.xticks(num_genes)
    plt.xlabel('Number of genes',fontsize=16)
    plt.ylabel('Accuracy',fontsize=16)
    plt.yticks(np.linspace(0,1,6))
    plt.legend(pca_variances)
    plt.title(ct)
        
    plt.figure(figsize=(5,5))
    color_key = color_key = dict(zip( num_genes, color_map(np.linspace(0.25,1,len(num_genes))) ) )
    for num in num_genes:
        plt.errorbar(pca_variances, mean_df.loc[num,:].values,
                 yerr = sem_df.loc[num,:].values,
                 c = color_key[num], linewidth = 2)
    plt.xticks(pca_variances)
    plt.xlabel('Percent Variance',fontsize=16)
    plt.ylabel('Accuracy',fontsize=16)
    plt.yticks(np.linspace(0,1,6))
    plt.legend(num_genes)
    plt.title(ct)
    
    
#%% 
confusion_matrices = glob.glob("results/classification/training_results/confusion_matrices/*.csv")
score_df = pd.DataFrame( [], index = num_gene_parameters, columns= pca_parameters )

for file in confusion_matrices:

    filename = os.path.basename(file).split(".")[0]
    num_genes = int(filename.split("_")[-2])
    pct_var = int(filename.split("_")[-1]) 
    
    cMat = pd.read_csv(file,index_col = 0 )
    num_correct = np.sum(np.diagonal( cMat ))
    num_cells = np.sum( np.sum( cMat ) )
    score_df.loc[num_genes,pct_var] = num_correct/num_cells
    
plt.figure(figsize=(5,5))
color_map = plt.get_cmap("Greys")
color_key = color_key = dict(zip( pca_parameters, color_map(np.linspace(0.25,1,len(pca_parameters))) ) )
for pct in pca_parameters:
    plt.plot(num_gene_parameters, score_df.loc[:,pct].values,
             c = color_key[pct], linewidth = 5)
plt.xscale('log')
plt.xticks(num_gene_parameters,fontsize=16)
plt.xlabel('Number of genes',fontsize=24)
plt.ylabel('Accuracy',fontsize=24)
plt.yticks(np.linspace(0,1,6),fontsize=16)
leg = plt.legend(pca_parameters,fontsize=20, title = "PCA % Variance")

plt.savefig("figures/classification/accuracy_v_genes.pdf")
plt.figure(figsize=(5,5))
color_key = color_key = dict(zip( num_gene_parameters, color_map(np.linspace(0.25,1,len(num_gene_parameters))) ) )
for num in num_gene_parameters:
    plt.plot(pca_parameters, score_df.loc[num,:].values,
             c = color_key[num], linewidth = 5)
plt.xticks(pca_parameters,fontsize=16)
plt.xlabel('Percent Variance',fontsize=24)
plt.ylabel('Accuracy',fontsize=24)
plt.yticks(np.linspace(0,1,6),fontsize=16)
plt.legend(num_gene_parameters, title = "Number of genes",fontsize=20,)
plt.savefig("figures/classification/accuracy_v_variance.pdf")


