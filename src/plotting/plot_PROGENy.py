from argparse import ArgumentParser
import matplotlib
from matplotlib import pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42 # Editable text when exporting pdfs
import numpy as np
import pandas as pd
import scipy as sp
import seaborn as sns

#%% Normalize count matrix using arcsinh transformation. Assumes each row is a sample/cell and columns correspond to genes
def normalize_arcsinh(count_matrix,cofactor = 1000):
    total_counts = np.squeeze(np.sum(count_matrix,axis = 1))
    normalized_counts = np.arcsinh(count_matrix / total_counts[:,np.newaxis] * cofactor)
    return(normalized_counts)
    

def plot_score_heatmap(pathway_scores, cell_type_labels, cmap = "coolwarm" ):
    
        sns.heatmap( pathway_scores,
                    center = 0, cmap = cmap)
        plt.figure()
        sns.heatmap( sp.stats.zscore(pathway_scores),
                     center = 0, cmap = cmap)
        
        cell_types = np.unique(cell_type_labels)
        num_cell_types = len(cell_types)

        pathway_average_scores = np.empty((num_cell_types,pathway_scores.shape[1],))
        
        for ct in range(num_cell_types):
            ct_idx = cell_type_labels.values == cell_types[ct]
            pathway_average_scores[ct,:] = np.mean( pathway_scores.loc[np.squeeze(ct_idx),:], axis=0)
                    
        
        plt.figure()
        sns.heatmap( pd.DataFrame(pathway_average_scores, columns = pathway_scores.columns, index = cell_types),
                     cmap = cmap, center = 0)
        plt.figure()
        sns.heatmap(pd.DataFrame(sp.stats.zscore(pathway_average_scores), columns = pathway_scores.columns, index = cell_types),
                    cmap = cmap, center = 0)
        

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
    
def compute_permuted_pathway_score(scRNA_expression, gene_symbols,  PROGENy_coefficients, num_permutations = 1000 ):

    normalized_expression = sp.stats.zscore( scRNA_expression )
    nan_idx = np.logical_not(np.all(np.isnan( normalized_expression ), axis = 0 ))

    normalized_expression = normalized_expression[:,nan_idx]
    gene_symbols = gene_symbols[nan_idx]

    # Pathway genes only 
    idx = [x in PROGENy_coefficients.index.values for x in gene_symbols]

    coeff_subset = PROGENy_coefficients.loc[ gene_symbols[idx] ] # only use coefficients found in expression
    coeff_subset = coeff_subset.reindex( gene_symbols[idx] ) # Reorder for multiplication
    
    permuted_scores = np.empty( (scRNA_expression.shape[0], coeff_subset.shape[1], num_permutations))
    
    for n in range(num_permutations):  
        expression_subset = normalized_expression[:,np.random.randint(0, normalized_expression.shape[1], np.sum(idx) )]
        permuted_scores[:,:,n] = np.matmul( expression_subset, coeff_subset)
        
        
    #df = pd.DataFrame(pathway_scores,columns=coeff_subset.columns)
    return(permuted_scores)

#%% Import data
PROGENy_scores = pd.read_csv("results/PROGENy/PROGENy_scores.csv",index_col=0)
predicted_cell_types = pd.read_csv("results/classification/processed/predicted_cell_types.csv",index_col = 0 )
count_matrix = np.array(sp.sparse.load_npz("results/classification/processed/processed_counts.npz").todense())
normalized_expression = normalize_arcsinh( count_matrix )
gene_symbols = np.squeeze(np.array(pd.read_csv("data/mouse_gene_symbols.csv",header=None)))
PROGENy_coefficients = pd.read_csv( "data/PROGENy_mouse_model_v2.csv", index_col = 0 )  

#%% Plot scores heatmap
plot_score_heatmap( PROGENy_scores, predicted_cell_types.predicted_cell_type, cmap = "coolwarm" )
plt.savefig("figures/PROGENy/cell_type_heatmap.pdf")

# Permuted scores
permuted_scores = compute_permuted_pathway_score(normalized_expression, gene_symbols, PROGENy_coefficients, num_permutations=1000 )

for n in range(PROGENy_scores.shape[1]):
    pathway_score = np.array(PROGENy_scores.iloc[:,n])
    plt.figure(num=None,figsize=(8,8))
    plt.scatter(PROGENy_scores.iloc[:,n],-np.log10(np.sum(np.abs(permuted_scores[:,n,:]) > np.abs(pathway_score[:,np.newaxis]),axis=1)/1000),
                s = 2, c = 'k' )#np.squeeze( predicted_cell_type == "CAF") )
    #plt.gca().set_yscale('log')
    plt.xlabel(PROGENy_scores.columns[n] + " PROGENy Score")
    plt.ylabel('-log10(p-value)')
    #plt.savefig( "figures/PROGENy/" + pathway_scores.columns[n] +"_volcano.pdf")
    plt.show()
    
    
#%% Plot number of pathway genes detected per cell 
for p in range(PROGENy_coefficients.shape[1]):
    
    gene_idx = PROGENy_coefficients.iloc[:,p] != 0
    pathway_genes = PROGENy_coefficients.index.values[gene_idx]
    
    idx = [x in pathway_genes for x in gene_symbols]
    expression_subset = normalized_expression[:,idx]
    plt.title("Number of genes in " + PROGENy_coefficients.columns[p] + " pathway: " + str(expression_subset.shape[1]))
    num_pathway_genes = np.sum( expression_subset != 0 , axis = 1)
    plt.hist(num_pathway_genes,bins=expression_subset.shape[1],color = 'k')
    plt.xlabel("Number of pathway genes detected")
    plt.ylabel("Number of cells")
    #plt.savefig("figures/PROGENy/" + PROGENy_coefficients.columns[p] + "_detected_genes.pdf")
    plt.show()
    
    
#%% Plot PROGENy gene coefficient vs. number of cells gene is detected 
for p in range(PROGENy_coefficients.shape[1]):
    
    
    gene_idx = PROGENy_coefficients.iloc[:,p] != 0
    pathway_genes = PROGENy_coefficients.index.values[gene_idx]
    
    idx = [x in pathway_genes for x in gene_symbols]
    expression_subset = normalized_expression[:,idx]
    
    print("Number of genes in " + PROGENy_coefficients.columns[p] + " pathway: " + str(expression_subset.shape[1]))
    num_cells_detected = np.sum( expression_subset != 0 , axis = 0)
    
    
    idx = [ gene in gene_symbols[idx] for gene in PROGENy_coefficients.index.values]
    plt.scatter(num_cells_detected, PROGENy_coefficients.loc[idx,PROGENy_coefficients.columns[p]], c = 'k')
    plt.xlabel( "Number of cells with gene expression")
    plt.ylabel( "PROGENy coefficient")
    #plt.savefig("figures/PROGENy/" + PROGENy_coefficients.columns[p] + "_gene_coefficients.pdf")
    plt.show()
    
    
#%% Plot number of pathway genes detected vs. PROGENy score
for p in range(PROGENy_coefficients.shape[1]):    
    gene_idx = PROGENy_coefficients.iloc[:,p] != 0
    pathway_genes = PROGENy_coefficients.index.values[gene_idx]
    
    idx = [x in pathway_genes for x in gene_symbols]
    expression_subset = normalized_expression[:,idx]
    
    print("Number of genes in " + PROGENy_coefficients.columns[p] + " pathway: " + str(expression_subset.shape[1]))
    num_pathway_genes = np.sum( expression_subset != 0 , axis = 1)
    plt.scatter(num_pathway_genes,PROGENy_scores.loc[:,PROGENy_coefficients.columns[p]],s=2, c = 'k')
    plt.xlabel("Number of measured pathway genes")
    plt.ylabel("PROGENy Score")
    #plt.savefig("figures/PROGENy/" + PROGENy_coefficients.columns[p] + "_detected_genes_v_score.pdf")
    plt.show()