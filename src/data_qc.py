import os
import pandas as pd
import scipy as sp
from scipy import sparse
import numpy as np
from sklearn.decomposition import PCA
from argparse import ArgumentParser, FileType
from sklearn.cluster import DBSCAN



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
    parser.add_argument("input_file", help = "input sparse expression matrix",
                        type=str,default="")
    parser.add_argument("gene_names", help = "gene names corresponding to rows of expression matrix",
                        type=str,default="")
    parser.add_argument("cell_labels", help = "meta data for each cell",
                        type=str,default="")    
    parser.add_argument("tsne_file", help = "tsne_coordinates of all cells",
                        type=str,default="") 
    parser.add_argument("mitochondrial_genes", help = "list of mitochondrial genes",
                        type=str,default="")
    parser.add_argument("-m","--mitochondrial_threshold", help = "output sparse expression matrix",
                        type=float,default=3)
    parser.add_argument("-n","--num_genes", help = "output sparse expression matrix",
                        type=int,default=750)
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
    scRNA_dense_expression = np.array(scRNA_expression.todense())

    gene_symbols = np.squeeze(np.array(pd.read_csv(argp.gene_names, header=None)))
    mito_genes = pd.read_table(argp.mitochondrial_genes)
    cell_labels = pd.read_csv(argp.cell_labels,header=None,index_col = 0)
    sequencing_metrics -= pd.read_csv( argp.metric_file, index_col = 0)
    
    
    
    
    tsne_coordinates = np.array(pd.read_csv(argp.tsne_file,header=None)


    # Find all mitochondrial genes in expression data
    converted_symbols = [ x.replace(".","-") for x in gene_symbols.index.values if "." in x]
    mt_gene_names = [ gene for gene in mito_genes.Symbol if gene in converted_symbols ]
    mt_gene_names = [ x.replace("-",".") for x in mt_gene_names ]
    mt_gene_idx = [ gene in mt_gene_names for gene in gene_symbols.index.values]
    mito_expression = np.array(scRNA_expression.todense()[:,mt_gene_idx])

    # Use only genes where greater than 10% of cells express the gene
    idx = np.squeeze(np.array(np.sum(mito_expression > 1, axis = 0 ) > mito_expression.shape[0] *.1))
    
    mean_expression = np.mean( mito_expression[:,idx],1)
    # Run pca using only mitochondrial genes
    pca = PCA()
    pca_transform = pca.fit_transform(mito_expression[:,idx]) 

    ## DBSCAN on t-SNE coordiantes ##
    tsne_coordinates = np.array(pd.read_csv("results/tSNE/filtered_tsne_coordinates.csv",header=None))
    density_cluster = DBSCAN(eps = 2, min_samples = 50)
    cluster_id = density_cluster.fit_predict(tsne_coordinates)
    
    
    #plt.figure(num=None,figsize=(8,8))
    #plt.scatter( tsne_coordinates[:,0], tsne_coordinates[:,1],s=2,c=cluster_id)
    #plt.axis("square")
    #plt.show()

    # Create dataframe 
    df = pd.DataFrame( np.concatenate( ( mito_expression,
                                        mean_expression[:,np.newaxis],
                                        pca_transform[:,0:2],
                                        num_genes_detected[:,np.newaxis]),axis=1),
                        columns = [gene for gene in gene_names.index.values[mt_gene_idx]]+["average","PC1","PC2","Number of genes detected"])
    df.to_csv("results/qc_info.csv",index=False)
    
    # Remove cells with high mitochondrial expression or low number of detected genes
    #idx = np.logical_and(mean_expression < argp.mitochondrial_threshold , num_genes_detected > argp.num_genes)    
    idx = mean_expression < argp.mitochondrial_threshold
    
    filtered_expression = scRNA_dense_expression[idx,:]
    filtered_labels = cell_labels.loc[idx]
    filtered_tsne_coordinates = tsne_coordinates.loc[idx,:]
    filtered_sparse_expression = sp.sparse.coo_matrix(filtered_expression)
    
    filtered_tsne_coordinates.to_csv("data/filtered/filtered_tsne_coordinates.csv",index=False,header=False)
    filtered_labels.to_csv("data/filtered/filtered_labels.csv",index=False,header=False)
    sp.sparse.save_npz("data/filtered/filtered_sparse.npz", filtered_sparse_expression)
        
if __name__ == '__main__':
    from sys import argv
    exit(main(argv))