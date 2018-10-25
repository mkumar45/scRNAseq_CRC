import numpy as np
import scipy as sp
from scipy import sparse
import pandas as pd
from argparse import ArgumentParser, FileType
from sys import exit

import matplotlib
from matplotlib import pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42 # Editable text when exporting pdfs
import seaborn as sns



from src import pcreode

#%%    
def product_mean():
    """
    Description
    -----------
    Create pairwise cell-type labels     
    Parameters
    ----------       
    cell_type_labels:
    cell-type labels for each row of scRNA_expression matrix
    Returns
    -------
    cell_type_pairs: array
    array of cell-type pairs
    """
    return()
#%%
def product_single_cell_type(ligand_expression,receptor_expression):

    return(interaction_scores)

#%% Score receptor-ligand interactions occuring between all cell types
def cell_cell_communication(scRNA_expression, gene_symbols, cell_type_labels, RL_list, group ):
   
    return(interaction_scores_list)

#%%
def get_cell_type_pairs( cell_type_labels ):
    """
    Description
    -----------
    Create pairwise cell-type labels     
    Parameters
    ----------       
    cell_type_labels:
    cell-type labels for each row of scRNA_expression matrix
    Returns
    -------
    cell_type_pairs: array
    array of cell-type pairs
    """
    return( cell_type_pairs )   
#%% Write sorted interaction scores to file
def sort_scores(interaction_scores,cell_type_labels,RL_list):
    return(df[::-1])

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
    parser.add_argument("gene_symbols", help = "gene symbols corresponding to columns of expression matrix",
                        type=str,default="")
    parser.add_argument("cell_type_labels", help = "cell type label for each row of expression matrix",
                        type=str,default="")
    parser.add_argument("sample_labels", help = "sample label for each cell (row) of expression matrix",
                        type=str,default="")
    parser.add_argument("-g","--group",
                        help = "use to run communication on different groups of cells"
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
    # Read data and convert to pandas dataframe
    scRNA_expression = np.array( sp.sparse.load_npz( argp.input_file ).todense() )
    gene_symbols = np.squeeze(np.array(pd.read_csv(argp.gene_symbols,header=None)))
    cell_type_labels = np.squeeze(np.array(pd.read_csv(argp.cell_type_labels,header=None)))
    sample_labels = np.squeeze(np.array(pd.read_csv(argp.sample_labels,index_col=0,header=None)))
    
    scRNA_expression = np.array( sp.sparse.load_npz( "data/combined/combined_sparse.npz" ).todense() )
    gene_symbols = np.squeeze(np.array(pd.read_csv("data/combined/combined_gene_names.csv",header=None)))
    cell_type_labels = np.squeeze(np.array(pd.read_csv("results/classification/predicted_cell_types.csv",header=None)))
    sample_labels = np.squeeze(np.array(pd.read_csv("data/combined/combined_labels.csv",index_col=0,header=None)))
    
    df = pd.DataFrame( data = scRNA_expression, index = None, columns = gene_symbols)
    
    # Subset cells to use for pCreode
    model_idx = [ x.split('_')[0] == "AOM" for x in sample_labels ]
    cell_types = np.unique(cell_type_labels)
    ct_idx = [ct in cell_types for ct in cell_type_labels]
    idx = np.squeeze(np.logical_and( model_idx, ct_idx))
    df_subset = df.loc[idx,:]
    
    #%% PCA    
    data_pca = pcreode.PCA( df_subset )
    data_pca.get_pca()
    data_pca.pca_plot_explained_var( xlim=(0,10))
    
    pca_test_data = data_pca.pca_set_components( 5)
    fig = plt.figure( figsize=(12,12))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    cc = 'k'
    ax1.scatter( pca_test_data[:,0], pca_test_data[:,1], alpha=0.5, s=25, c=cc)
    ax2.scatter( pca_test_data[:,2], pca_test_data[:,1], alpha=0.5, s=25, c=cc)
    ax3.scatter( pca_test_data[:,2], pca_test_data[:,3], alpha=0.5, s=25, c=cc)
    ax4.scatter( pca_test_data[:,4], pca_test_data[:,3], alpha=0.5, s=25, c=cc)
    ax1.set_xlabel("PC1", fontsize=15), ax1.set_ylabel("PC2", fontsize=15)
    ax2.set_xlabel("PC3", fontsize=15), ax2.set_ylabel("PC2", fontsize=15)
    ax3.set_xlabel("PC3", fontsize=15), ax3.set_ylabel("PC4", fontsize=15)
    ax4.set_xlabel("PC5", fontsize=15), ax4.set_ylabel("PC4", fontsize=15)
    
    plt.figure(figsize=(9,9))
    for ct in cell_types:
        ct_idx = cell_type_labels[idx] == ct
        plt.scatter(pca_test_data[ct_idx,0], pca_test_data[ct_idx,1], alpha=0.5, s=25)
    plt.legend(cell_types)
    pca_reduced_data = data_pca.pca_set_components(5)

    #%% Density calculation
    dens = pcreode.Density( pca_reduced_data)
    radius_guess = dens.nearest_neighbor_hist( )
    
    density = dens.get_density( radius = radius_guess)
    #density_4 = dens.get_density( radius=0.4)
    #density_15 = dens.get_density( radius=1.5)
    dens.density_hist( n_bins=100)
    
    fig = plt.figure( figsize=(16,8))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    cc = density
    ax1.scatter( pca_reduced_data[:,0], pca_reduced_data[:,1], alpha=0.5, s=25, c=cc)
    ax2.scatter( pca_reduced_data[:,2], pca_reduced_data[:,1], alpha=0.5, s=25, c=cc)
    ax1.set_xlabel("PC1", fontsize=15), ax1.set_ylabel("PC2", fontsize=15)
    ax2.set_xlabel("PC3", fontsize=15), ax2.set_ylabel("PC2", fontsize=15)
    
    
    #%% Down sample
    noise = 10.0
    target = 25.0
    
    downed, downed_ind = pcreode.Down_Sample( pca_reduced_data, density, noise, target)
    fig = plt.figure( figsize=(16,8))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    cc = density[downed_ind]
    ax1.scatter( downed[:,0], downed[:,1], alpha=0.5, s=25, c=cc)
    ax2.scatter( downed[:,2], downed[:,1], alpha=0.5, s=25, c=cc)
    ax1.set_xlabel("PC1", fontsize=15), ax1.set_ylabel("PC2", fontsize=15)
    ax2.set_xlabel("PC3", fontsize=15), ax2.set_ylabel("PC2", fontsize=15)
    
    
    
    fig = plt.figure( figsize=(16,8))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    cc = density
    ax1.scatter( pca_reduced_data[:,0], pca_reduced_data[:,1], alpha=0.5, s=25, c=cc)
    ax2.scatter( pca_reduced_data[:,2], pca_reduced_data[:,1], alpha=0.5, s=25, c=cc)
    ax1.set_xlabel("PC1", fontsize=15), ax1.set_ylabel("PC2", fontsize=15)
    ax2.set_xlabel("PC3", fontsize=15), ax2.set_ylabel("PC2", fontsize=15)
    #%% Obtain and score p-Creode graphs
    file_path = "figures/pCreode/"
    out_graph, out_ids = pcreode.pCreode( data=pca_reduced_data, density=density, noise=noise, 
                                          target=target, file_path=file_path, num_runs=10)
    
    # Score graphs
    graph_ranks = pcreode.pCreode_Scoring( data=pca_reduced_data, file_path=file_path, num_graphs=10)
    
    #%% Plotting
    gid = graph_ranks[0] # this will select the first graph ID in the ranking from above
    analysis = pcreode.Analysis( file_path=file_path, graph_id=gid, data=pca_reduced_data, density=density, noise=noise)
    seed = 1 
    
    ct_labels = cell_type_labels[idx]
    analysis.plot_save_qual_graph( seed=seed, overlay= ct_labels.astype('S'), file_out='cell_type_overlay')

    for gene_symbol in ["Cdh5","Pecam1",
                        "Cd3d","Cd3e","Cd3g",
                        "Cd68","Csf1r",
                        "Cd14","S100a9","S100a8",
                        "Pdgfra","Col5a1","Tnc",
                        "Vim","Acta2","Des","Pdgfrb"]:
        analysis.plot_save_graph( seed=seed, overlay=df_subset[gene_symbol], file_out=gene_symbol, upper_range=1.25)



#%%
if __name__ == '__main__':
    from sys import argv
    exit(main(argv))
