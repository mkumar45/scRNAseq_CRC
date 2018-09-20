import numpy as np
import scipy as sp
from scipy import sparse
import pandas as pd
from argparse import ArgumentParser, FileType


#%%    
def product_mean(ligand_expression,receptor_expression):
    average_ligand = np.transpose(np.squeeze([np.mean(x,axis=1) for x in ligand_expression]))
    average_receptor = np.transpose(np.squeeze([np.mean(x,axis=1) for x in receptor_expression]))
    
    num_interactions = average_ligand.shape[0]
    num_cell_types = average_ligand.shape[1]
    interaction_scores = np.empty((num_interactions,num_cell_types**2))
        
    for n in range(num_interactions):
        score_mat = np.outer( average_receptor[n,:], np.transpose(average_ligand[n,:]) )
        interaction_scores[n,:] = np.reshape(score_mat,(num_cell_types**2,))
    return(interaction_scores)
#%%
def product_single_cell_type(ligand_expression,receptor_expression):
    average_ligand = np.transpose(np.squeeze([np.mean(x,axis=1) for x in ligand_expression]))
    average_receptor = np.transpose(np.squeeze([np.mean(x,axis=1) for x in receptor_expression]))
    
    num_interactions = average_ligand.shape[0]
    num_cell_types = average_ligand.shape[1]
    interaction_scores = np.empty((num_interactions,num_cell_types))
        
    for n in range(num_interactions):
        print(average_receptor[n,:])
        print( np.mean(np.transpose(average_ligand[n,:])))
        score_vec = np.dot( average_receptor[n,:], np.mean(np.transpose(average_ligand[n,:])) )
        print(score_vec)
        interaction_scores[n,:] = score_vec#np.reshape(score_mat,(num_cell_types**2,))
    return(interaction_scores)

#%% Score receptor-ligand interactions occuring between all cell types
def cell_cell_communication(scRNA_expression, gene_symbols, cell_type_labels, RL_list, group ):
    """
    Description
    -----------     
    Calculate receptor-ligand interactions scores between all cell-types 
    Parameters
    ----------
    

    scRNA_expression :
    gene_names :
    cell_type_labels :
    RL_list :
    Returns
    -------
    
    """ 
        
    # Use only interactions with both ligand and receptor measured
    ligand_idx = RL_list.Ligand_GeneSymbol.isin( gene_symbols )
    receptor_idx = RL_list.Receptor_GeneSymbol.isin( gene_symbols )
    interaction_idx = np.logical_and(ligand_idx, receptor_idx)
    RL_list = RL_list.loc[interaction_idx,:]
    RL_list.drop_duplicates()
      
    #cell_type_pairs = get_cell_type_pairs( cell_type_labels )
    

    # Computue scores
    if group is None: # for all cells 
        interaction_scores =  compute_scores(scRNA_expression, gene_symbols, cell_type_labels, RL_list )
    else: # compute scores for each unique group
        indexes = np.unique(group, return_index=True)[1]
        unique_groups = [np.array(group)[index] for index in sorted(indexes)] 
        
        interaction_scores_list = []
        for grp in range(len(unique_groups)):
            idx = np.squeeze( group.values == unique_groups[grp] )  
            expression_subset = scRNA_expression[ idx,:]                      
            interaction_scores = compute_scores( expression_subset, gene_symbols, cell_type_labels.loc[idx,:], RL_list  )
            #group_pairs = get_cell_type_pairs( cell_type_labels.loc[idx,:] ) # get cell-type pairs for current group
            interaction_scores_list.append(  pd.DataFrame( interaction_scores,
                                                           columns = RL_list.Ligand_GeneSymbol + "_" + RL_list.Receptor_GeneSymbol,
                                                           index = group.index.values[idx]) ) 
    return(interaction_scores_list)
            #missing_pairs = [ pair for pair in cell_type_pairs if pair not in group_pairs]
            #for pair in missing_pairs:
            #    df[pair]= 0
            #df = df[cell_type_pairs] #reorder
            #sort_df = sort_scores(interaction_scores, np.unique(cell_type_labels.loc[idx,:]), RL_list )

            # Write output files
            #df.to_csv(output_dir+unique_groups[grp]+"_interaction_scores.csv",sep= ",")
            #sort_df.to_csv( output_dir+unique_groups[grp]+"_sorted_scores.csv",sep= ",")
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

    cell_types = np.unique(cell_type_labels)
    num_cell_types = len(cell_types)
    cell_type_1 = np.repeat(np.arange(num_cell_types),num_cell_types)
    cell_type_2 = np.tile(np.arange(num_cell_types),num_cell_types)
    
    cell_type_pairs = cell_types[cell_type_1]+ "_" + cell_types[cell_type_2] 
    return( cell_type_pairs )

          
#%% 
def import_RL( filepath ):
    RL_list = pd.read_csv( filepath )
    RL_list.rename(columns={"Ligand_ApprovedSymbol":"Ligand_GeneSymbol",
                            "Receptor_ApprovedSymbol":"Receptor_GeneSymbol"},
                            inplace=True)
    RL_list[["Ligand_GeneSymbol","Receptor_GeneSymbol"]]
    #RL_list = RL_homologs(RL_list,input_species = "human")
    return(RL_list[["Ligand_GeneSymbol","Receptor_GeneSymbol"]])

#%% Wrapper function to calculate scores
def compute_scores( scRNA_expression, gene_symbols, cell_type_labels, RL_list ):
        
    #ligand_expression = cell_type_specific_expression(scRNA_expression, gene_symbols, RL_list.Ligand_GeneSymbol,cell_type_labels)
    #receptor_expression = cell_type_specific_expression(scRNA_expression, gene_symbols,RL_list.Receptor_GeneSymbol,cell_type_labels)
    
    
    num_cells = scRNA_expression.shape[0]
    num_interactions = RL_list.shape[0]
    ligand_expression = np.empty((num_cells,num_interactions))
    receptor_expression = np.empty((num_cells,num_interactions))

    for n in range(num_interactions):
        ligand_expression[:,n] = np.squeeze( scRNA_expression[:, gene_symbols == RL_list.Ligand_GeneSymbol.values[n] ])
        receptor_expression[:,n] = np.squeeze(scRNA_expression[:, gene_symbols == RL_list.Receptor_GeneSymbol.values[n] ] )
    
    interaction_scores = receptor_expression* np.mean( ligand_expression, axis = 0)
    #interaction_scores = pd.DataFrame(interaction_scores, index = cell_type_labels.index.values)
    #interaction_scores.columns = RL_list.Ligand_GeneSymbol + "_" + RL_list.Receptor_GeneSymbol
    #interaction_scores = score_function(ligand_expression,receptor_expression)
    #null_scores = self.calculate_null_scores(RL_list,score_function,num_null_reps)
    return(interaction_scores)

#%% Create list with cell-type specific expression for each gene in gene_symbols
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
        
        gene_idx[ return_genes == sym ] = np.nonzero( gene_symbols == sym )[0]
        
    ct_expression = [[] for n in range(num_cell_types)]

    for ct in range(num_cell_types):
        ct_idx = cell_type_labels.values == cell_types[ct]
        
        ct_expression[ct] = scRNA_expression[ np.squeeze(ct_idx), gene_idx.astype(int)] # cast to int for indexing
    return(ct_expression)
    
   
#%% Write sorted interaction scores to file
def sort_scores(interaction_scores,cell_type_labels,RL_list):
    cell_types = np.unique(cell_type_labels)
    num_cell_types = len(cell_types)
    
    
    
    vec = interaction_scores.reshape(interaction_scores.size,1)
    arg_sort = np.argsort(vec,axis=0)
    sort_idx = np.unravel_index(arg_sort,interaction_scores.shape)
    
    row_idx = np.squeeze(sort_idx[0])
    col_idx = np.squeeze(sort_idx[1])
    #row_idx, col_idx = sort_indices(interaction_scores)

    
    #if score_function == hf.product_mean:
    df = pd.DataFrame( {"Receptor": RL_list.iloc[row_idx].Receptor_GeneSymbol,
                        "Ligand": RL_list.iloc[row_idx].Ligand_GeneSymbol,
                        "Receptor_cell_type": cell_types[col_idx // num_cell_types], # integer divsion
                        "Ligand _cell_type": cell_types[col_idx % num_cell_types], # modulus
                        "Score": interaction_scores[row_idx,col_idx]}) 
    #elif score_function == hf.product_single_cell_type: 
    #    df = pd.DataFrame( {"Receptor": RL_list.iloc[row_idx].Receptor_GeneSymbol,
    #            "Ligand": RL_list.iloc[row_idx].Ligand_GeneSymbol,
    #            "Receptor_cell_type": cell_types[col_idx],
    #            "Score": interaction_scores[row_idx,col_idx]}) 
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
    parser.add_argument("gene_names", help = "gene names corresponding to rows of expression matrix",
                        type=str,default="")
    parser.add_argument("cell_type_labels", help = "cell-type for each cell in dataset",
                        type=str,default="")
    parser.add_argument("RL_file",
                        help = "file containing receptor ligand list",
                        type=str,default="")
    parser.add_argument("output_file", help = "directory to save output files",
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
    # Read data
    scRNA_expression = np.array( sp.sparse.load_npz( argp.input_file ).todense() )
    gene_symbols = np.squeeze(np.array(pd.read_csv(argp.gene_names,header=None)))
    RL_list = import_RL( argp.RL_file )
    cell_type_labels = pd.read_csv(argp.cell_type_labels,header=None)
    #scRNA_expression = np.array(sp.sparse.load_npz( "data/combined/combined_sparse.npz" ).todense())
    #gene_symbols = np.squeeze( np.array( pd.read_csv("data/combined/combined_gene_names.csv",header=None) ))
    #RL_list = import_RL( "data/mouse_receptor_ligand.csv" )
    #cell_type_labels = pd.read_csv("results/classification/predicted_cell_types.csv",header=None)
    #group = pd.read_csv("data/combined/combined_labels.csv",index_col = 0,header=None)
    
    if argp.group != None:
        group = pd.read_csv( argp.group, index_col = 0, header = None )
    else:
        group = None
    
    interaction_scores = cell_cell_communication(scRNA_expression, gene_symbols, cell_type_labels, RL_list, group )
    interaction_scores = pd.concat(interaction_scores)
    interaction_scores.to_csv(argp.output_file,sep= ",")
    

if __name__ == '__main__':
    from sys import argv
    exit(main(argv))