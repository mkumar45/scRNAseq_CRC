require(viper) 
run_Dorothea = function( expression_data, regulon_object, gene_symbols )
{
  
  sparse_matrix = readRDS(expression_data)
  viper_regulon = readRDS(regulon_object)
  predicted_cell_type = read.csv("results/classification/processed/predicted_cell_types.csv",row.names = 1)
  gene_symbols = read.table(gene_symbols)
  colnames(sparse_matrix) = gene_symbols$V1
  
  # Limit to only genes in viper regulons
  #viper_genes = list()
  #for (r in 1:length(viper_regulon))
  #{
   # viper_genes = union(viper_genes,names(viper_regulon[[r]]$tfmode))

  #}
  #viper_genes = unlist(viper_genes)
  #gene_idx = gene_symbols$V1 %in% viper_genes
  
  # Use only high confidnece TFs
  split_names = strsplit( names(viper_regulon),"_")
  confidence_level = sapply(split_names,`[`,2)
  TF_names = sapply(split_names,'[',1)
  TF_idx = confidence_level %in% c("A")
  
  # Estimate TF activities
  count_matrix = as.matrix(sparse_matrix)
  total_counts = rowSums(count_matrix)
  normalized_counts = asinh(count_matrix / total_counts * 1000)
  TF_activities = viper(eset = t(scale(normalized_counts,center=TRUE,scale=TRUE)),
                        regulon = viper_regulon[TF_idx], nes = T, method = 'rank', minsize = 4, eset.filter = F,mvws = 1)
  
  rownames(TF_activities) = TF_names[TF_idx]
  colnames(TF_activities) = rownames(predicted_cell_type)
  # Save results
  write.csv(t(TF_activities), file = 'results/Dorothea/TF_activities_confidence_A.csv', row.names = TRUE)
}

args = commandArgs(trailingOnly=TRUE)
run_Dorothea( args[1], args[2], args[3] )