run_Dorothea = function( expression_data, regulon_object )
{
  
  
  viper_regulon = readRDS("data/TFregulons/Robjects_VIPERformat/mouse_dorothea2_regulon_v1.rds")
  sparse_matrix = readRDS("data/combined/combined_sparse.rds")
  
  gene_symbols = read.table("data/combined/combined_gene_names.csv")
  colnames(sparse_matrix) = gene_symbols$V1
  
  

  # Limit to only genes in viper regulons
  viper_genes = list()
  for (r in 1:length(viper_regulon))
  {
    viper_genes = union(viper_genes,names(viper_regulon[[r]]$tfmode))

  }
  viper_genes = unlist(viper_genes)
  gene_idx = gene_symbols$V1 %in% viper_genes
  
  
  # Use only high confidnece TFs
  split_names = strsplit( names(viper_regulon),"_")
  confidence_level = sapply(split_names,`[`,2)
  TF_names = sapply(split_names,'[',1)
  TF_idx = confidence_level %in% c("A","B","C")
  
  # Estimate TF activities
  TF_activities = viper(eset = t(scale(as.matrix(sparse_matrix)[,gene_idx],center=TRUE,scale=TRUE)),
                        regulon = viper_regulon[TF_idx], nes = T, method = 'none', minsize = 4, eset.filter = F)
  
  
  # Save results
  write.csv(t(TF_activities), file = 'results/TF_activities_confidence_ABC.csv', col.names = TRUE, row.names = FALSE)
  
  
}