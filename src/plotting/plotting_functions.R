#############################################
###### Test for significatnly different #####
# interactions between two subsets of cells #  
#############################################

#arcsinh_normalize = function( X, cofactor )
#{
#  total_counts = rowSums(count_matrix)
#  normalized_counts = np.arcsinh(count_matrix / total_counts[:,np.newaxis] * cofactor)
#}

interaction_test = function( scores_1, scores_2,  output_file = NA )
{
  # Use only numeric columns
  numeric = unlist(lapply(scores_1,is.numeric))
  scores_1 = scores_1[,numeric]
  scores_2 = scores_2[,numeric]
  
  #
  num_tests = ncol((scores_1))
  p_value = matrix(NA, nrow = num_tests, ncol= 1)
  fold_change = matrix(NA, nrow = num_tests, ncol= 1)
  SNR = matrix(NA, nrow = num_tests, ncol= 1)
  for (n in 1:num_tests)
  {
    test = wilcox.test( scores_1[,n], scores_2[,n] )
    p_value[n] = test$p.value
    fold_change[n] = mean(scores_1[,n])/mean(scores_2[,n])
    SNR[n] = ( median(scores_1[,n]) - median(scores_2[,n]) ) /  ( sd(scores_1[,n]) + sd(scores_2[,n]) )

    
    #if (test$p.value < 0.01 & !is.nan(p_value[n]) & abs(log2(fold_change[n])) > 2 )
    #{
      #print(colnames(interaction_scores)[n])
    #}
  }
  
  df = data.frame( colnames(scores_1)[1:num_tests],fold_change, SNR,p_value)
  colnames(df) = c("name","fold_change","SNR","p_value")
  if (!is.na(output_file) )
  {
    df_subset = subset(df, p_value < 0.01 & abs(log2(fold_change)) > 2)
    
    RL_names = strsplit(as.character(df_subset$name),"_")
    df_subset$ligand =sapply(RL_names,`[`,1)
    df_subset$receptor =sapply(RL_names,`[`,2)
    
    write.csv(df_subset[order(abs(log2(df_subset$fold_change))),],output_file,row.names=FALSE)
  }
  return(df)
}
###########################################################
# Volcano plot of interactions  different b/w APC and AOM #  
###########################################################
volacano_plot = function( df, metric = "fold_change", FC = 2, q_val = .01,
                          xmin = NA, xmax = NA, ymin = NA, ymax = NA, add_text = TRUE )
{
  if ( metric == "fold_change")
  {
    df = subset(df, fold_change != Inf)
    plt = ggplot( df , aes( y = -log10(FDR), x = log2(fold_change) ) ) + geom_point() + 
      geom_point( data = subset(df, abs( log2(fold_change)) > FC & FDR < q_val),aes(color = "red")) +
      theme_bw() +
      theme( aspect.ratio = 1)
  } else
  {
    plt = ggplot( df , aes( y = -log10(FDR), x = SNR ) ) + geom_point() + 
      geom_point( data = subset(df, abs( SNR ) > FC & FDR < q_val),aes(color = "red")) +
      theme_bw() +
      theme( aspect.ratio = 1)
  }
  
  if (add_text)
  {
    plt = plt + geom_text( data = subset(df, abs( SNR ) > FC & FDR < q_val),aes(color = "red",label = name))
  }
    

  
  if ( !is.na(xmin))
  {
    print( max(df$fold_change) )
    plt = plt + scale_x_continuous( limits = c(xmin,xmax))
  }
  if ( !is.na(ymin))
  {
    print( max(df$fold_change) )
    plt = plt + scale_y_continuous( limits = c(ymin,ymax))
  }
  return(plt)
}



#########################################
# Plot receptor/ligand expression #
#########################################
plot_expression = function( expression_subset, gene_symbol, by = "sample_type"  )
{
  bplot = ggplot( expression_subset, 
                  aes_string( x = by, y = gene_symbol ),
                  color = by) + 
          geom_jitter(height = 0.01) + 
          theme_bw() + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                aspect.ratio = 1, text = element_text(size=20)) 
     # position = position_jitterdodge()+ 
    #scale_y_continuous( limits = c(-0.25,1), breaks = c(-0.05,0,0.25,0.5,1)) +
    #stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.5) +
  return(bplot)
  
}

plot_multi_gene_expression = function(count_matrix, gene_symbols, jitter = 0.05 )
{
  tall_df = data.frame()
  for (n in 1:length(gene_symbols))
  {
    df = count_matrix[,c(gene_symbols[n],"cell_type","sample_type","tissue")] 
    df$gene_symbol = gene_symbols[n]
    names(df) = c("counts","cell_type","sample_type","tissue","gene_symbol")
    
    tall_df = rbind(tall_df,df) 
  }
  #tall_df$normalized_counts = log10(tall_df$counts)
  plt = ggplot( tall_df, aes_string( x = "gene_symbol", y = "counts", color="sample_type")) + 
    geom_point(size=0.01,position = position_jitterdodge(jitter.height = jitter)) + 
    #scale_y_continuous( limits = c(-0.25,1), breaks = c(-0.05,0,0.25,0.5,1)) +
    #stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.5) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, size = 12,hjust = 1))
  #plt+ stat_compare_means()# method = "wilcox.test", comparisons = hide.ns = FALSE, paired = FALSE, label = "p.signif" )
  print(plt)
}

normalize_counts = function( count_matrix, cofactor = 1000)
{
  total_counts = rowSums(count_matrix)
  normalized_counts = asinh(count_matrix / total_counts * cofactor)
  return(normalized_counts)
}