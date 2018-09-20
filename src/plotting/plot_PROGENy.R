require("ggplot2")
require("circlize")
#require("ComplexHeatmap")
setwd("Z:/users/mkumar/scRNAseq_CRC//")

pathway_scores = read.csv("results/PROGENy/PROGENy_scores.csv",header=TRUE )
pathway_names = colnames(pathway_scores)[-1]
sample_labels = read.csv("data/combined/combined_labels.csv",header=FALSE)
cell_type_labels = read.csv("results/classification/predicted_cell_types.csv",header=FALSE)
split_sample = strsplit(gsub("S4","S2",as.character(sample_labels$V2)),"_")
pathway_scores$sample_type =sapply(split_sample,`[`,1)
pathway_scores$replicate_number = sapply(split_sample,`[`,3)
pathway_scores$cell_type = cell_type_labels$V1


cell_types = levels(unique(pathway_scores$cell_type))
num_cell_types = length(cell_types)
num_pathways = length(pathway_names)

p_val = matrix(NA, nrow = num_cell_types, ncol= num_pathways)
test_stat = matrix(NA, nrow = num_cell_types, ncol= num_pathways)
fold_change = matrix(NA, nrow = num_cell_types, ncol= num_pathways)


pathway_scores = subset(pathway_scores, (cell_type != "APC" & cell_type != "AOM") | ( cell_type == "APC" & sample_type == "APC") | ( cell_type == "AOM" & sample_type == "AOM")  )


for (p in 1:num_pathways)
{
  pdf(paste("figures/PROGENy/",pathway_names[p],".pdf"),width=11,height=8.5,paper='special') # open pdf 
  
  bplot = ggplot( pathway_scores, aes_string( x = "cell_type", y = pathway_names[p], color="sample_type")) + geom_point(size=0.25,position = position_jitterdodge()) + 
    #scale_y_continuous( limits = c(-0.25,1), breaks = c(-0.05,0,0.25,0.5,1)) +
    stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.5) +
    theme_bw() 
  print(bplot)
  dev.off()
}
  ## Stastical tests
  for (n in 1:num_cell_types)
  {   
    #apc_scores = subset(pathway_scores, cell_type == cell_types[n] & sample_type == "Apc")[,pathway_names[p]]
    #aom_scores = subset(pathway_scores, cell_type == cell_types[n] & sample_type == "AOM")[,pathway_names[p]]
    
    #print(length(apc_scores))
    #print(length(aom_scores))
    #print(length(apc_scores) > 0 & length(aom_scores) > 0)
    
    if ( length(apc_scores) > 0 & length(aom_scores) > 0)
    {
      #test = wilcox.test( apc_scores, aom_scores )
      #print(paste("pathway:",pathway_names[p],".","cell_type:", cell_types[n]))
      #fold_change[n,p] = log2(median( apc_scores)/median(aom_scores))
      #test_stat[n,p] = test$statistic
      #p_val[n,p] = test$p.value
      
    }
  }
}
#test = wilcox.test( subset(pathway_scores, cell_type == "CAF" & sample_type == "Apc")$TGFb,
#                    subset(pathway_scores, cell_type == "CAF" & sample_type == "AOM")$TGFb)


# # Heatmap of fold-change
# rownames(fold_change) = cell_types
# colnames(fold_change) = pathway_names
#      
# 
# cVals = seq( from = -10, to = 10, by = .5 )
# col_func = colorRampPalette(c("red","white","blue")) # Use red-white color scale for all heatmaps
# color_map = colorRamp2(cVals,col_func(length(cVals)))    
# Heatmap(fold_change, na_col = "grey",cluster_rows = FALSE,cluster_columns = FALSE, 
#         col = color_map, name = "hMap")
# #### Annotate significant interactions
# num_row = dim(fold_change)[1]
# num_col = dim(fold_change)[2]
# significant_idx = arrayInd( which( p_val < .001 ), dim(fold_change)) # indices of significant scores
# sig_row_idx = (num_row + 1) - significant_idx[,1]  # row_idx starts from bottom left
# sig_col_idx = significant_idx[,2]
# # Convert to units for figure
# x_coord = unit((sig_col_idx-0.5)/num_col, "npc")
# y_coord = unit((sig_row_idx-0.5)/num_row, "npc")
# width = unit(1/num_col, "npc")
# height = unit(1/num_row, "npc")
# #decorate_heatmap_body("hMap",  code = {dec = grid.circle(x_coord, y_coord, r = unit(0.25, "mm"),
# #                                                             default.units = "npc", name = NULL,
# #                                                             gp=gpar(col = NULL,alpha=1,lwd=1,lty="solid",fill="black"), draw = TRUE, vp = NULL) })
# 
# 
# 
# 
# 
# 
# 
# # Volcano plot
# df = data.frame( as.vector(fold_change), as.vector(-log10(p_val)), rep(colnames(fold_change),each=dim(fold_change)[1]),rep(rownames(fold_change),dim(fold_change)[2]) )
# colnames(df)= c("fold_change","p_value","pathway","cell_type")
# 
# plt = ggplot( data = df, aes( y = p_value, x = fold_change, color = pathway)) + geom_point() + 
#       scale_y_continuous(limits = c(-1,10))
#       theme_bw()
# print(plt)