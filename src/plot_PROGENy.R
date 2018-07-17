require("ggplot2")
setwd("Z:/users/mkumar/Lau/")

pathway_scores = read.csv("results/PROGENy_scores.csv",header=TRUE )
pathway_names = colnames(pathway_scores)[-1]
sample_labels = read.csv("data/combined_labels.csv",header=FALSE)
cell_type_labels = read.csv("results/predicted_cell_types.csv",header=FALSE)
split_sample = strsplit(gsub("S4","S2",as.character(sample_labels$V1)),"_")
pathway_scores$sample_type =sapply(split_sample,`[`,1)
pathway_scores$replicate_number = sapply(split_sample,`[`,3)
pathway_scores$cell_type = cell_type_labels$V1




##
bplot = ggplot( pathway_scores, aes( x = cell_type, y = TGFb, color=sample_type)) + geom_point(binaxis = "y",stackdir = "center",size=0.25,binwidth = 0.02,width = 0.05,position = position_jitterdodge()) + 
  #scale_y_continuous( limits = c(-0.25,1), breaks = c(-0.05,0,0.25,0.5,1)) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.5) +
  theme_bw() 
print(bplot)

test = wilcox.test( subset(pathway_scores, cell_type == "CAF" & sample_type == "Apc")$TGFb,
                    subset(pathway_scores, cell_type == "CAF" & sample_type == "AOM")$TGFb)





for (p in 1:num_pathways)
{

  bplot = ggplot( pathway_scores, aes_string( x = "cell_type", y = pathway_names[p], color="sample_type")) + geom_point(binaxis = "y",stackdir = "center",size=0.25,binwidth = 0.02,width = 0.05,position = position_jitterdodge()) + 
    #scale_y_continuous( limits = c(-0.25,1), breaks = c(-0.05,0,0.25,0.5,1)) +
    stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.5) +
    theme_bw() 
  print(bplot)
}

## Stastical tests
cell_types = levels(unique(pathway_scores$cell_type))
num_cell_types = length(cell_types)
num_pathways = length(pathway_names)

diff = matrix(NA, nrow = num_cell_types, ncol= num_pathways)
for (n in 1:num_cell_types){
  for (p in 1:num_pathways){

    
    apc_scores = subset(pathway_scores, cell_type == cell_types[n] & sample_type == "Apc")[,pathway_names[p]]
    aom_scores = subset(pathway_scores, cell_type == cell_types[n] & sample_type == "AOM")[,pathway_names[p]]
    
    if ( length(apc_scores) > 0 & length(aom_scores) > 0)
    {
      test = wilcox.test( apc_scores, aom_scores )
      
      if (test$p.value < 0.001)
      {
        print(paste("pathway:",pathway_names[p],".","cell_type:", cell_types[n]))
        diff[n,p] = 1
        
      }
    } 
    
  }
}

rownames(diff) = cell_types
colnames(diff) = pathway_names
         
Heatmap(diff, na_col = "grey",cluster_rows = FALSE,cluster_columns = FALSE)