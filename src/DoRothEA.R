# Load requeired packages
require(viper) 
require(ggplot2)
require(ComplexHeatmap)
setwd("Z:/users/mkumar/scRNAseq_CRC/")

# Load TF regulon genesets in VIPER format
#load('data/TFregulons/Robjects_VIPERformat/consensus/BEST_viperRegulon.rdata')
viper_regulon = readRDS("data/TFregulons/Robjects_VIPERformat/mouse_dorothea2_regulon_v1.rds")
TF_activities = read.csv("results/TFactivities.csv", row.names = 1)
RL_list = read.csv("data/mouse_receptor_ligand.csv")

# Clean TF names & explore object
names(viper_regulon) = sapply(strsplit(names(viper_regulon), split = ' - '), head, 1)


# Explore the regulons object
names(viper_regulon)[1:10]
viper_regulon[[1]]

viper_genes = list()
for (r in 1:length(viper_regulon))
{
  viper_genes = union(viper_genes,names(viper_regulon[[r]]$tfmode))
}
viper_genes = unlist(viper_genes)
pct_ligand = list() 
pct_receptor = list()
TF_names = list()
for ( r in 1:length(viper_regulon))
{
  regulon_genes = names(viper_regulon[[r]]$tfmode)
  
  if ( "Tgfb2" %in% regulon_genes )
  {
    TF_names = append(TF_names, names(viper_regulon)[r])
  }
  
  pct_ligand[r] = sum( regulon_genes %in% RL_list$Ligand_ApprovedSymbol) / length(regulon_genes)
  pct_receptor[r] = sum( regulon_genes  %in% RL_list$Receptor_ApprovedSymbol) / length(regulon_genes)
}






TF_names = unlist(TF_names)
Heatmap( TF_activities[TF_names,predicted_cell_type=="CAF"],cluster_columns = TRUE)
cor(t(TF_activities[TF_names,predicted_cell_type=="CAF"]),expression_data[predicted_cell_type=="CAF","Tgfbr2"])

df = data.frame( unlist(pct_ligand), unlist(pct_receptor))
colnames(df) = c("pct_ligand","pct_receptor")
ggplot( data =df, aes(x = pct_ligand, y = pct_receptor)) + geom_point() + theme_bw() + coord_fixed(ratio=1)


##########################################################################################
## Example 1: Computing single-sample TF activities from a normalized gene expression matrix 
##########################################################################################
# Load expression matrix: NOTE that genes at comparable scales (e.g. zscores)
#load('data/expression/example_expressionMatrix_zscores.rdata')

sparse_matrix = readRDS("data/filtered/filtered_sparse.rds")
gene_symbols = read.table("data/combined/combined_gene_names.csv")
colnames(sparse_matrix) = gene_symbols$V1
idx = gene_symbols$V1 %in% viper_genes


# Estimate TF activities
TF_activities = viper(eset = t(scale(as.matrix(sparse_matrix)[,idx],center=TRUE,scale=TRUE)),
                      regulon = viper_regulon, nes = T, method = 'none', minsize = 4, eset.filter = F)
# Save results

write.csv(TF_activities, file = 'results/TFactivities.csv')