library("Matrix")

make_sparse = function( expression_data, gene_symbols = "data/mouse_gene_symbols.csv", output_file )
{
  gene_symbols = read.table(gene_symbols)
  sparse_matrix = Matrix( data = as.matrix(read.csv(expression_data,header=FALSE)) , sparse = TRUE)
  colnames(sparse_matrix) = gene_symbols$V1
  saveRDS( sparse_matrix, output_file)
}

args = commandArgs(trailingOnly=TRUE)
make_sparse( args[1], args[2], args[3] )
