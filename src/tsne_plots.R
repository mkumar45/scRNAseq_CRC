require("ggplot2")
setwd("Z:/users/mkumar/Lau/")

tsne_coordinates = read.csv("results/combined_sparse_tsne.csv",header=FALSE)
colnames(tsne_coordinates) = c("tSNE1","tSNE2")
cell_type_predictions = read.csv("results/predicted_cell_types.csv",header=FALSE)
colnames(cell_type_predictions) = "cell_type"
prediction_probabilities = read.csv("results/prediction_probabilities.csv",header=TRUE)
max_probability = apply( prediction_probabilities, MARGIN = 1, max)
replicate_name = read.csv("data/combined_labels.csv",header=FALSE)
split_sample = strsplit(gsub("S4","S2",as.character(replicate_name$V1)),"_")
sample_type =sapply(split_sample,`[`,1)

#df = cbind( tsne_coordinates, cell_type_predictions, prediction_probabilities, max_probability, sample_type, replicate_name)
df = cbind( tsne_coordinates, sample_type, replicate_name)

p1 = ggplot(data = df, aes( x = tSNE1, y = tSNE2, color = cell_type )) + geom_point() + theme_bw()
print(p1)
p2 = ggplot(data = df, aes( x = tSNE1, y = tSNE2, color = sample_type )) + geom_point() + theme_bw()
print(p2)


p3 = ggplot(data = df, aes( x = tSNE1, y = tSNE2, color = max_probability )) + geom_point() + theme_bw() + coord_fixed()
print(p3)

