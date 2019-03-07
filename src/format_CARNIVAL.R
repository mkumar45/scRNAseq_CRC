source('Z:/users/mkumar/scRNAseq_CRC/src/map_UNIPROT.R')
source("Z:/users/mkumar/scRNAseq_CRC/src/plotting/plotting_functions.R")



write_PKN = function(PKN_file,output_file)
{
  #' Import and format prior knowledge network (PKN)
  #'
  #' Filter PKN to contain only signed and directed network in format required for CARNIVAL
  #'
  #' @param PKN_file path to file containing prior knowledge network
  #' @param output_file file name for writing results
  #' @return Writes a .csv file cotaining formated PKN for CARNIVAL
  #' @export 
  omnipath_network = read.table( file = PKN_file,sep = '\t', header = TRUE)
  directed_network = subset(omnipath_network, is_directed == 1 & (is_stimulation == 1 | is_inhibition == 1))
  directed_network$sign[directed_network$is_stimulation == 1] = 1
  directed_network$sign[directed_network$is_inhibition == 1] = -1
  network = directed_network[,c("source","sign","target")]
  write.csv( network, output_file, row.names = FALSE)
}







format_CARNIVAL = function( data_directory = "data/CARNIVAL/",
                            PKN_file = "data/omnipath_PKN.sif",
                            homolog_file = "data/HOM_MouseHumanSequence.rpt",
                            PROGENy_path = "results/PROGENy/PROGENy_scores.csv",
                            TF_path ="results/Dorothea/TF_activities_confidence_AB.csv",
                            sample_file = "results/classification/processed/sample_labels.csv",
                            predicted_cell_types = "results/classification/processed/predicted_cell_types.csv"
                          )
{
  #' Create and format input files for CARNIVAL
  #'
  #' Format TF inferred activities as measurement file for CARNIVAL and
  #' format receptor ligand interactions as input file for CARNVIAL
  #'
  #' @param data_directory path to directory to write formatted CARNIVAL inputs
  #' @param PKN_file path to file containing prior knowledge network
  #' @param homolog_file path to file containing homolog information from mouse gene informatics (character)
  #' @param PROGENy_path path to file containing PROGENy pathway inferred activity
  #' @TF_path path to file containing DoRothEA inferred TF activity
  #' @param sample_file path to file containing sample/cell meta data. For subsetting. 
  #' @param predicted_cell_type path to file containing predicted cell type labels. For subsetting
  #' @return Saves input files for running CARNIVAL in "data_directory"
  #' @export 
  dir.create(data_directory)
  
  write_PKN(PKN_file, paste0( data_directory, "PKN.csv" ) )

  
  # Import data
  PROGENy_scores = read.csv(PROGENy_path, header=TRUE, row.names = 1)
  TF_activities = read.csv(TF_path,row.names=1)
  num_pathways = dim(PROGENy_scores)[2] 
  num_TF = dim(TF_activities)[2]
  
  sample_labels = read.csv(sample_file, row.names = 1, header=TRUE)
  colnames(sample_labels) = "label"
  cell_type_labels = read.csv(predicted_cell_types,row.names = 1,header=TRUE)
  
  
  split_sample = strsplit(gsub("S4","S2",as.character(sample_labels$label)),"_")
  
  PROGENy_scores$sample_type = sapply(split_sample,`[`,1)
  PROGENy_scores$tissue =sapply(split_sample,`[`,2)
  PROGENy_scores$replicate_number = sapply(split_sample,`[`,3)
  PROGENy_scores$cell_type = cell_type_labels$predicted_cell_type
  
  TF_activities$sample_type =sapply(split_sample,`[`,1)
  TF_activities$tissue =sapply(split_sample,`[`,2)
  TF_activities$replicate_number = sapply(split_sample,`[`,3)
  TF_activities$cell_type = cell_type_labels$predicted_cell_type
  

  #### Measurements are transcription factor activities
  # CAF activities
  #TF_activites = read.csv("results/Dorothea/TF_activities_confidence_AB.csv",header=TRUE,row.names = 1)
  #sample_labels = read.csv("results/classification/processed/sample_labels.csv", row.names = 1, header=TRUE)
  #predicted_cell_types = read.csv("results/classification/processed/predicted_cell_types.csv",row.names = 1,header=TRUE)
  
  # CAF_TF_activities
  CAF_TF_activites = subset( TF_activities, cell_type == "CAF" & sample_type != "WT"  & tissue == "tumor")
  numeric = unlist(lapply(CAF_TF_activites,is.numeric))
  average_CAF_TF = colMeans(CAF_TF_activites[,numeric])
  
  
  # Tumor activities
  AOM_tumor_TF = subset( TF_activities, cell_type == "Tumor" & sample_type == "AOM"  & tissue == "tumor")
  APC_tumor_TF = subset( TF_activities, cell_type == "Tumor" & sample_type == "APC"  & tissue == "tumor")
  tumor_TF = interaction_test(AOM_tumor_TF, APC_tumor_TF)

  ### Convert mouse gene names to human UNIPROT IDs for consistency with PKN
  TF_names = map_UNIPROT( colnames(TF_activities[,numeric]), homolog_file = "data/HOM_MouseHumanSequence.rpt" )
  colnames(TF_names) = c("mouse_symbol","mouse_UNIPROT","human_symbol","human_UNIPROT")
  write.csv(TF_names, paste0(data_directory,"TF_names.csv"),row.names=FALSE)
  
  average_CAF_TF = average_CAF_TF[TF_names$mouse_symbol]
  names(average_CAF_TF) = TF_names$human_UNIPROT
  
  tumor_TF = subset(tumor_TF, name %in% TF_names$mouse_symbol) 
  tumor_TF$name = TF_names$human_UNIPROT
  
  # Write TF files
  ## Use top up/down regulated TFs
  sorted_CAF_TF_activity = sort(average_CAF_TF)
  sorted_tumor_TF_activity = sort(tumor_TF$SNR)
  names(sorted_tumor_TF_activity) = tumor_TF$name
  write.table( sorted_CAF_TF_activity[abs(sorted_CAF_TF_activity) > 1.5], paste0(data_directory,"CAF_TF.csv"), sep = ",",col.names = FALSE,row.names = TRUE )
  write.table( sorted_tumor_TF_activity[abs(sorted_tumor_TF_activity) > 1.5], paste0(data_directory,"tumor_TF.csv"), sep = ",",col.names = FALSE,row.names = TRUE )
  
  #measWeights = as.matrix( sorted_TF_activity[ abs(sorted_TF_activity) > 1 ] )
  #measWeights = as.matrix( sorted_TF_activity[ abs(sorted_TF_activity) > 0.5 ] )
  
  #measWeights <- as.matrix(read_delim(file = measFile, delim = "\t", escape_double = FALSE, trim_ws = TRUE))
  #measurements <- t(sign(measWeights)) # Extracted sign of measurement for ILP fitting
  #measWeights <- t(abs(measWeights)) # Weights are all positives
  
  #data <- measurements
  #pknList <- as.data.frame(directed_network[,c("source","sign","target")])
  #colnames(pknList) <- c("Node1", "Sign", "Node2")
  
  # PROGENy scores
  #APC_PROGENy_scores = PROGENy_scores[ (predicted_cell_types == "APC" & mouse_type == "APC" & tissue_type == "tumor"),]
  #AOM_PROGENy_scores = PROGENy_scores[ (predicted_cell_types == "AOM" & mouse_type == "AOM"),]
  #CAF_PROGENy_scores = colMeans(PROGENy_scores[ (predicted_cell_types == "CAF" & mouse_type == "AOM"),])
  
  
  
  

  AOM_tumor_PROGENy = subset( PROGENy_scores, cell_type == "Tumor" & sample_type == "AOM"  & tissue == "tumor")
  APC_tumor_PROGENy = subset( PROGENy_scores, cell_type == "Tumor" & sample_type == "APC"  & tissue == "tumor")
  tumor_PROGENy_diff = interaction_test(AOM_tumor_PROGENy, APC_tumor_PROGENy)
  
  # Fibroblasts
  CAF_PROGENy = subset( PROGENy_scores, cell_type == "CAF" & sample_type == "AOM"  & tissue == "tumor")
  numeric = unlist(lapply(CAF_PROGENy,is.numeric))
  average_CAF_PROGENy = colMeans(CAF_PROGENy[,numeric])
   #colMeans(AOM_PROGENy_scores) - colMeans(APC_PROGENy_scores)
  edgeWeights = data.frame( tumor_PROGENy_diff$SNR  )
  rownames(edgeWeights) = tumor_PROGENy_diff$name
  edgeWeights = t(edgeWeights)
  CAF_PROGENy_scores = t(CAF_PROGENy_scores)
  rownames(edgeWeights) = "score"
  
  write.table( edgeWeights, paste0(data_directory,"tumor_PROGENY.csv"), sep = ",",col.names = TRUE, row.names = FALSE )
  write.table( t(average_CAF_PROGENy), paste0(data_directory,"CAF_PROGENy.csv"), sep = ",",col.names = TRUE, row.names = FALSE )
  
  ## Set input nodes as receptors or ligands from 
  RL_interactions = read.csv("results/communication/tumor_diff_interactions.csv",header=TRUE)
  RL_symbols = strsplit(as.character(RL_interactions$name),"_")
  receptor_symbols = sapply(RL_symbols,`[`,2)
  receptor_IDs = map_UNIPROT( unique(receptor_symbols), homolog_file = "data/HOM_MouseHumanSequence.rpt" )
  tumor_inputs = matrix(NaN, nrow = 1, ncol = length(receptor_IDs$human_UNIPROT))
  colnames(tumor_inputs) = receptor_IDs$human_UNIPROT
  
  #RL_interactions = read.csv("results/communication/CAF_interactions.csv",header=TRUE)
  #RL_symbols = strsplit(as.character(RL_interactions$name),"_")
  #receptor_symbols = sapply(RL_symbols,`[`,2)
  #receptor_IDs = map_UNIPROT( unique(receptor_symbols), homolog_file = "data/HOM_MouseHumanSequence.rpt" )
  #CAF_inputs = matrix(NaN, nrow = 1, ncol = length(receptor_IDs$human_UNIPROT))
  #colnames(CAF_inputs) = receptor_IDs$human_UNIPROT

  #Write input files
  write.csv(tumor_inputs, paste0( data_directory, "tumor_input_receptors.csv"),row.names = FALSE)
  #write.csv(CAF_inputs, paste0( data_directory, "CAF_input_receptors.csv"),row.names = FALSE)
  
  
  #MappedPertNode <- AddPerturbationNode(network)
  #inputs <- MappedPertNode$inputs
  #network <- MappedPertNode$network
  
  
  #if (!is.null(inputFile)) {
  #  inputs = read.table(file = "C:/Users/Manu/Dropbox (MIT)/lab/projects/julio/data/input_TGFb.txt",sep="\t",header=TRUE)
  #  inputs <- read.table(file = inputFile, sep = "\t", header = TRUE)
  #}
}