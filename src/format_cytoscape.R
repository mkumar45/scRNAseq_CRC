format_cytoscape = function(  results_dir = "CARNIVAL/",  
                         homolog_file = "data/HOM_MouseHumanSequence.rpt", 
                         RL_interaction_file = "data/mouse_receptor_ligand.csv",
                         TF_regulon_file = "data/TFregulons/Robjects_VIPERformat/mouse_dorothea2_regulon_v1.rds"
                      ) 
{
  #' Format CARNIVAL outputs for visualing in cytoscae
  #'
  #' Create edge list and node list that can be imported directly into cytoscape
  #'
  #' @param results_dir path to directory containing CARNIVAL output results
  #' @param homolog_file path to file containing homolog information from mouse gene informatics (character)
  #' @param RL_interaction_file path to file containing receptor/ligand interactions for annotating nodes (character)
  #' @param TF_regulon_file path to file containing DoRothEA transcription factor regulons to use for VIPER inference for annotating nodes
  #' @return Writes two .csv files cotaining edge list and node list 
  #' @export 
  weighted_network = read.table(paste0(results_dir,"weightedModel_1.txt"),header=TRUE,stringsAsFactors = FALSE)
  node_attributes = read.table(paste0(results_dir,"nodesAttributes_1.txt"),header=TRUE,fill=TRUE,stringsAsFactors = FALSE)
  
  homolog_table = read.table(homolog_file,sep="\t",header=TRUE)
  
  ## Map UNIPROT IDs to Gene Symbols ##
  unique_UNIPROT = unique( union( as.character(weighted_network$Node1),as.character(weighted_network$Node2)))
  unique_SYMBOL = list()
  mapped_network = weighted_network
  mapped_nodes = node_attributes
  for (uniprot in unique_UNIPROT)
  {
    idx = homolog_table$SWISS_PROT.IDs == uniprot
    if (any(idx))  
    {
      symbol = as.character(homolog_table[idx,]$Symbol)  
  
      
    } else 
    { 
      symbol = uniprot
      print(paste(uniprot, "not mapped"))
    }
    mapped_network$Node1[ weighted_network$Node1 == uniprot] = symbol
    mapped_network$Node2[ weighted_network$Node2 == uniprot] = symbol
    mapped_nodes$Node[ node_attributes$Node == uniprot] = symbol
    
    unique_SYMBOL = append(unique_SYMBOL,symbol)
  }
  
  ## Annotate nodes as receptor,ligand, or transcription factor
  RL_interactions = read.csv(RL_interaction_file,header=TRUE)
  TF_regulons = readRDS(TF_regulon_file)
  split_names = strsplit( as.character(RL_interactions$Pair_Name),"_")
  ligand_name = sapply(split_names,'[',1)
  ligand_name_mouse = RL_interactions$Ligand_ApprovedSymbol
  receptor_name = sapply(split_names,'[',2)
  receptor_name_mouse = RL_interactions$Receptor_ApprovedSymbol
    
  split_TF = strsplit( as.character(names(TF_regulons)),"_")
  mouse_TF_names = sapply(split_TF,'[',1)
   
  mapped_network$Interaction_type = "PPI"
  for (sym in unique_SYMBOL)
  {
    if ( sym %in% ligand_name )  { mapped_nodes$nodesP[mapped_nodes$Node==sym] = "L"}
    if ( sym %in% receptor_name )  { mapped_nodes$nodesP[mapped_nodes$Node==sym] = "R"}
    if ( sym %in% TF_names$human_symbol )  
    { 
      mapped_nodes$nodesP[mapped_nodes$Node==sym] = "TF"
      #mouse_TF = TF_names$mouse_symbol[ TF_names$human_symbol == sym]
      #idx = mouse_TF_names == as.character(mouse_TF)
      #regulon = TF_regulons[[which(idx)]]
      
      #RL_idx = (names(regulon$tfmode) %in% ligand_name_mouse | names(regulon$tfmode) %in% receptor_name_mouse)
      #if (any(RL_idx))
      #{
      #  TF_network = data.frame(matrix(NA,length(regulon$tfmode[RL_idx]),4))
      #  TF_network[,1] = as.character(sym) # source is the transcription factor
      #  TF_network[,2] = regulon$tfmode[RL_idx] # sign of interaction
      #  TF_network[,3] = names(regulon$tfmode[RL_idx]) # target
      #  TF_network[,4] = 100 # Arbitrary Weight
      #  TF_network[,5] = "TF" # Inferred TF interaction, not PPI
      #  colnames(TF_network) = colnames(mapped_network)
      #  mapped_network = rbind(mapped_network,TF_network)
      #}
    }
  }
  
  ## Write files
  write.csv(mapped_network,"results/CARNIVAL/weightedModel_GeneSymbols.csv",row.names = FALSE,quote = FALSE)
  write.csv(mapped_nodes,"results/CARNIVAL/nodeAttributes_GeneSymbols.csv",row.names = FALSE,quote = FALSE)
}