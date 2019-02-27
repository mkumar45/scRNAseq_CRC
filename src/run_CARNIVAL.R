# Load necessary packages and functions
library(readr)
library(tidyr)
library(igraph)
library(devtools)

run_CARNIVAL = function( PKN_file, # file containing prior knowledge netwwork
                         measurement_file, # file containing inferred TF activity to be included in network
                         inputs_file, # file containing known perturbed nodes to be included in network 
                         weights_file = NULL, # file containing PROGENy scores to weight signaling pathways for inclusing in network
                         result_dir = "CARNIVAL/", # specify a name for result directory; if NULL, then date and time will be used by default
                         parallelCR = F, # running parallelised version?
                         inverseCR = T, # running inverse causal reasoning version (i.e. no inputs and starting from TFs)
                         nodeID = "uniprot", # specify whether the nodes are in 'uniprot' or 'gene' ID format
                         UP2GS = T, # convert UniProtIDs to Gene Symbols in the plotting step?
                         DOTfig = T, #  write DOT figure? (can be visualised by e.g. GraphViz)
                         Export_all = F, #  export all ILP variables or not; if F, only predicted node values and sif file will be written
                         mipGAP = 0.05, # in proportion to the best estimated integer solution
                         poolrelGAP = 0.0001, # in relative to the pool of best solution
                         limitPop = 500, # limit the number of populated solutions after identified best solution
                         poolCap = 100, # limit the pool size to store populated solution
                         poolIntensity = 4, # (for populating) select search intensity [0 default/ 1 to 4]
                         poolReplace = 2, # select replacement strategy of the pool solution [0 default/ 1 to 2]
                         alphaWeight = 1, # constant coefficient for fitting error in the objective function in case TF activities are not assigned [default 1]
                         betaWeight = 0.2, # relative coefficient of model size to fitting error in the objective function [default 0.2]
                         timelimit = 300 # set time limit for cplex optimisation (in seconds)
                        ) 
{
    
  load_all("CARNIVAL-master/") # load CARNIVAL package
  load(file = system.file("inst/progenyMembers.RData",package="CARNIVAL"))
  
  # Assign parallelisation parameters
  if (parallelCR)
  {
    library(doParallel)
    argsJob=commandArgs(trailingOnly = TRUE)
    repIndex <- as.numeric(argsJob[1])
    condition <- as.character(argsJob[2]) #Can additionally be used to loop over conditions of interest
  } else 
  {
    repIndex=1;condition=1
  }
  
  # Create a directory to store results
  current_dir <- getwd()
  dir.create("results",showWarnings = FALSE)
  setwd(paste(current_dir,"/results",sep=""))
  if (is.null(result_dir)) {
    dir_name <- paste("results_",Sys.time(),sep="")
  } else {
    dir_name <- result_dir
  }
  dir.create(dir_name); setwd(current_dir) 

  # Input processing
  network <- read.table(file = PKN_file, sep = ",", header = TRUE )
  #measWeights <- as.matrix(read_delim(file = measurement_file, delim = ",", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
  measWeights <- t(as.matrix(read.table(file = measurement_file, sep = ",",row.names = 1)))
  measurements <- sign(measWeights) # Extracted sign of measurement for ILP fitting
  measWeights <- abs(measWeights) # Weights are all positives
  
  if (!is.null(inputs_file)) {
    inputs <- read.table(file = inputs_file, sep = ",", header = TRUE)
  }
  if (!is.null(weights_file))
  {
    edgeWeights <- read_delim(file = weightFile, delim = ",", escape_double = FALSE, trim_ws = TRUE)
    scores <- assignPROGENyScores(progeny = edgeWeights, progenyMembers = progenyMembers, id = nodeID)
  } else 
  {
    scores <- NULL
  }
  
  # Adding perturbation node for the case of inverse causal reasoning
  if (inverseCR) 
  {
    MappedPertNode <- AddPerturbationNode(network)
    inputs <- MappedPertNode$inputs
    network <- MappedPertNode$network
  }
  data <- measurements
  pknList <- as.data.frame(network)
  colnames(pknList) <- c("Node1", "Sign", "Node2")
  print(head(pknList))
  setwd(paste0(current_dir,"/CARNIVAL-master/R/")) # temporary shift to src directory
  
  # Remove intermediate cplex files (if any)
  AllFiles <- list.files()
  CloneFiles <- which(grepl(pattern = "clone",x = AllFiles,fixed = T))
  if (length(CloneFiles)>0) {
    for (counter in 1:length(CloneFiles)) {
      file.remove(AllFiles[CloneFiles[counter]])
    }
  }
  
  # Remove redundant files prior to optimisation
  if (file.exists(paste0("testFile_",condition,"_",repIndex,".lp"))) {file.remove(paste0("testFile_",condition,"_",repIndex,".lp"))} # might be useful for debugging
  if (file.exists(paste0("results_cplex_",condition, "_",repIndex,".txt"))) {file.remove(paste0("results_cplex_",condition,"_",repIndex,".txt"))}
  if (file.exists("cplex.log")) {file.remove("cplex.log")}
  if (file.exists(paste0("cplexCommand_", condition,"_",repIndex,".txt"))) {file.remove(paste0("cplexCommand_", condition,"_",repIndex,".txt"))}
  
  # Write constraints as ILP inputs
  ptm <- proc.time()
  print("Writing constraints...")
  variables <- writeLPFile(data=measurements,pknList=pknList,inputs=inputs,betaWeight=betaWeight,
                           scores=scores,mipGAP=mipGAP,poolrelGAP=poolrelGAP,limitPop=limitPop,
                           poolCap=poolCap,poolIntensity=poolIntensity,poolReplace=poolReplace,
                           timelimit=timelimit,measWeights=measWeights,
                           repIndex=repIndex,condition = condition)
  Elapsed_1 <- proc.time() - ptm
  
  # Solve ILP problem with cplex, remove temp files, and return to the main directory
  ptm <- proc.time()
  print("Solving LP problem...")
  #system("cmd.exe",input=paste0("\"",getwd(), "/cplex -f cplexCommand_", condition,"_",repIndex,".txt","\""))
  system("cmd.exe",input=paste0("cplex -f cplexCommand_", condition,"_",repIndex,".txt","\""))
  Elapsed_2 <- proc.time() - ptm
  
  # Move result files to result folder and remove redundant files after the optimisation
  if (file.exists(paste0("testFile_",condition,"_",repIndex,".lp"))) 
    {
    file.remove(paste0("testFile_",condition,"_",repIndex,".lp"))
    } # might be useful for debugging
  if (file.exists(paste0("results_cplex_",condition, "_",repIndex,".txt"))) 
  {
    file.copy(from = paste0("results_cplex_",condition,"_",repIndex,".txt"),to = paste(current_dir,"/results/",dir_name,"/results_cplex.txt",sep=""),overwrite = TRUE)
    file.remove(paste0("results_cplex_",condition,"_",repIndex,".txt"))
    }
  if (file.exists("cplex.log"))
  {
    file.copy(from = "cplex.log",to = paste0(current_dir,"/results/",dir_name,"/cplex_",condition,"_",repIndex,".log"), overwrite = TRUE); file.remove("cplex.log")
    }
  if (file.exists(paste0("cplexCommand_", condition,"_",repIndex,".txt"))) 
  {
    file.remove(paste0("cplexCommand_", condition,"_",repIndex,".txt"))
  }
  
  
  setwd(current_dir)
  
  # Write result files in the results folder
  ptm <- proc.time()
  print("Writing result files...")
  if (file.exists(paste0("results/",dir_name,"/results_cplex.txt"))) {
    for(i in 1:length(variables)){
      res <- exportResult(cplexSolutionFileName = paste0("results/",dir_name,"/results_cplex.txt"),
                          variables = variables, pknList = pknList, conditionIDX = i,
                          dir_name = dir_name, inputs=inputs,measurements=measurements,
                          Export_all = Export_all,writeIndividualResults = T)
      # res <- files2res(counterlist) # retrieve results from previously generated result files
    }
    if (!is.null(res)) {
      if (UP2GS) {res <- Uniprot2GeneSymbol(res)}
      if (DOTfig) {WriteDOTfig(res=res,dir_name=dir_name,
                               inputs=inputs,measurements=measurements,UP2GS=UP2GS)}
      # if (DOTfig) {WriteDOTfig(res=res,idxModel=c(1,2),dir_name=dir_name,
      #                             inputs=inputs,measurements=measurements,UP2GS=UP2GS)}
      save(res,file = paste0(dir_name,"/results_CARNIVAL.Rdata"))
    }
  } else {
    print("No result to be written")
  }
  Elapsed_3 <- proc.time() - ptm
  
  #file.remove(paste0(dir_name,"/results_cplex.txt")) # optional; remove cplex results (to save space)
  
  # Logged computational time
  ElapsedAll <- as.data.frame(matrix(t(c(Elapsed_1[3],Elapsed_2[3],Elapsed_3[3])),3,1))
  rownames(ElapsedAll) <- c("WriteConstraints:","CplexSolving:","ExportResults:")
  write.table(x = ElapsedAll,file = paste0("results/",dir_name,"/elapsed_time.txt"),col.names = F,row.names = T,quote = F)
  
  # --- End of script --- #
  
}
args = commandArgs(trailingOnly=TRUE)
run_CARNIVAL( args[1], args[2], args[3] )