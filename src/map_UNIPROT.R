map_UNIPROT = function( input_SYMBOLS, 
                        homolog_file = "data/HOM_MouseHumanSequence.rpt"
                      )
{
  #' Map list of mouse gene symbols to mouse UNIPROT ID, human gene symbols, and human UNIPROT IDs
  #'
  #' CARNIVAL prior knowldegde network (PKN) uses UNIPROT IDs so gene symbols need to be mapped before running carnvial for consistency
  #'
  #' @param input_SYMBOLS mouse gene symbols to convert (list)
  #' @param homolog_file Path to file containing homolog information from mouse gene informatics (character)
  #' @return A data frame containing mouse and human gene symbols and UNIPROT IDs for every mouse symbol found in MGI database
  #' @export 
  homolog_table = read.table(homolog_file,sep="\t",header=TRUE)
  human_SYMBOLS = list()
  human_UNIPROT = list()
  mouse_UNIPROT = list()
  mouse_SYMBOLS = list()
  for (sym in input_SYMBOLS)
  {
    idx = homolog_table$Symbol == sym # Find mouse gene symbol in homolog table
    mouse_prot = as.character(homolog_table[idx,]$SWISS_PROT.IDs) # get corresponding UNIPROT ID
    
    if (identical(mouse_prot,character(0)))
    {
      print(paste0("No uniprot ID found for " ,sym))
      next
    }
    mouse_SYMBOLS = append(mouse_SYMBOLS,sym[1])
    mouse_UNIPROT = append(mouse_UNIPROT,mouse_prot)

    
    homolog_id = homolog_ta?map_ble[idx,]$HomoloGene.ID # Get homolog ID correspond to mouse gene
    
    df = subset(homolog_table, HomoloGene.ID == homolog_id & Common.Organism.Name == "human"  ) # Find human genes with same homolog ID
    
    human_name = as.character(df$Symbol[1])
    human_prot = as.character(df$SWISS_PROT.IDs[1]) 
    
    human_UNIPROT = append(human_UNIPROT,human_prot)
    human_SYMBOLS = append(human_SYMBOLS,human_name)
  }
  mouse_SYMBOLS = unlist(mouse_SYMBOLS)
  human_SYMBOLS = unlist(human_SYMBOLS)
  mouse_UNIPROT = unlist(mouse_UNIPROT)
  human_UNIPROT = unlist(human_UNIPROT)
  
  df = data.frame( mouse_SYMBOLS, human_SYMBOLS, mouse_UNIPROT, human_UNIPROT)
  return(df)
}