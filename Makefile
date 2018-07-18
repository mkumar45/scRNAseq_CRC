# Makefile for cell-cell communication analysis

.PHONY : help
help : Makefile
	@sed -n 's/^##//p' $<

# $@ = target of the current rule
# $^ = all dependencies
# $< = first dependices of current rule
# % = wildcard (use in dependencies)
# $* = replaced by the stem the rule matched (use in actions)



### Define variables + scripts for pre-processing ### 
LANGUAGE=python3
PROCESS_SRC=src/process_raw_data.py
PROCESS_EXE=$(LANGUAGE) $(PROCESS_SRC)

TSNE_SRC=src/run_tsne.py
TSNE_EXE=$(LANGUAGE) $(TSNE_SRC)

### Variables + scripts for classification ###
CLASSIFY_SRC=code/cell_type_classification.py
CLASSIFY_EXE=$(LANGUAGE) $(CLASSIFY_SRC)

# Define filenames
TUMOR_SAMPLES=$(wildcard data/raw/*.csv)
PROCESSED_FILES=$(patsubst data/raw/%.csv, data/processed/%_sparse.npz, $(TUMOR_SAMPLES))
SPARSE_FILES=$(patsubst data/raw/%.csv, data/processed/%_sparse.npz, $(TUMOR_SAMPLES))
GENE_NAMES=$(patsubst data/raw/%.csv, data/processed/%_gene_names.csv, $(TUMOR_SAMPLES))


.PHONY : variables
variables: 
	@echo TUMOR_SAMPLES: $(TUMOR_SAMPLES)
	@echo SPARSE_FILES: $(SPARSE_FILES)
	@echo GENE_NAMES: $(GENE_NAMES)
	

#	
# Format all tumor dataset 
#
## processed_data : read in csv files and save sparse representation and gene names for each sample. 
.PHONY : processed_data
processed_data : $(SPARSE_FILES) $(GENE_NAMES) 

data/processed/ :
	mkdir -p $@
data/combined/ :
	mkdir -p $@
data/filtered/ :
	mkdir -p $@
	
## data/processed/%_sparse.npz data/processed/%_gene_names.csv : read individual samples and process into sparse format
data/processed/%_sparse.npz data/processed/%_gene_names.csv : data/raw/%.csv $(PROCESS_SRC) | data/processed/
	$(PROCESS_EXE) $< data/processed/$*_sparse.npz data/processed/$*_gene_names.csv 

## data/combined_sparse.npz data/combined_labels.csv : combine all samples into a single sparse representation
data/combined/combined_sparse.npz data/combined/combined_labels.csv : src/combine_data.py $(SPARSE_FILES) $(GENE_NAMES) data/combined/
	$(LANGUAGE) $< data/processed/ $@
	
## data/combined_sparse.npz data/combined_labels.csv : combine all samples into a single sparse representation
data/filtered/filtered_sparse.npz data/filtered/filtered_labels.csv results/qc_info.csv: src/data_qc.py  data/combined/combined_sparse.npz data/mitochondrial_genes_MGI.txt data/filtered/
	$(LANGUAGE) $< data/combined/combined_sparse.npz data/combined/combined_gene_names.csv data/combined/combined_labels.csv results/combined_sparse_tsne.csv data/mitochondrial_genes_MGI.txt 	

#
# t-SNE
#	
## results/%_tsne.csv : run tsne using sparse representation as input
results/%_tsne.csv : data/%.npz src/bhtsne/bh_tsne
	$(TSNE_EXE) $< $@ -n 1000 

## src/bhtsne/bh_tsne : compile tsne executable
src/bhtsne/bh_tsne : src/bhtsne/sptree.cpp src/bhtsne/tsne.cpp src/bhtsne/tsne_main.cpp
	g++ src/bhtsne/sptree.cpp src/bhtsne/tsne.cpp src/bhtsne/tsne_main.cpp -o src/bhtsne/bh_tsne -O2

#
# Cell type classification ###
#
.PHONY : classification
classification : $(GMM_FILES) $(BIC_FILES)
results/classification/%_classification : data/tumor_data/%_processed $(CLASSIFY_SRC)
	@echo classification
	$(CLASSIFY_EXE) $< $@ 


