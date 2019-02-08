# Makefile for cell-cell communication analysis

.PHONY : help
help : Makefile
	@sed -n 's/^##//p' $<

# $@ = target of the current rule
# $^ = all dependencies
# $< = first dependices of current rule
# % = wildcard (use in dependencies)
# $* = replaced by the stem the rule matched (use in actions)

# Define variables + scripts for pre-processing #
LANGUAGE=python3

# Define filenames
COUNT_FILES=$(wildcard data/counts/*.tsv)
QC_FILES=$(patsubst data/counts/%.tsv, data/qc/%_sparse.npz, $(COUNT_FILES))

.PHONY : variables
variables: 
	@echo COUNT_FILES: $(COUNT_FILES)
	@echo QC_FILES: $(QC_FILES)
	@echo METRIC_FILES: $(METRIC_FILES)
	@echo SPARSE_FILES: $(SPARSE_FILES)
	@echo SPARSE_R_FILES: $(SPARSE_R_FILES)
	@echo GENE_NAMES: $(GENE_NAMES)

	
.PHONY : quality_control	
quality_control : $(QC_FILES)
	
## data/qc/%_sparse.npz : import individual samples and process into sparse format
data/qc/%_sparse.npz : src/qc.py data/counts/%.tsv
	mkdir -p results/qc figures/qc data/qc
	$(LANGUAGE) src/qc.py data/counts/$*.tsv
	#sbatch src/qc.sh $* for running on mit luria cluster
	
# Combination of data and quality control is done in notebook 

# Cell type classification #
.PHONY : classification 
classification : training_data trained_classifier  prediction
training_data : results/classification/sparse_training_counts.npz
results/classification/sparse_training_counts.npz :  src/create_training_data.py data/combined/processed_counts.npz data/cell_type_markers.csv 
	$(LANGUAGE) $< data/combined/processed_counts.npz data/mouse_gene_symbols.csv data/cell_type_markers.csv
trained_classifier : training_data src/train_classifier.py 
#results/classification/sparse_training_counts.npz results/classification/training_labels.csv
	mkdir -p results/classification/training_results
	sh src/test_parameters.sh
prediction: results/classification/predicted_cell_types.csv
results/classification/predicted_cell_types.csv : src/predict_cell_type.py trained_classifier training_data
	$(LANGUAGE) src/predict_cell_type.py results/classification/sparse_training_counts.npz results/classification/training_labels.csv data/combined/processed_counts.npz -n 1024 -pca 90
	

#
# Pathway scores
results/PROGENy/PROGENy_scores.csv : src/compute_progeny_pathway_scores.py data/combined/combined_sparse.npz data/combined/combined_gene_names.csv results/classification/predicted_cell_types.csv data/PROGENy_mouse_model_v2.csv
	mkdir -p results/PROGENy/
	$(LANGUAGE) $< data/combined/combined_sparse.npz data/combined/combined_gene_names.csv results/classification/predicted_cell_types.csv data/PROGENy_mouse_model_v2.csv $@

#
# Interaction scores
results/communication/interaction_scores.csv : src/cell_cell_communication.py data/combined/combined_sparse.npz data/combined/combined_gene_names.csv results/classification/predicted_cell_types.csv data/mouse_receptor_ligand.csv  data/combined/combined_labels.csv
	mkdir -p results/communication/
	$(LANGUAGE) $< data/combined/combined_sparse.npz data/combined/combined_gene_names.csv results/classification/predicted_cell_types.csv data/mouse_receptor_ligand.csv  $@ -g data/combined/combined_labels.csv
#
# CARNIVAL
results/CARNIVAL :
	mkdir results/CARNIVAL/
	Rscript format_CARNIVAL
	Rscript run_CARNIVAL
	Rscript map_UNIPROT
	
# figures
figures/tsne_cell_type.pdf figures/tsne_sample.pdf : src/plotting/tsne_plots.py results/tSNE/combined_sparse_tsne.csv
	$(LANGUAGE) $< results/tSNE/combined_sparse_tsne.csv
	