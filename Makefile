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
PYTHON=python3

# Define filenames
COUNT_FILES=$(wildcard data/filtered_counts/*.csv)
SPARSE_FILES=$(patsubst data/filtered_counts/%.csv, data/sparse/%.npz, $(COUNT_FILES))

.PHONY : variables
variables: 
	@echo COUNT_FILES: $(COUNT_FILES)
	#@echo SPARSE_FILES: $(SPARSE_FILES)

	
#.PHONY : quality_control	
#quality_control : $(QC_FILES)

.PHONY : all
all : sparse_samples classification processed_predictions

.PHONY : sparse_samples
sparse_samples : $(SPARSE_FILES)


## data/qc/%_sparse.npz : import individual samples and process into sparse format
data/sparse/%.npz : src/qc.py data/filtered_counts/%.csv data/mouse_gene_symbols_unfiltered.csv
	mkdir -p results/UMAP results/DPC figures/UMAP data/sparse 
	$(PYTHON) src/format_sparse.py data/filtered_counts/$*.csv data/mouse_gene_symbols_unfiltered.csv
	#sbatch src/qc.sh $* for running on mit luria cluster

### data/qc/%_sparse.npz : import individual samples and process into sparse format
#data/qc/%_sparse.npz : src/qc.py data/counts/%.tsv
#	mkdir -p results/qc figures/qc data/qc
#	$(PYTHON) src/qc.py data/counts/$*.tsv
#	#sbatch src/qc.sh $* for running on mit luria cluster
	
# Combination of data and quality control is done in notebook 

# Cell type classification #
.PHONY : classification 
classification : training_data prediction
training_data : results/classification/sparse_training_counts.npz
results/classification/sparse_training_counts.npz :  src/create_training_data.py data/combined/processed_counts.npz data/cell_type_markers.csv 
	$(PYTHON) $< data/combined/processed_counts.npz data/mouse_gene_symbols.csv data/cell_type_markers.csv
trained_classifier : training_data src/train_classifier.py 
#results/classification/sparse_training_counts.npz results/classification/training_labels.csv
	mkdir -p results/classification/training_results
	sh src/test_parameters.sh
prediction: results/classification/predicted_cell_types.csv
results/classification/predicted_cell_types.csv : src/predict_cell_type.py results/classification/sparse_training_counts.npz
	$(PYTHON) src/predict_cell_type.py results/classification/sparse_training_counts.npz results/classification/training_labels.csv data/combined/processed_counts.npz -n 1024 -pca 90

.PHONY : processed_predictions
processed_predictions : results/classification/processed/processed_counts.npz
results/classification/processed/processed_counts.npz : src/process_predictions.py results/classification/predicted_cell_types.csv
	mkdir -p results/classification/processed/
	$(PYTHON) src/process_predictions.py
	
#
# Interaction scores
results/communication/interaction_scores.csv : src/cell_cell_communication.py results/classification/processed/processed_counts.npz
#data/combined/combined_sparse.npz data/combined/combined_gene_names.csv results/classification/predicted_cell_types.csv data/mouse_receptor_ligand.csv  data/combined/combined_labels.csv
	mkdir -p results/communication/
	$(PYTHON) src/cell_cell_communication.py results/classification/processed/processed_counts.npz data/mouse_gene_symbols.csv results/classification/processed/predicted_cell_types.csv data/mouse_receptor_ligand.csv  $@ -g results/classification/processed/sample_labels.csv

#
# Pathway scores
results/PROGENy/PROGENy_scores.csv : src/compute_PROGENy_scores.py results/classification/processed/processed_counts.npz data/mouse_gene_symbols.csv results/classification/processed/predicted_cell_types.csv data/PROGENy_mouse_model_v2.csv
	mkdir -p results/PROGENy/
	$(PYTHON) $< results/classification/processed/processed_counts.npz data/mouse_gene_symbols.csv results/classification/processed/predicted_cell_types.csv data/PROGENy_mouse_model_v2.csv $@
  
# TF scores
results/Dorothea/TF_activities.csv : src/run_Dorothea.R results/classification/processed/processed_counts.RDS data/TFregulons/Robjects_VIPERformat/mouse_dorothea2_regulon_v1.rds data/mouse_gene_symbols.csv
	Rscript src/run_Dorothea.R results/classification/processed/processed_counts.RDS data/TFregulons/Robjects_VIPERformat/mouse_dorothea2_regulon_v1.rds data/mouse_gene_symbols.csv

#
# CARNIVAL
results/CARNIVAL :
	mkdir results/CARNIVAL/
	Rscript src/format_CARNIVAL.R
	Rscript src/run_CARNIVAL.R
	Rscript src/map_UNIPROT.R
	
# figures
figures/tsne_cell_type.pdf figures/tsne_sample.pdf : src/plotting/tsne_plots.py results/tSNE/combined_sparse_tsne.csv
	$(PYTHON) $< results/tSNE/combined_sparse_tsne.csv
	
results/classification/processed/processed_counts.RDS : results/classification/processed/processed_counts.npz src/convert_sparse.py src/make_sparse_RDS.R  
	$(PYTHON) src/convert_sparse.py $< results/classification/processed/processed_counts.csv
	Rscript src/make_sparse_RDS.R results/classification/processed/processed_counts.csv $@

%.RDS : %.npz src/convert_sparse.py src/make_sparse_RDS.R  data/mouse_gene_symbols.csv
	$(PYTHON) src/convert_sparse.py $< $*.csv
	Rscript src/make_sparse_RDS.R $*.csv data/mouse_gene_symbols.csv $@
