#!/bin/bash
#SBATCH -N 1                      # Number of nodes. You must always set -N 1 unless you receive special instruction from the system admin
#SBATCH -n 16                      # Number of CPUs. Equivalent to the -pe whole_nodes 1 option in SGE
#SBATCH --mail-type=END           # Type of email notification- BEGIN,END,FAIL,ALL. Equivalent to the -m option in SGE 
#SBATCH --mail-user=mkumar@mit.edu  # Email to which notifications will be sent.

module load python3 

for pct_var in $(seq 30 20 90); do
	for num_gene in $( seq 7 11 ); do
		sbatch  src/train_classifier.sh $pct_var $((2**$num_gene))
	done
done