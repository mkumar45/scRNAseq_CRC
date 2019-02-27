#!/bin/bash
#SBATCH -N 1                      # Number of nodes. You must always set -N 1 unless you receive special instruction from the system admin
#SBATCH -n 16                      # Number of CPUs. Equivalent to the -pe whole_nodes 1 option in SGE
#SBATCH --mail-type=END           # Type of email notification- BEGIN,END,FAIL,ALL. Equivalent to the -m option in SGE 
#SBATCH --mail-user=mkumar@mit.edu  # Email to which notifications will be sent.

module load python3 
python3 src/qc.py "data/counts/$1.tsv"
