#!/bin/bash
#SBATCH -J filter-and-combine-ds
#SBATCH -p standard
#SBATCH -t 4-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12GB
#SBATCH --account dmobley_lab
#SBATCH --export ALL
#SBATCH --constraint=fastscratch
#SBATCH -o filter-and-combine-ds.out
#SBATCH -e filter-and-combine-ds.err

date
hostname

source ~/.bashrc
conda activate nagl2-conf-gen 

python filter-and-combine-ds.py download-opt                                           \
    --core-opt-dataset      "OpenFF NAGL2 Training Optimization Dataset Part 1 v4.0"                      \
    --core-opt-dataset      "OpenFF NAGL2 Training Optimization Dataset Part 2 v4.0"                   \
    --output                "filtered_and_combined_nagl2_opt.json"             \
    --confrmsd    "confrmsd_filtered_and_combined_nagl2_opt.json" \
    --verbose


date
