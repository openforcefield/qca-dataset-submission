#!/bin/bash
#SBATCH -J generate-combined-dataset
#SBATCH -p standard
#SBATCH -t 4-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12GB
#SBATCH --account dmobley_lab
#SBATCH --export ALL
#SBATCH --constraint=fastscratch
#SBATCH -o generate-combined-dataset.out
#SBATCH -e generate-combined-dataset.err

date
hostname

source ~/.bashrc
conda activate nagl2-conf-gen 

python generate-combined-dataset.py 

date
