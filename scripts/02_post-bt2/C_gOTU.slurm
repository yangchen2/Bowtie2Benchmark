#!/bin/bash -l
# AUTHOR: Yang Chen (yac027@ucsd.edu)
# DATE: 10/01/25
# DESCRIPTION: Map collapsed BIOM tables (phylum, genus, species) to gOTUs
#              Produces *_collapsed_gOTU.biom in each TBL_DIR.

#SBATCH --job-name=C_gOTU
#SBATCH --output=/ddn_scratch/yac027/03_Bowtie2Benchmark/slurm_out/Test1/C_gOTU_%A.out
#SBATCH --error=/ddn_scratch/yac027/03_Bowtie2Benchmark/slurm_out/Test1/C_gOTU_%A.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --partition=short
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yac027@ucsd.edu

# Activate conda environment
source ~/.bashrc
source activate qiime2-2023.2

# Path to script
SCRIPT=/ddn_scratch/yac027/03_Bowtie2Benchmark/scripts/02_post-bt2/C_gOTU.py

echo "Starting gOTU remapping at $(date)"
python $SCRIPT
echo "Finished gOTU remapping at $(date)"

