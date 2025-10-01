#!/bin/bash -l
# AUTHOR: Yang Chen (yac027@ucsd.edu)
# DATE: 10/01/25
# DESCRIPTION: Collapse partitioned BIOM tables (phylum, genus, species) 
#              into unified sample tables across all TBL_DIRS.

#SBATCH --job-name=B_collapse_tbls
#SBATCH --output=/ddn_scratch/yac027/Bowtie2Benchmark/slurm_out/Test1/B_collapse_tbls_%A.out
#SBATCH --error=/ddn_scratch/yac027/Bowtie2Benchmark/slurm_out/Test1/B_collapse_tbls_%A.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --partition=short
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yac027@ucsd.edu

# Activate conda environment
source activate microbiome_base

# Paths
SCRIPT=/ddn_scratch/yac027/Bowtie2Benchmark/scripts/B_collapse_tbls.py

echo "Starting table collapsing at $(date)"
python $SCRIPT
echo "Finished table collapsing at $(date)"

