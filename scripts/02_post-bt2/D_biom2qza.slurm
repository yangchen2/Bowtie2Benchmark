#!/bin/bash -l
#SBATCH --job-name=biom2qza
#SBATCH --output=/ddn_scratch/yac027/03_Bowtie2Benchmark/slurm_out/Test1/D_biom2qza_%A_%a.out
#SBATCH --error=/ddn_scratch/yac027/03_Bowtie2Benchmark/slurm_out/Test1/D_biom2qza_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=short
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yac027@ucsd.edu

# Activate QIIME2 environment
source ~/.bashrc
source activate qiime2-2023.2

# Base directory where your BIOM tables live
BASE=/ddn_scratch/yac027/03_Bowtie2Benchmark/tables/Test1

# Loop through each subdirectory
for dir in $BASE/*; do
    if [ -d "$dir" ]; then
        echo "[INFO] Processing directory: $dir"

        # Loop through each biom file
        for biom in "$dir"/*.biom; do
            [ -e "$biom" ] || continue  # skip if none
            base=$(basename "$biom" .biom)
            out="$dir/${base}.qza"

            echo "[INFO] Converting $biom -> $out"
            qiime tools import \
              --input-path "$biom" \
              --type 'FeatureTable[Frequency]' \
              --input-format BIOMV210Format \
              --output-path "$out"
        done
    fi
done

echo 'Finished!'
