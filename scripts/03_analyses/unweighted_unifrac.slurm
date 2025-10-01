#!/bin/bash -l
#SBATCH --job-name=unifrac_unw
#SBATCH --output=/ddn_scratch/yac027/Bowtie2Benchmark/slurm_out/Test1/unifrac_unw_%x_%j.out
#SBATCH --error=/ddn_scratch/yac027/Bowtie2Benchmark/slurm_out/Test1/unifrac_unw_%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=short
#SBATCH --mail-user=yac027@ucsd.edu
#SBATCH --mail-type=END,FAIL

# Activate QIIME2 environment
source activate mmvec_qiime2-2020.6_barnacle2

# Path to reference tree
TREE=/projects/wol/qiyun/wol2/databases/qiime2/phylogeny.qza

# Base output directory
OUT_BASE=/ddn_scratch/yac027/Bowtie2Benchmark/analyses/beta_div/unweighted_unifrac/Test1

# Loop through each input directory
for dir in /ddn_scratch/yac027/Bowtie2Benchmark/tables/Test1/*/; do
    dirname=$(basename "$dir")
    echo "Processing directory: $dirname"
    
    # Make output directory for this run
    OUT_DIR="${OUT_BASE}/${dirname}"
    mkdir -p "$OUT_DIR"
    cd "$dir" || continue

    # Loop through each collapsed_gOTU biom table
    for biom in *_collapsed_gOTU.biom; do
        [[ -f "$biom" ]] || continue
        base=$(basename "$biom" .biom)

        echo "Converting $biom to QZA"
        qiime tools import \
          --input-path "$biom" \
          --type 'FeatureTable[Frequency]' \
          --input-format BIOMV210Format \
          --output-path "${OUT_DIR}/${base}.qza"

        echo "Running unweighted UniFrac on $biom"
        qiime diversity beta-phylogenetic \
          --p-metric unweighted_unifrac \
          --i-table "${OUT_DIR}/${base}.qza" \
          --i-phylogeny "$TREE" \
          --o-distance-matrix "${OUT_DIR}/${base}_unweighted_unifrac.qza"

        echo "Exporting UniFrac distance matrix for $biom"
        qiime tools export \
          --input-path "${OUT_DIR}/${base}_unweighted_unifrac.qza" \
          --output-path "${OUT_DIR}"

        mv "${OUT_DIR}/distance-matrix.tsv" "${OUT_DIR}/${base}_unweighted_unifrac.tsv"
    done

    echo "Finished directory: $dirname"
done

