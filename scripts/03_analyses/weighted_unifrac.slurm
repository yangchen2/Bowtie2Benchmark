#!/bin/bash -l
# AUTHOR: Yang Chen (yac027@ucsd.edu)
# DATE: 10/01/25
# DESCRIPTION: Run weighted UniFrac on collapsed_gOTU qza tables across all Test1 directories

#SBATCH --job-name=unifrac_w
#SBATCH --output=/ddn_scratch/yac027/03_Bowtie2Benchmark/slurm_out/Test1/unifrac_w_%x_%j.out
#SBATCH --error=/ddn_scratch/yac027/03_Bowtie2Benchmark/slurm_out/Test1/unifrac_w_%x_%j.err
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

# Base input directory (tables with qza files)
IN_BASE=/ddn_scratch/yac027/03_Bowtie2Benchmark/tables/Test1

# Base output directory
OUT_BASE=/ddn_scratch/yac027/03_Bowtie2Benchmark/analyses/beta_div/weighted_unifrac/Test1

# Loop through each input directory
for dir in ${IN_BASE}/*/; do
    dirname=$(basename "$dir")
    echo "[INFO] Processing directory: $dirname"
    
    # Make output directory for this run
    OUT_DIR="${OUT_BASE}/${dirname}"
    mkdir -p "$OUT_DIR"

    # Loop through each collapsed_gOTU qza table in tables dir
    for qza in "$dir"/*_collapsed_gOTU.qza; do
        [[ -f "$qza" ]] || continue
        base=$(basename "$qza" .qza)

        echo "[INFO] Running weighted UniFrac on $qza"
        qiime diversity beta-phylogenetic \
          --p-metric weighted_unifrac \
          --i-table "$qza" \
          --i-phylogeny "$TREE" \
          --o-distance-matrix "${OUT_DIR}/${base}_weighted_unifrac.qza"

        echo "[INFO] Exporting UniFrac distance matrix for $qza"
        qiime tools export \
          --input-path "${OUT_DIR}/${base}_weighted_unifrac.qza" \
          --output-path "${OUT_DIR}"

        mv "${OUT_DIR}/distance-matrix.tsv" "${OUT_DIR}/${base}_weighted_unifrac.tsv"
    done

    echo "[INFO] Finished directory: $dirname"
done

echo "[DONE] All weighted UniFrac runs completed."

