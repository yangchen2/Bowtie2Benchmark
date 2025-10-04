#!/bin/bash -l
#SBATCH --job-name=pcoa-viz
#SBATCH --output=/ddn_scratch/yac027/03_Bowtie2Benchmark/slurm_out/Test1/pcoa-viz_%A_%a.out
#SBATCH --error=/ddn_scratch/yac027/03_Bowtie2Benchmark/slurm_out/Test1/pcoa-viz_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=short
#SBATCH --mail-user=yac027@ucsd.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-23   # 24 total jobs (0 through 23)

set -e  # stop immediately if any command fails

# Activate qiime2 environment
source activate mmvec_qiime2-2020.6_barnacle2

# Path to metadata
metadata=/ddn_scratch/yac027/03_Bowtie2Benchmark/metaG_fastq_skin/metadata.tsv

# Build an array of all 24 directories (12 unw + 12 w)
dirs=(
    $(ls -d /ddn_scratch/yac027/03_Bowtie2Benchmark/analyses/beta_div/unweighted_unifrac/Test1/*/)
    $(ls -d /ddn_scratch/yac027/03_Bowtie2Benchmark/analyses/beta_div/weighted_unifrac/Test1/*/)
)

# Pick the directory for this SLURM task
dir=${dirs[$SLURM_ARRAY_TASK_ID]}
echo "[INFO] Processing directory: $dir"

# Detect which distance matrix exists and set variables
if [ -f "${dir}/species_collapsed_gOTU_unweighted_unifrac.qza" ]; then
    metric="unweighted_unifrac"
    infile="${dir}/species_collapsed_gOTU_unweighted_unifrac.qza"
elif [ -f "${dir}/species_collapsed_gOTU_weighted_unifrac.qza" ]; then
    metric="weighted_unifrac"
    infile="${dir}/species_collapsed_gOTU_weighted_unifrac.qza"
else
    echo "[ERROR] No distance matrix qza found in $dir"
    exit 1
fi

echo "[INFO] Metric detected: $metric"
echo "[INFO] Input file: $infile"

# Output filenames
pcoa_out="${dir}/pcoa_${metric}.qza"
viz_out="${dir}/emperor_${metric}.qzv"

# Step 1: Run PCoA
qiime diversity pcoa \
  --i-distance-matrix "$infile" \
  --o-pcoa "$pcoa_out"
echo "[INFO] Finished PCoA for $dir"

# Step 2: Run Emperor plot
qiime emperor plot \
  --i-pcoa "$pcoa_out" \
  --m-metadata-file "$metadata" \
  --o-visualization "$viz_out"
echo "[INFO] Finished Emperor plot for $dir"

echo "[INFO] Job complete"

