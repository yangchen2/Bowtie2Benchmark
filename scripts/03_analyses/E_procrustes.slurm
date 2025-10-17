#!/usr/bin/env bash
# AUTHOR: Yang Chen (yac027@ucsd.edu)
# DATE: 10/16/25
# DESCRIPTION: Run pairwise Procrustes analyses and export disparity results
#              for both weighted and unweighted UniFrac data (15 groups = 105 pairs).

#SBATCH --job-name=procrustes
#SBATCH --output=../../slurm_out/E_procrustes_%A_%a.out
#SBATCH --error=../../slurm_out/E_procrustes_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=short
#SBATCH --mail-user=yac027@ucsd.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-209   # 0–104 = unweighted; 105–209 = weighted (total 210 tasks)

set -euo pipefail
source ~/.bashrc
conda activate qiime2-2023.2

################################################################################
# SECTION 0: Determine which metric this array job belongs to
################################################################################

if (( SLURM_ARRAY_TASK_ID < 105 )); then
    METRIC="unweighted_unifrac"
    METRIC_ID=$SLURM_ARRAY_TASK_ID
else
    METRIC="weighted_unifrac"
    METRIC_ID=$((SLURM_ARRAY_TASK_ID - 105))
fi

echo "[$(date)] Running for metric: $METRIC (Array index: $METRIC_ID)"

################################################################################
# SECTION 1: Pairwise Procrustes Analyses
################################################################################

BASE_DIR="../../analyses/beta_div/${METRIC}/Test1"
OUT_DIR="../../analyses/procrustes/${METRIC}"
mkdir -p "$OUT_DIR"

# Collect all PCoA QZA files (should be 15 total)
mapfile -t PCOA_FILES < <(find "$BASE_DIR" -type f -iname "pcoa*.qza" | sort -V)
N_FILES=${#PCOA_FILES[@]}

echo "Found $N_FILES PCoA files for $METRIC:"
printf '  %s\n' "${PCOA_FILES[@]}"

TOTAL_PAIRS=105
echo "Assuming total pairwise combinations for $METRIC: $TOTAL_PAIRS"

if (( METRIC_ID >= TOTAL_PAIRS )); then
    echo "Task ID $METRIC_ID exceeds total pairs ($TOTAL_PAIRS). Exiting."
    exit 0
fi

# Generate pair indices explicitly (15 choose 2 = 105)
PAIRS=()
for ((i=0; i<N_FILES; i++)); do
  for ((j=i+1; j<N_FILES; j++)); do
    PAIRS+=("${i},${j}")
  done
done

pair="${PAIRS[$METRIC_ID]}"
IFS=',' read -r i j <<< "$pair"
FILE1="${PCOA_FILES[$i]}"
FILE2="${PCOA_FILES[$j]}"

NAME1=$(basename "$(dirname "$FILE1")")
NAME2=$(basename "$(dirname "$FILE2")")
OUT_NAME="${OUT_DIR}/${NAME1}_vs_${NAME2}_procrustes"

# Skip if output directory already exists
if [[ -d "$OUT_NAME" ]]; then
    echo "[$(date)] Skipping $NAME1 vs $NAME2 for $METRIC (already exists: $OUT_NAME)"
else
    echo "[$(date)] Running Procrustes: $NAME1 vs $NAME2 for $METRIC"

    if [[ ! -f "$FILE1" || ! -f "$FILE2" ]]; then
        echo "Missing PCoA input(s): $FILE1 or $FILE2. Skipping."
        exit 0
    fi

    qiime diversity procrustes-analysis \
      --i-reference "$FILE1" \
      --i-other "$FILE2" \
      --output-dir "$OUT_NAME" || {
        echo "[$(date)] Procrustes failed for $NAME1 vs $NAME2 ($METRIC)"
        exit 1
      }

    echo "[$(date)] Finished: $NAME1 vs $NAME2 for $METRIC"
fi

################################################################################
# SECTION 2: Export disparity_results.qza for this metric
################################################################################

PROCRUSTES_DIR="../../analyses/procrustes/${METRIC}"

mapfile -t DIRS < <(find "$PROCRUSTES_DIR" -mindepth 1 -maxdepth 1 -type d | sort)
NUM_DIRS=${#DIRS[@]}
echo "Found $NUM_DIRS Procrustes result directories for $METRIC."

if (( METRIC_ID < NUM_DIRS )); then
    DIR="${DIRS[$METRIC_ID]}"
    echo "[$(date)] Exporting disparity results for: $DIR"

    QZA_PATH="${DIR}/disparity_results.qza"
    OUT_EXP="${DIR}/disparity_exported"

    if [[ -f "$QZA_PATH" ]]; then
        mkdir -p "$OUT_EXP"
        qiime tools export \
            --input-path "$QZA_PATH" \
            --output-path "$OUT_EXP"
        echo "[$(date)] Exported disparity_results.qza for $DIR"
    else
        echo "No disparity_results.qza found in $DIR — skipping export."
    fi
else
    echo "Metric index $METRIC_ID exceeds directory count ($NUM_DIRS); skipping export."
fi

echo "[$(date)] Completed processing for $METRIC (Task $SLURM_ARRAY_TASK_ID)"

