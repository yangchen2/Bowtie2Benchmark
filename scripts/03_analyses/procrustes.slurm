#!/bin/bash -l
# AUTHOR: Yang Chen (yac027@ucsd.edu)
# DATE: 10/13/25
# DESCRIPTION: Pairwise Procrustes analyses among weighted/unweighted PCoA QZA files from Test1

#SBATCH --job-name=procrustes
#SBATCH --output=../../slurm_out/procrustes_unweighted_%A_%a.out
#SBATCH --error=../../slurm_out/procrustes_unweighted_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=short
#SBATCH --mail-user=yac027@ucsd.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-66

# Define directories
BASE_DIR="../../analyses/beta_div/unweighted_unifrac/Test1"
OUT_DIR="../../analyses/procrustes/unweighted_unifrac"
mkdir -p "$OUT_DIR"

# Collect all pcoa files
mapfile -t PCOA_FILES < <(find "$BASE_DIR" -maxdepth 2 -type f -name "pcoa*.qza" | sort)
N_FILES=${#PCOA_FILES[@]}

echo "Found $N_FILES PCoA files:"
printf '  %s\n' "${PCOA_FILES[@]}"

# Count how many pairs
PAIR_INDEX=0
declare -A PAIRS
for ((i=0; i<N_FILES; i++)); do
  for ((j=i+1; j<N_FILES; j++)); do
    PAIRS[$PAIR_INDEX]="${i},${j}"
    ((PAIR_INDEX++))
  done
done
TOTAL_PAIRS=$PAIR_INDEX

echo "Total pairwise combinations: $TOTAL_PAIRS"

# Run each pair as an array
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
  echo "Running in single-task mode: all pairs will be processed."
  TASK_IDS=$(seq 0 $((TOTAL_PAIRS-1)))
else
  echo "Running as array job: task ID $SLURM_ARRAY_TASK_ID"
  TASK_IDS=$SLURM_ARRAY_TASK_ID
fi

# Main loop for going through pairwise pcoa and running Procrustes
for TASK_ID in $TASK_IDS; do
  pair="${PAIRS[$TASK_ID]}"
  IFS=',' read -r i j <<< "$pair"
  FILE1="${PCOA_FILES[$i]}"
  FILE2="${PCOA_FILES[$j]}"

  NAME1=$(basename "$(dirname "$FILE1")")
  NAME2=$(basename "$(dirname "$FILE2")")

  OUT_NAME="${OUT_DIR}/${NAME1}_vs_${NAME2}_procrustes"

  echo "[$(date)] Running Procrustes: $NAME1 vs $NAME2"

  qiime diversity procrustes-analysis \
    --i-reference "$FILE1" \
    --i-other "$FILE2" \
    --output-dir "$OUT_NAME"

  echo "[$(date)] Finished: $NAME1 vs $NAME2"
done

echo "All Procrustes analyses completed."

