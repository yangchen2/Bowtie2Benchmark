#!/bin/bash -l
# AUTHOR: Yang Chen (yac027@ucsd.edu)
# DATE: 10/14/25
# DESCRIPTION: Aggregate Procrustes M² and p-values for weighted and unweighted UniFrac analyses.

#SBATCH --job-name=agg_procrustes_M2pval
#SBATCH --output=../../slurm_out/agg_procrustes_M2pval_%j.out
#SBATCH --error=../../slurm_out/agg_procrustes_M2pval_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=short
#SBATCH --mail-user=yac027@ucsd.edu
#SBATCH --mail-type=END,FAIL

# === Paths ===
BASE_PROCRUSTES="../../analyses/procrustes"
OUT_DIR="/ddn_scratch/yac027/03_Bowtie2Benchmark/analyses/procrustes"
OUT_FILE="${OUT_DIR}/procrustes_all_M2_pval.tsv"
TMP_DIR=$(mktemp -d)

echo "[$(date)] Starting aggregation of Procrustes M² and p-values..."
mkdir -p "$OUT_DIR"

################################################################################
# STEP 1: Extract M² and p-value for each metric
################################################################################

for METRIC in unweighted weighted; do
    BASE="${BASE_PROCRUSTES}/${METRIC}_unifrac"
    TMP_FILE="${TMP_DIR}/${METRIC}_M2pval.tsv"

    echo -e "comparison\t${METRIC}_unf_procrustes_M2\t${METRIC}_unf_procrustes_pval" > "$TMP_FILE"

    echo "Processing $METRIC UniFrac..."
    for d in "${BASE}"/*/; do
        TSV="${d}/disparity_exported/ProcrustesStatistics.tsv"
        if [[ -f "$TSV" ]]; then
            # extract true M² and p-value columns from 'results' row
            M2=$(awk '/^results/{print $2}' "$TSV")
            PVAL=$(awk '/^results/{print $3}' "$TSV")
            echo -e "$(basename "$d")\t${M2}\t${PVAL}" >> "$TMP_FILE"
        else
            echo -e "$(basename "$d")\tNA\tNA" >> "$TMP_FILE"
        fi
    done
    echo "[$(date)] Finished $METRIC UniFrac extraction."
done

################################################################################
# STEP 2: Merge weighted + unweighted tables
################################################################################

UNW="${TMP_DIR}/unweighted_M2pval.tsv"
W="${TMP_DIR}/weighted_M2pval.tsv"

echo "Merging weighted and unweighted tables..."
join -t $'\t' -a1 -a2 -e "NA" -o auto \
  <(sort "$UNW") <(sort "$W") > "$OUT_FILE"

################################################################################
# STEP 3: Ensure header is at the top
################################################################################

HEADER="comparison\tunweighted_unf_procrustes_M2\tunweighted_unf_procrustes_pval\tweighted_unf_procrustes_M2\tweighted_unf_procrustes_pval"
{ echo -e "$HEADER"; grep -v '^comparison' "$OUT_FILE"; } > "${OUT_FILE}.tmp" && mv "${OUT_FILE}.tmp" "$OUT_FILE"

echo "[$(date)] Final file written to: $OUT_FILE"
rm -rf "$TMP_DIR"
echo "[$(date)] Temporary files removed."

