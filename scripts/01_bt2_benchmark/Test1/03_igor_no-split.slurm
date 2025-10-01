#!/bin/bash -l
# AUTHOR: Yang Chen (yac027@ucsd.edu)
# DATE: 10/01/25
# DESCRIPTION: Bowtie2 alignment (Igor's branch, genome-independent mapping with -l 100 -d)
#              of skin metagenome data against WoLr2 database (no partitioning)

#SBATCH --job-name=03_igor_no-split
#SBATCH --output=/ddn_scratch/yac027/Bowtie2Benchmark/slurm_out/Test1/03_igor_no-split_%A_%a.out
#SBATCH --error=/ddn_scratch/yac027/Bowtie2Benchmark/slurm_out/Test1/03_igor_no-split_%A_%a.err
#SBATCH --array=0-74
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --partition=short

# Set relevant paths
REF_INDEX="/scratch/qp-woltka/WoLr2/WoLr2"
FASTQ_DIR="/ddn_scratch/yac027/Bowtie2Benchmark/metaG_fastq_skin"
OUTPUT_DIR="/ddn_scratch/yac027/Bowtie2Benchmark/sams/Test1/03_igor_no-split"
RUNTIME_DIR="/ddn_scratch/yac027/Bowtie2Benchmark/runtimes/Test1"
RUNTIME_TSV="${RUNTIME_DIR}/runtimes.tsv"
mkdir -p "$OUTPUT_DIR" "$RUNTIME_DIR"

# Collect FASTQ R1 files
FASTQ_FILES=($FASTQ_DIR/*_R1.fastq.gz)
NUM_FILES=${#FASTQ_FILES[@]}

if [[ $SLURM_ARRAY_TASK_ID -ge $NUM_FILES ]]; then
    echo "SLURM task ID exceeds available FASTQ files. Exiting."
    exit 1
fi

# Assign FASTQ files
R1=${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}
R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
SAMPLE_NAME=$(basename "$R1" _R1.fastq.gz)

# Debug info
echo "Processing sample: $SAMPLE_NAME"
echo "R1: $R1"
echo "R2: $R2"
echo "Using reference index: $REF_INDEX"

# ----------------------------
# Run Bowtie2 (Igor branch) + timing
# ----------------------------
START=$(date +%s)

cd /ddn_scratch/yac027/Igor_bowtie2_test/bowtie2
./bowtie2 -p 16 \
  --no-exact-upfront --no-1mm-upfront \
  -x "$REF_INDEX" \
  -1 "$R1" -2 "$R2" \
  --very-sensitive \
  --seed 42 \
  -a \
  --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" \
  --score-min "L,0,-0.05" \
  --no-unal \
  -l 100 -d \
  --quiet \
  -S "$OUTPUT_DIR/${SAMPLE_NAME}.sam"

END=$(date +%s)
RUNTIME=$((END - START))
RUNTIME_HMS=$(printf "%02d:%02d:%02d" $((RUNTIME/3600)) $((RUNTIME%3600/60)) $((RUNTIME%60)))

echo "Finished processing $SAMPLE_NAME in ${RUNTIME_HMS}"

# ----------------------------
# Update runtimes.tsv
# ----------------------------
if [[ ! -f "$RUNTIME_TSV" ]]; then
    echo -e "Sample\t${SLURM_JOB_NAME}" > "$RUNTIME_TSV"
fi

if ! head -n1 "$RUNTIME_TSV" | grep -qw "${SLURM_JOB_NAME}"; then
    sed -i "1s/$/\t${SLURM_JOB_NAME}/" "$RUNTIME_TSV"
    awk 'NR==1{print;next}{print $0 "\t"}' "$RUNTIME_TSV" > "${RUNTIME_TSV}.tmp" && mv "${RUNTIME_TSV}.tmp" "$RUNTIME_TSV"
fi

if grep -q "^${SAMPLE_NAME}" "$RUNTIME_TSV"; then
    awk -v sample="$SAMPLE_NAME" -v job="${SLURM_JOB_NAME}" -v time="$RUNTIME_HMS" '
        BEGIN{OFS="\t"}
        NR==1{for(i=1;i<=NF;i++) if($i==job) col=i}
        NR>1{if($1==sample){$col=time}}
        {print}
    ' "$RUNTIME_TSV" > "${RUNTIME_TSV}.tmp" && mv "${RUNTIME_TSV}.tmp" "$RUNTIME_TSV"
else
    COLS=$(head -n1 "$RUNTIME_TSV" | awk '{print NF}')
    NEWROW="$SAMPLE_NAME"
    for ((i=2; i<=COLS; i++)); do
        if [[ $(head -n1 "$RUNTIME_TSV" | cut -f$i) == "${SLURM_JOB_NAME}" ]]; then
            NEWROW="${NEWROW}\t${RUNTIME_HMS}"
        else
            NEWROW="${NEWROW}\t"
        fi
    done
    echo -e "$NEWROW" >> "$RUNTIME_TSV"
fi

