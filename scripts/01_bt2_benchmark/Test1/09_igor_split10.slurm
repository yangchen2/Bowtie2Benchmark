#!/bin/bash -l
# AUTHOR: Yang Chen (yac027@ucsd.edu)
# DATE: 10/01/25
# DESCRIPTION: Igor's Bowtie2 genome-independent read mapping with skin metaG data using WoLr2 10-partition index

#SBATCH --job-name=09_igor_split10
#SBATCH --output=/ddn_scratch/yac027/03_Bowtie2Benchmark/slurm_out/Test1/09_igor_split10_%A_%a.out
#SBATCH --error=/ddn_scratch/yac027/03_Bowtie2Benchmark/slurm_out/Test1/09_igor_split10_%A_%a.err
#SBATCH --array=0-74
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --partition=short
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yac027@ucsd.edu

# Parameters
PARTITION_COUNT=10
PARTITION_DIR="/ddn_scratch/y1weng/71_db_splitting/index10"
FASTQ_DIR="/ddn_scratch/yac027/03_Bowtie2Benchmark/metaG_fastq_skin"
OUTPUT_DIR="/ddn_scratch/yac027/03_Bowtie2Benchmark/sams/Test1/09_igor_split10"
RUNTIME_DIR="/ddn_scratch/yac027/03_Bowtie2Benchmark/runtimes/Test1"
RUNTIME_TSV="${RUNTIME_DIR}/${SLURM_JOB_NAME}.tsv"

mkdir -p "$OUTPUT_DIR" "$RUNTIME_DIR"

# Collect all R1 files
FASTQ_FILES=($FASTQ_DIR/*_R1.fastq.gz)
NUM_FILES=${#FASTQ_FILES[@]}

if [[ $SLURM_ARRAY_TASK_ID -ge $NUM_FILES ]]; then
    echo "SLURM task ID exceeds available FASTQ files. Exiting." >&2
    exit 1
fi

R1=${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}
R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
SAMPLE_NAME=$(basename "$R1" _R1.fastq.gz)

echo "Running Igor Bowtie2 on sample $SAMPLE_NAME"

START=$(date +%s)

cd /ddn_scratch/yac027/Igor_bowtie2_test/bowtie2
for ((i=1; i<=PARTITION_COUNT; i++)); do
    REF_INDEX="${PARTITION_DIR}/split_${i}"
    OUTFILE="${OUTPUT_DIR}/${SAMPLE_NAME}_index${i}.sam"

    echo "Mapping to index partition $i..."

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
      -S "$OUTFILE"

    if [[ $? -ne 0 ]]; then
        echo "Bowtie2 failed on sample $SAMPLE_NAME, partition $i." >&2
    else
        echo "Finished partition $i"
    fi
done

END=$(date +%s)
RUNTIME=$((END - START))
RUNTIME_HMS=$(printf "%02d:%02d:%02d" $((RUNTIME/3600)) $((RUNTIME%3600/60)) $((RUNTIME%60)))

echo "Completed all partitions for $SAMPLE_NAME in ${RUNTIME_HMS}"

# ----------------------------
# Update per-job runtimes TSV
# ----------------------------
if [[ ! -f "$RUNTIME_TSV" ]]; then
    echo -e "Sample\tRuntime" > "$RUNTIME_TSV"
fi

awk -v sample="$SAMPLE_NAME" -v time="$RUNTIME_HMS" '
BEGIN{OFS="\t"}
NR==1{print; next}
{
    if($1==sample){$2=time; found=1}
    print
}
END{
    if(!found) print sample, time
}
' "$RUNTIME_TSV" > "${RUNTIME_TSV}.tmp" && mv "${RUNTIME_TSV}.tmp" "$RUNTIME_TSV"

