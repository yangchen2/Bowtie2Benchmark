#!/bin/bash -l
# AUTHOR: Yang Chen (yac027@ucsd.edu)
# DATE: 10/01/25
# DESCRIPTION: Igor's Bowtie2 genome-independent read mapping of skin metaG data against WoLr2 database (no partitioning)

#SBATCH --job-name=03_igor_no-split
#SBATCH --output=/ddn_scratch/yac027/03_Bowtie2Benchmark/slurm_out/Test1/03_igor_no-split_%A_%a.out
#SBATCH --error=/ddn_scratch/yac027/03_Bowtie2Benchmark/slurm_out/Test1/03_igor_no-split_%A_%a.err
#SBATCH --array=0-74
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --partition=short
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yac027@ucsd.edu

# Activate Igor's Bowtie2 environment
source activate bowtie2

# Set relevant paths
REF_INDEX="/scratch/qp-woltka/WoLr2/WoLr2"
FASTQ_DIR="/ddn_scratch/yac027/03_Bowtie2Benchmark/metaG_fastq_skin"
OUTPUT_DIR="/ddn_scratch/yac027/03_Bowtie2Benchmark/sams/Test1/03_igor_no-split"
RUNTIME_DIR="/ddn_scratch/yac027/03_Bowtie2Benchmark/runtimes/Test1"
RUNTIME_TSV="${RUNTIME_DIR}/${SLURM_JOB_NAME}.tsv"
mkdir -p "$OUTPUT_DIR" "$RUNTIME_DIR"

# Collect FASTQ R1 files
FASTQ_FILES=($FASTQ_DIR/*_R1.fastq.gz)
NUM_FILES=${#FASTQ_FILES[@]}

if [[ $SLURM_ARRAY_TASK_ID -ge $NUM_FILES ]]; then
    echo "SLURM task ID exceeds available FASTQ files. Exiting."
    exit 1
fi

# Assign FASTQ files
FASTQ1=${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}
FASTQ2="${FASTQ1/_R1.fastq.gz/_R2.fastq.gz}"
SAMPLE_NAME=$(basename "$FASTQ1" _R1.fastq.gz)

# Debug info
echo "Processing sample: $SAMPLE_NAME"
echo "FASTQ1: $FASTQ1"
echo "FASTQ2: $FASTQ2"
echo "Using reference index: $REF_INDEX"

# ----------------------------
# Run Igor's Bowtie2 + timing
# ----------------------------
START=$(date +%s)

bowtie2 -p 16 \
  -a -l 100 -d \
  --no-exact-upfront --no-1mm-upfront \
  --very-sensitive \
  --seed 42 \
  --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" \
  --score-min "L,0,-0.05" \
  --no-unal \
  --quiet \
  -x "$REF_INDEX" \
  -1 "$FASTQ1" -2 "$FASTQ2" \
  -S "$OUTPUT_DIR/${SAMPLE_NAME}.sam" \
  2>&1

END=$(date +%s)
RUNTIME=$((END - START))
RUNTIME_HMS=$(printf "%02d:%02d:%02d" $((RUNTIME/3600)) $((RUNTIME%3600/60)) $((RUNTIME%60)))

echo "Finished processing $SAMPLE_NAME in ${RUNTIME_HMS}"

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

