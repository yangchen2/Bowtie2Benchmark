#!/bin/bash -l
#SBATCH --job-name=woltka
#SBATCH --output=/ddn_scratch/yac027/Bowtie2Benchmark/slurm_out/Test1/woltka_%A_%a.out
#SBATCH --error=/ddn_scratch/yac027/Bowtie2Benchmark/slurm_out/Test1/woltka_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-user=yac027@ucsd.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-11   # adjust based on number of dirs

# Activate conda environment
source activate mopp

# Define list of SAM directories
SAM_DIRS=(
  "01_k16_no-split"
  "02_a_no-split"
  "03_igor_no-split"
  "04_k16_split4"
  "05_a_split4"
  "06_igor_split4"
  "07_k16_split10"
  "08_a_split10"
  "09_igor_split10"
  "04a_k16_split4a"
  "05a_a_split4a"
  "06a_igor_split4a"
)

# Pick the directory for this array task
SAM_DIR_BASE=${SAM_DIRS[$SLURM_ARRAY_TASK_ID]}

# Extract numeric prefix (everything before the first "_")
JOB_PREFIX=$(echo "$SAM_DIR_BASE" | cut -d'_' -f1)

# Reset job name to prefix number (overrides static #SBATCH name)
JOB_NAME="woltka_${JOB_PREFIX}"
scontrol update JobId=$SLURM_JOB_ID JobName=$JOB_NAME

# Full paths
BASE_IN=/ddn_scratch/yac027/Bowtie2Benchmark/sams/Test1
BASE_OUT=/ddn_scratch/yac027/Bowtie2Benchmark/tables/Test1

SAM_DIR=${BASE_IN}/${SAM_DIR_BASE}
OUTPUT_PATH=${BASE_OUT}/${SAM_DIR_BASE}

# Paths for wolr2
TAXONOMY_DIR=/projects/wol/qiyun/wol2/taxonomy
MAP_FILE=${TAXONOMY_DIR}/taxid.map
NODES_FILE=${TAXONOMY_DIR}/nodes.dmp
NAMES_FILE=${TAXONOMY_DIR}/names.dmp

# Log
echo "[$JOB_PREFIX] Processing: $SAM_DIR"
mkdir -p "$OUTPUT_PATH"

# Run Woltka classify
woltka classify \
  --input "$SAM_DIR" \
  --output "$OUTPUT_PATH" \
  --map "$MAP_FILE" \
  --nodes "$NODES_FILE" \
  --names "$NAMES_FILE" \
  --rank phylum,genus,species

echo "[$JOB_PREFIX] FINISHED: $SAM_DIR"
