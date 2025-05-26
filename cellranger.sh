#!/bin/bash
#SBATCH --job-name=cellranger
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00

set -euo pipefail

# Load Environment 
module purge
module load cellranger/3.1.0

# Ensure output directory exists and move there
mkdir -p /path_to_output #your output directory
cd /path_to_output

FASTQ_DIR=""
SAMPLE=""

# Define input files (paired-end)
FQ1="$FASTQ_DIR/${SAMPLE}_1.fastq"
FQ2="$FASTQ_DIR/${SAMPLE}_2.fastq"
FQ1_GZ="$FQ1.gz"
FQ2_GZ="$FQ2.gz"

# Compress if needed
if [ -f "$FQ1" ] && [ ! -f "$FQ1_GZ" ]; then
  echo "Compressing $FQ1"
  gzip "$FQ1"
fi
if [ -f "$FQ2" ] && [ ! -f "$FQ2_GZ" ]; then
  echo "Compressing $FQ2"
  gzip "$FQ2"
fi

# Check for existence
if [ ! -f "$FQ1_GZ" ] || [ ! -f "$FQ2_GZ" ]; then
  echo "Error: Both paired-end files (${SAMPLE}_1.fastq.gz and ${SAMPLE}_2.fastq.gz) are required in $FASTQ_DIR"
  echo "This script only supports paired-end data for CellRanger. Exiting."
  exit 1
fi

# Rename for Cell Ranger compatibility
CR_FQ1="${FASTQ_DIR}/${SAMPLE}_S1_L001_R1_001.fastq.gz"
CR_FQ2="${FASTQ_DIR}/${SAMPLE}_S1_L001_R2_001.fastq.gz"

if [ ! -f "$CR_FQ1" ]; then
  echo "Renaming $FQ1_GZ to $CR_FQ1"
  cp "$FQ1_GZ" "$CR_FQ1"
fi
if [ ! -f "$CR_FQ2" ]; then
  echo "Renaming $FQ2_GZ to $CR_FQ2"
  cp "$FQ2_GZ" "$CR_FQ2"
fi

# Run Cell Ranger
# point your ref directory here
cellranger count \
  --id=$SAMPLE \
  --transcriptome= path_to_ref/GRCh38-2024 \
  --fastqs=$FASTQ_DIR \
  --sample=$SAMPLE \
  --localcores=16 \
  --localmem=64