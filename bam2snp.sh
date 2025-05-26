#!/bin/bash
#SBATCH --job-name=cell_snp_pipeline
#SBATCH --output=snv_result/%x_%j.out
#SBATCH --error=snv_result/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00

set -euo pipefail

# Load Environment 
module purge
# Change this to your miniconda path if needed
source ~/miniconda3/etc/profile.d/conda.sh
conda activate cellsnp_env   # Or your own env with pysam, samtools, cellsnp-lite

# === USER CONFIG ===
SAMPLE="$1"
BAMDIR="dir_to_your_path/${SAMPLE}/outs"
INPUT_BAM="${BAMDIR}/possorted_genome_bam.bam"
PER_CELL_DIR="${BAMDIR}/bampercell"
SNVOUT="dir_to_your_output"
SNPVCF="dir_to_your_reference"
BARCODE_TAG="CB"      # "CB" for 10x, "XC" for Drop-seq
MIN_READS=100 #threshold

mkdir -p "$PER_CELL_DIR"
mkdir -p "$SNVOUT"

#1.Split BAM into per-cell BAMs 
echo "Splitting BAM using splitbam.py..."
#make sure you point the right directory here
python splitbam.py \ 
    "$INPUT_BAM" "$PER_CELL_DIR" "$BARCODE_TAG" "$MIN_READS"

#1.2: Index all per-cell BAMs
echo "Indexing all per-cell BAMs..."
find "$PER_CELL_DIR" -name '*.bam' | while read bamfile; do
    samtools index "$bamfile"
done

#1.3: Prepare bam_list.txt and sample_list.txt
echo "Generating input lists for cellsnp-lite..."
find "$PER_CELL_DIR" -name '*.bam' | sort > "${SNVOUT}/${SAMPLE}_bam_list.txt"
awk -F'/' '{print $NF}' "${SNVOUT}/${SAMPLE}_bam_list.txt" | sed 's/.bam$//' > "${SNVOUT}/${SAMPLE}_sample_list.txt"

#2: Run cellsnp-lite
echo "Running cellsnp-lite..."
cellsnp-lite \
    -S "${SNVOUT}/${SAMPLE}_bam_list.txt" \
    -i "${SNVOUT}/${SAMPLE}_sample_list.txt" \
    -O "${SNVOUT}/cellsnp_${SAMPLE}" \
    -R "$SNPVCF" \
    -p 8 \
    --cellTAG None \
    --UMItag None \
    --gzip

# === Step 3: Cleanup ===
echo "Cleaning up intermediate files..."

# Remove input lists
rm -f "${SNVOUT}/${SAMPLE}_bam_list.txt"
rm -f "${SNVOUT}/${SAMPLE}_sample_list.txt"

# Remove per-cell BAMs (optional – **WARNING: only do this if you’re 100% done!**)
# rm -rf "$PER_CELL_DIR"

echo "Cleanup complete!"
echo "Pipeline completed for $SAMPLE!"