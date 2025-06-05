#!/bin/bash
##
## scRNA-seq_Variant_Calling_Pipeline.sh
##
## This script automates a single-cell RNA sequencing analysis workflow,
## from raw FASTQ files through alignment to variant calling.
## It supports both 10x Genomics (via Cell Ranger) and Drop-Seq/SMART-Seq2
## (via STARsolo) protocols, followed by quality control with FastQC and
## variant calling with cellsnp-lite.
##
## Usage: sbatch scRNA-seq_Variant_Calling_Pipeline.sh <SAMPLE_ID> <PROTOCOL>
##
## <SAMPLE_ID>: A unique identifier for your sample (e.g., SRR12298114).
##              This ID will be used for naming output directories and files.
## <PROTOCOL>: The sequencing protocol used. Case-insensitive.
##             - "10x"       : For 10x Genomics data (uses Cell Ranger).
##             - "Drop-Seq"  : For Drop-Seq data (uses STARsolo).
##             - "SMART-Seq2": For SMART-Seq2 data (uses STARsolo).
##             - "Plate"     : Alias for SMART-Seq2/plate-based data (uses STARsolo).
##
## Prerequisites:
##   - SLURM Workload Manager: This script is designed to be submitted as a Slurm job.
##   - Software Modules: The following modules must be available on your system:
##     - cellranger/<VERSION> (e.g., cellranger/3.1.0)
##     - star/<VERSION>     (e.g., star/2.7.11b)
##     - fastqc/<VERSION>   (e.g., fastqc/0.12.1)
##     - samtools/<VERSION> (e.g., samtools/1.10)
##   - Conda Environment: A Conda environment named 'cellsnp_env' must exist and
##     contain 'cellsnp-lite'.
##   - Reference Data: All specified reference genome, transcriptome, SNP VCF,
##     and whitelist files must be correctly indexed (e.g., FASTA with .fai).
##   - FASTQ File Naming:
##     - For '10x' protocol: FASTQ files should be in the 10x-compatible format
##       (e.g., SAMPLE_ID_S1_L001_R1_001.fastq.gz, SAMPLE_ID_S1_L001_R2_001.fastq.gz,
##       and SAMPLE_ID_S1_L001_I1_001.fastq.gz for the index read).
##     - For 'Drop-Seq'/'SMART-Seq2'/'Plate' protocols: FASTQ files should be
##       named as SAMPLE_ID_1.fastq.gz and SAMPLE_ID_2.fastq.gz (or .fq.gz/.fastq).
##
## Recommended SLURM Resources:
##   - --cpus-per-task: 16 (for STARsolo/Cell Ranger multithreading)
##   - --mem: 64G (adjust based on genome size and read depth, for human data)
##
## Output Directory Structure:
##   - ${BASE_PROJECT_DIR}/logs/              : Slurm job output/error logs.
##   - ${BASE_PROJECT_DIR}/bam_data/         : Contains alignment outputs.
##     - fastqc_reports/${SAMPLE_ID}/ : FastQC reports.
##     - ${SAMPLE_ID}_analysis/       : Cell Ranger output directory (for 10x).
##     - ${SAMPLE_ID}/                : STARsolo output directory (for Drop-Seq/SMART-Seq2).
##   - ${BASE_PROJECT_DIR}/snv_result/       : Contains variant calling results.
##     - cellsnp_${SAMPLE_ID}/        : cellsnp-lite output files.
##
##

# Request resources for the job using SLURM directives
#SBATCH --job-name=scRNA_variant_calling  # Name of your job
#SBATCH --cpus-per-task=16              # Number of CPU cores per task
#SBATCH --mem=64G                       # Amount of RAM (e.g., 64GB)
#SBATCH --time=120:00:00                # Maximum job run time (HH:MM:SS)
#SBATCH --output=/%x_%j.out # Standard output file
#SBATCH --error=/%x_%j.err   # Standard error file
#SBATCH --mail-user= # Email address for notifications
#SBATCH --mail-type=BEGIN,END,FAIL      # Email notifications for job start, end, and failure

# Exit immediately if a command exits with a non-zero status.
# Treat unset variables as an error when substituting.
# Exit if any command in a pipeline fails.
set -euo pipefail

# === USER CONFIGURATION (EDIT THESE PATHS AND OPTIONS) ===
# This section defines all key directories, versions, and parameters
# that you might need to modify for your specific setup or data.

# Mandatory arguments (passed when submitting the job)
SAMPLE_ID="$1"       # Example: SRR12298114
PROTOCOL="$2"        # Example: 10x, Drop-Seq, SMART-Seq2, Plate

# Base Project Directory: All other paths are relative to this.
BASE_PROJECT_DIR=""

# FASTQ Directory Configuration:
# Define the base directory for your raw FASTQ files.
# Then, specify the exact subdirectories where Cell Ranger or STARsolo
# should look for the FASTQ files based on the protocol.
FASTQ_BASE_DIR="${BASE_PROJECT_DIR}/"
# For 10x Genomics (Cell Ranger): FASTQs are expected directly in this directory
# (e.g., ${FASTQ_DIR_10X}/SRR12298114_S1_L001_R1_001.fastq.gz).
FASTQ_DIR_10X="${FASTQ_BASE_DIR}/"
# For Drop-Seq/SMART-Seq2/Plate (STARsolo): FASTQs are expected in sample-specific
# subdirectories within this base (e.g., ${FASTQ_DIR_DROPSEQ_BASE}/SRR12298114/SRR12298114_1.fastq.gz).
FASTQ_DIR_DROPSEQ_BASE="${FASTQ_BASE_DIR}/"

# Software Version Configuration:
# Ensure these versions match the modules available on your HPC system.
CELLRANGER_VERSION="3.1.0" # Version of Cell Ranger module to load
STAR_VERSION="2.7.11b"     # Version of STAR module to load
FASTQC_VERSION="0.12.1"    # Version of FastQC module to load
SAMTOOLS_VERSION="1.10"    # Version of Samtools module to load

# Reference Data Directories and Files:
# IMPORTANT: Verify these paths and ensure all reference files are correctly indexed.
REF_DIR="${BASE_PROJECT_DIR}/ref"
# For Cell Ranger (10x): Path to the 10x Genomics-formatted reference transcriptome.
CELLRANGER_REF_TRANSCRIPTOME="${REF_DIR}/GRCh38-2024"
# For STARsolo: Path to the STAR-indexed genome directory.
STAR_GENOME_DIR="${REF_DIR}/GRCh38_Ensembl_111/star_index"
# For STARsolo 10x protocol: Path to the 10x Cell Barcode whitelist file.
TENX_WHITELIST="${REF_DIR}/3M-february-2018.txt"

# Reference Files for cellsnp-lite variant calling:
SNP_VCF="${REF_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf.gz" # Human SNP VCF file
REFERENCE_FASTA="${REF_DIR}/Homo_sapiens_assembly38.fasta"   # Human reference FASTA. Must have .fai index!

# Conda environment configuration for cellsnp-lite:
# Base path to your Miniconda/Anaconda installation.
CONDA_BASE="~/miniconda3"
# Name of the Conda environment containing cellsnp-lite.
CELLSNP_CONDA_ENV="cellsnp_env"

# Output Directories (these will be created if they don't exist):
BAM_DATA_DIR="${BASE_PROJECT_DIR}/bam_data"       # Directory for alignment outputs (BAMs, Cell Ranger/STARsolo results)
SNV_RESULT_DIR="${BASE_PROJECT_DIR}/snv_result"   # Directory for variant calling results
LOG_DIR="${BASE_PROJECT_DIR}/logs"                # Centralized directory for Slurm job output/error logs

# === DO NOT EDIT BELOW THIS LINE UNLESS YOU KNOW WHAT YOU'RE DOING ===
# This section contains the core logic of the pipeline.

echo "[$(date)] Starting scRNA-seq pipeline for sample: ${SAMPLE_ID}, protocol: ${PROTOCOL}."
echo "Job running on host: $(hostname)"
echo "Base project directory: ${BASE_PROJECT_DIR}"

# Create primary output directories if they don't exist
mkdir -p "$BAM_DATA_DIR"
mkdir -p "$SNV_RESULT_DIR"
mkdir -p "$LOG_DIR"

# Variables to be set dynamically based on protocol
FASTQ_INPUT_DIR=""      # The specific directory where FASTQs for the current sample are located
ALIGNMENT_OUT_DIR=""    # The output directory for the alignment step (Cell Ranger or STARsolo)
QC_OUT_DIR="${BAM_DATA_DIR}/fastqc_reports/${SAMPLE_ID}" # Dedicated FastQC output directory per sample

# Create the FastQC output directory
mkdir -p "$QC_OUT_DIR"

# --- Protocol-specific FASTQ input and alignment output directory assignment ---
# This block sets the correct FASTQ_INPUT_DIR and ALIGNMENT_OUT_DIR based on the chosen protocol.
if [[ "$PROTOCOL" =~ [Tt]en[Xx]|[1]0[xX] ]]; then
    echo "[$(date)] Protocol: 10x Genomics (will use Cell Ranger)."
    FASTQ_INPUT_DIR="${FASTQ_DIR_10X}"
    # Cell Ranger creates its own sample-specific subdirectory (e.g., bam_data/SRR12298114_analysis)
    ALIGNMENT_OUT_DIR="${BAM_DATA_DIR}"
elif [[ "$PROTOCOL" =~ [Dd]rop.?[Ss]eq|[Ss][Mm][Aa][Rr][Tt]|[Pp]late ]]; then
    echo "[$(date)] Protocol: STARsolo-compatible (${PROTOCOL})."
    # STARsolo expects FASTQs in a sample-specific subdirectory
    FASTQ_INPUT_DIR="${FASTQ_DIR_DROPSEQ_BASE}/${SAMPLE_ID}"
    # STARsolo will output directly into a sample-specific subdirectory within bam_data
    ALIGNMENT_OUT_DIR="${BAM_DATA_DIR}/${SAMPLE_ID}"
else
    echo "ERROR: Protocol '${PROTOCOL}' not recognized. Supported: '10x', 'Drop-Seq', 'SMART-Seq2', 'Plate'."
    exit 1
fi

echo "FASTQ input directory: ${FASTQ_INPUT_DIR}"
echo "Alignment output directory: ${ALIGNMENT_OUT_DIR}"

# --- Universal FASTQ file detection (R1/R2) for FastQC and STARsolo ---
# This section tries to find the paired-end FASTQ files (R1 and R2) for the sample.
# It checks for common extensions (.fastq.gz, .fq.gz, .fastq).
# Note: For the '10x' protocol, an additional check for the I1 (Index) file
# will be performed specifically in the Cell Ranger section below.
R1=""
R2=""

# Check for .fastq.gz files
if [[ -e "${FASTQ_INPUT_DIR}/${SAMPLE_ID}_1.fastq.gz" && -e "${FASTQ_INPUT_DIR}/${SAMPLE_ID}_2.fastq.gz" ]]; then
    R1="${FASTQ_INPUT_DIR}/${SAMPLE_ID}_1.fastq.gz"
    R2="${FASTQ_INPUT_DIR}/${SAMPLE_ID}_2.fastq.gz"
# Check for .fq.gz files
elif [[ -e "${FASTQ_INPUT_DIR}/${SAMPLE_ID}_1.fq.gz" && -e "${FASTQ_INPUT_DIR}/${SAMPLE_ID}_2.fq.gz" ]]; then
    R1="${FASTQ_INPUT_DIR}/${SAMPLE_ID}_1.fq.gz"
    R2="${FASTQ_INPUT_DIR}/${SAMPLE_ID}_2.fq.gz"
# Check for uncompressed .fastq files
elif [[ -e "${FASTQ_INPUT_DIR}/${SAMPLE_ID}_1.fastq" && -e "${FASTQ_INPUT_DIR}/${SAMPLE_ID}_2.fastq" ]]; then
    R1="${FASTQ_INPUT_DIR}/${SAMPLE_ID}_1.fastq"
    R2="${FASTQ_INPUT_DIR}/${SAMPLE_ID}_2.fastq"
else
    echo "ERROR: Paired-end FASTQ files (e.g., ${SAMPLE_ID}_1.fastq.gz) not found for ${SAMPLE_ID} in ${FASTQ_INPUT_DIR}"
    exit 1
fi

echo "Detected R1 FASTQ file: ${R1}"
echo "Detected R2 FASTQ file: ${R2}"

---
## Step 1: Quality Control with FastQC

echo "[$(date)] Starting FastQC for ${SAMPLE_ID}..."
# Load the FastQC module
module load fastqc/"$FASTQC_VERSION"

# Run FastQC on both read files
fastqc "$R1" "$R2" --outdir "$QC_OUT_DIR"

echo "[$(date)] FastQC finished. Reports are located in: ${QC_OUT_DIR}"
# Unload the FastQC module to keep the environment clean
module unload fastqc/"$FASTQC_VERSION"

---
## Step 2: FASTQ to BAM Alignment (Cell Ranger or STARsolo)

# Initialize variables for the final BAM and barcode file paths
# These paths will be determined by the chosen alignment tool.
INPUT_BAM=""                # Path to the final aligned BAM file
BARCODE_FILE_UNCOMPRESSED="" # Path to the uncompressed cell barcode file

if [[ "$PROTOCOL" =~ [Tt]en[Xx]|[1]0[xX] ]]; then
    echo "[$(date)] Running Cell Ranger 'count' for ${SAMPLE_ID}..."
    # Load the Cell Ranger module
    module load cellranger/"$CELLRANGER_VERSION"
    echo "[$(date)] Cell Ranger executable path: $(which cellranger)" # Verify module load

    # Change to the BAM_DATA_DIR, as Cell Ranger will create its output
    # directory (${SAMPLE_ID}_analysis) directly within it.
    mkdir -p "${BAM_DATA_DIR}"
    cd "${BAM_DATA_DIR}"

    # --- Cell Ranger Specific FASTQ Verification (including I1 read) ---
    # Cell Ranger expects FASTQ files to be in a specific 10x Genomics format,
    # typically generated by bcl2fastq or similar tools. This includes the
    # R1 (Read 1), R2 (Read 2), and I1 (Index 1) files, named with the
    # '_S1_L001_R1_001.fastq.gz' pattern.
    EXPECTED_CR_R1="${FASTQ_INPUT_DIR}/${SAMPLE_ID}_S1_L001_R1_001.fastq.gz"
    EXPECTED_CR_R2="${FASTQ_INPUT_DIR}/${SAMPLE_ID}_S1_L001_R2_001.fastq.gz"
    EXPECTED_CR_I1="${FASTQ_INPUT_DIR}/${SAMPLE_ID}_S1_L001_I1_001.fastq.gz" # Crucial for 10x!

    echo "[$(date)] Checking for Cell Ranger-ready FASTQ files for sample: ${SAMPLE_ID} in ${FASTQ_INPUT_DIR}..."
    if [ ! -f "$EXPECTED_CR_R1" ]; then
        echo "Error: Cell Ranger R1 file not found: ${EXPECTED_CR_R1}"
        echo "Please ensure FASTQ files are in 10x-compatible format (e.g., generated by bcl2fastq or similar prep)."
        exit 1
    fi
    if [ ! -f "$EXPECTED_CR_R2" ]; then
        echo "Error: Cell Ranger R2 file not found: ${EXPECTED_CR_R2}"
        echo "Please ensure FASTQ files are in 10x-compatible format (e.g., generated by bcl2fastq or similar prep)."
        exit 1
    fi
    if [ ! -f "$EXPECTED_CR_I1" ]; then
        echo "Error: Cell Ranger I1 (Index) file not found: ${EXPECTED_CR_I1}"
        echo "This file is crucial for Cell Ranger's 10x Genomics analysis."
        echo "Please ensure FASTQ files are in 10x-compatible format (e.g., generated by bcl2fastq or similar prep)."
        exit 1
    fi
    echo "[$(date)] All required Cell Ranger FASTQ files (R1, R2, I1) found."

    # Execute Cell Ranger 'count' command
    # --id: The name of the output directory (e.g., SRR12298114_analysis)
    # --transcriptome: Path to the 10x Genomics reference transcriptome
    # --fastqs: Directory containing the Cell Ranger-formatted FASTQ files
    # --sample: The sample prefix used in the FASTQ filenames
    # --localcores: Number of cores to use (from Slurm allocation)
    # --localmem: Amount of memory to use (from Slurm allocation)
    cellranger count \
        --id="${SAMPLE_ID}_analysis" \
        --transcriptome="$CELLRANGER_REF_TRANSCRIPTOME" \
        --fastqs="$FASTQ_INPUT_DIR" \
        --sample="$SAMPLE_ID" \
        --localcores="$SLURM_CPUS_PER_TASK" \
        --localmem="${SLURM_MEM}" # SLURM_MEM includes 'G' or 'M' suffix

    echo "[$(date)] Cell Ranger pipeline finished for ${SAMPLE_ID}."

    # Define paths to Cell Ranger outputs needed for cellsnp-lite
    INPUT_BAM="${BAM_DATA_DIR}/${SAMPLE_ID}_analysis/outs/possorted_genome_bam.bam"
    BARCODE_FILE_UNCOMPRESSED="${BAM_DATA_DIR}/${SAMPLE_ID}_analysis/outs/filtered_gene_bc_matrices/GRCh38-2024/barcodes.tsv"

    # Unload the Cell Ranger module
    module unload cellranger/"$CELLRANGER_VERSION"

elif [[ "$PROTOCOL" =~ [Dd]rop.?[Ss]eq|[Ss][Mm][Aa][Rr][Tt]|[Pp]late ]]; then
    echo "[$(date)] Running STARsolo for ${SAMPLE_ID} with ${PROTOCOL} protocol..."
    # Load the STAR module
    module load star/"$STAR_VERSION"

    # Create the sample-specific output directory for STARsolo and change into it.
    # STARsolo will place its output files directly here.
    mkdir -p "$ALIGNMENT_OUT_DIR"
    cd "$ALIGNMENT_OUT_DIR"

    # Base STARsolo command components
    STAR_COMMAND="STAR --genomeDir \"$STAR_GENOME_DIR\" \\
         --runThreadN \"$SLURM_CPUS_PER_TASK\" \\
         --readFilesIn \"$R2\" \"$R1\" \\
         --readFilesCommand zcat \\
         --outSAMtype BAM SortedByCoordinate \\
         --outSAMattributes NH HI AS nM CB UB GN GX"

    # Add protocol-specific parameters for STARsolo
    if [[ "$PROTOCOL" =~ [Tt]en[Xx]|[1]0[xX] ]]; then
        # This block defines STARsolo parameters for 10x Genomics chemistry.
        # Note: If PROTOCOL is "10x", the script prioritizes Cell Ranger above.
        # This section would be used if you explicitly used "10x_STARsolo" or similar.
        local CB_START=1; local CB_LEN=16; local UMI_START=17; local UMI_LEN=12
        STAR_COMMAND+=" \\
         --soloType CB_UMI_Simple \\
         --soloCBstart $CB_START --soloCBlen $CB_LEN \\
         --soloUMIstart $UMI_START --soloUMIlen $UMI_LEN \\
         --soloCBwhitelist \"$TENX_WHITELIST\" \\
         --soloFeatures Gene \\
         --soloCellFilter EmptyDrops_CR \\
         --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \\
         --soloUMIfiltering MultiGeneUMI_CR \\
         --soloUMIdedup 1MM_CR \\
         --outFileNamePrefix \"${SAMPLE_ID}_10x_\""
        INPUT_BAM="${ALIGNMENT_OUT_DIR}/${SAMPLE_ID}_10x_Aligned.sortedByCoord.out.bam"
        BARCODE_FILE_UNCOMPRESSED="${ALIGNMENT_OUT_DIR}/${SAMPLE_ID}_10x_Solo.out/Gene/filtered/barcodes.tsv"

    elif [[ "$PROTOCOL" =~ [Dd]rop.?[Ss]eq ]]; then
        # STARsolo parameters for Drop-Seq chemistry.
        local CB_START=1; local CB_LEN=12; local UMI_START=13; local UMI_LEN=8
        STAR_COMMAND+=" \\
         --soloType CB_UMI_Simple \\
         --soloCBstart $CB_START --soloCBlen $CB_LEN \\
         --soloUMIstart $UMI_START --soloUMIlen $UMI_LEN \\
         --soloCBwhitelist None \\
         --soloFeatures Gene \\
         --outFileNamePrefix \"${SAMPLE_ID}_DropSeq_\""
        INPUT_BAM="${ALIGNMENT_OUT_DIR}/${SAMPLE_ID}_DropSeq_Aligned.sortedByCoord.out.bam"
        BARCODE_FILE_UNCOMPRESSED="${ALIGNMENT_OUT_DIR}/${SAMPLE_ID}_DropSeq_Solo.out/Gene/filtered/barcodes.tsv"

    elif [[ "$PROTOCOL" =~ [Ss][Mm][Aa][Rr][Tt]|[Pp]late ]]; then
        # STARsolo parameters for SMART-Seq2 or generic plate-based sequencing.
        # These protocols typically do not have molecular tags or barcodes,
        # so Solo features are omitted, and a barcode file is generally not generated.
        STAR_COMMAND+=" \\
         --outFileNamePrefix \"${SAMPLE_ID}_SMARTSeq2_\""
        INPUT_BAM="${ALIGNMENT_OUT_DIR}/${SAMPLE_ID}_SMARTSeq2_Aligned.sortedByCoord.out.bam"
        BARCODE_FILE_UNCOMPRESSED="" # No barcode file typically needed for cellsnp-lite with SMART-Seq2
    fi

    echo "[$(date)] Executing STARsolo command:"
    echo "$STAR_COMMAND"
    eval "$STAR_COMMAND" # Execute the constructed STARsolo command
    echo "[$(date)] STARsolo finished for ${SAMPLE_ID}."

    # Unload the STAR module
    module unload star/"$STAR_VERSION"

else
    # This branch should ideally not be reached due to initial protocol check.
    echo "ERROR: Alignment protocol logic error. Check script."
    exit 1
fi

# After alignment, return to the base project directory for consistency
cd "$BASE_PROJECT_DIR"

# --- Verify the existence of the aligned BAM file ---
if [ ! -s "${INPUT_BAM}" ]; then
    echo "FATAL ERROR: Aligned BAM file not found or is empty: ${INPUT_BAM}"
    echo "Please check if the alignment step (Cell Ranger or STARsolo) ran successfully."
    exit 1
fi
echo "[$(date)] Aligned BAM file found: ${INPUT_BAM}"

# --- Verify the existence of the barcode file (if applicable) ---
if [[ "$PROTOCOL" =~ [Tt]en[Xx]|[1]0[xX]|[Dd]rop.?[Ss]eq ]]; then
    if [ ! -s "${BARCODE_FILE_UNCOMPRESSED}" ]; then
        echo "FATAL ERROR: Filtered barcodes file (uncompressed) not found or is empty: ${BARCODE_FILE_UNCOMPRESSED}"
        echo "This file is required for cellsnp-lite with 10x Genomics or Drop-Seq protocols."
        exit 1
    fi
    echo "[$(date)] Barcode file found: ${BARCODE_FILE_UNCOMPRESSED}"
fi

---
## Step 3: Variant Calling with cellsnp-lite

echo "[$(date)] Starting variant calling with cellsnp-lite for ${SAMPLE_ID}..."

# Ensure a clean environment before loading the conda environment
module purge
# Activate the conda environment containing cellsnp-lite
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "$CELLSNP_CONDA_ENV"

# --- Debugging: Verify cellsnp-lite executable path ---
echo "[$(date)] Verifying cellsnp-lite environment and executable..."
which cellsnp-lite || { echo "ERROR: cellsnp-lite not found in PATH. Check conda environment '$CELLSNP_CONDA_ENV'."; exit 1; }
echo "--------------------------------------------------------"

# --- Verify reference files for cellsnp-lite ---
if [ ! -s "${SNP_VCF}" ]; then
    echo "FATAL ERROR: SNP VCF file not found or is empty: ${SNP_VCF}"
    exit 1
fi
if [ ! -s "${REFERENCE_FASTA}" ]; then
    echo "FATAL ERROR: Reference FASTA file not found or is empty: ${REFERENCE_FASTA}"
    exit 1
fi
if [ ! -s "${REFERENCE_FASTA}.fai" ]; then
    echo "FATAL ERROR: Reference FASTA index (.fai) not found: ${REFERENCE_FASTA}.fai"
    echo "Please index your FASTA with 'samtools faidx ${REFERENCE_FASTA}'."
    exit 1
fi
echo "[$(date)] All cellsnp-lite reference files verified."

# --- Ensure BAM is indexed ---
# cellsnp-lite requires the input BAM file to be indexed (.bai).
if [ ! -f "${INPUT_BAM}.bai" ]; then
    echo "[$(date)] Indexing input BAM file: ${INPUT_BAM}..."
    # Load Samtools module to perform indexing
    module load samtools/"$SAMTOOLS_VERSION"
    samtools index "$INPUT_BAM"
    # Unload Samtools after use
    module unload samtools/"$SAMTOOLS_VERSION"
else
    echo "[$(date)] BAM file is already indexed."
fi

# --- Gzip the barcodes.tsv if not already gzipped (for 10x/Drop-Seq protocols) ---
# cellsnp-lite expects gzipped barcode files when provided.
BARCODE_FILE_COMPRESSED=""
if [[ "$PROTOCOL" =~ [Tt]en[Xx]|[1]0[xX]|[Dd]rop.?[Ss]eq ]]; then
    BARCODE_FILE_COMPRESSED="${BARCODE_FILE_UNCOMPRESSED}.gz"
    if [ ! -s "${BARCODE_FILE_COMPRESSED}" ]; then
        echo "[$(date)] Gzipping barcodes.tsv: ${BARCODE_FILE_UNCOMPRESSED}..."
        gzip -c "${BARCODE_FILE_UNCOMPRESSED}" > "${BARCODE_FILE_COMPRESSED}"
    else
        echo "[$(date)] Barcodes.tsv is already gzipped: ${BARCODE_FILE_COMPRESSED}"
    fi
fi

# --- Construct and execute the cellsnp-lite command ---
CELLSNP_CMD="cellsnp-lite \\
    -s \"$INPUT_BAM\" \\
    -O \"${SNV_RESULT_DIR}/cellsnp_${SAMPLE_ID}\" \\
    -R \"$SNP_VCF\""

# Add barcode file and cell/UMI tags only if applicable for the protocol
if [[ "$PROTOCOL" =~ [Tt]en[Xx]|[1]0[xX]|[Dd]rop.?[Ss]eq ]]; then
    CELLSNP_CMD+=" \\
    -b \"$BARCODE_FILE_COMPRESSED\" \\
    --cellTAG CB --UMItag UB"
fi

# Add common filtering and output parameters
CELLSNP_CMD+=" \\
    --minMAPQ 30 \\
    --minCOUNT 1 \\
    --minMAF 0.01 \\
    --genotype \\
    --refseq \"$REFERENCE_FASTA\" \\
    --gzip \\
    -p \"$SLURM_CPUS_PER_TASK\"" # Use allocated CPU cores

echo "[$(date)] Running cellsnp-lite command:"
echo "$CELLSNP_CMD"
echo "--------------------------------------------------------"

eval "$CELLSNP_CMD" # Execute the constructed cellsnp-lite command

echo "[$(date)] cellsnp-lite pipeline completed for ${SAMPLE_ID}."
echo "[$(date)] Variant calling results are located in: ${SNV_RESULT_DIR}/cellsnp_${SAMPLE_ID}"

# Deactivate the conda environment to clean up
echo "[$(date)] Deactivating conda environment: ${CELLSNP_CONDA_ENV}..."
conda deactivate

echo "[$(date)] All pipeline steps finished successfully for ${SAMPLE_ID}."