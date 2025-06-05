# Single-Cell Variant Calling Pipeline

This repository contains a modular pipeline for variant calling on single-cell RNA-seq data, using Cell Ranger, a Python-based BAM splitting tool, and [cellsnp-lite](https://github.com/single-cell-genetics/cellsnp-lite).

## Features

- Compatible with 10x Genomics (Cell Ranger) or Drop-seq BAMs.
- Automatically splits pooled BAM into per-cell BAMs.
- Efficient SNP/variant calling per cell using cellsnp-lite.
- User-friendly: prompts for inputs, no hardcoded personal paths.

## Requirements

- [Cell Ranger](https://support.10xgenomics.com/)
- Python 3.6+ with `pysam`
- [samtools](http://www.htslib.org/)
- [cellsnp-lite](https://github.com/single-cell-genetics/cellsnp-lite)
- Unix/Linux shell

## Reference File Preparation

You will need:

- **Reference genome** directory (e.g., from 10x Genomics or built via `cellranger mkref`)
- **SNP VCF** file (e.g., dbSNP or 1000 Genomes, bgzipped and tabix-indexed)

**Example commands (for human):**

```bash
# Download Cell Ranger reference (GRCh38)
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -xzf refdata-gex-GRCh38-2020-A.tar.gz

# Download and index SNP VCF
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.39.gz
tabix -p vcf GCF_000001405.39.gz
```

- For other species, obtain reference and SNP VCFs from Ensembl, NCBI, or appropriate consortia.
- Set the reference genome path in Cell Ranger (`--transcriptome=...`).
- Provide the SNP VCF path to the pipeline when prompted.


## Processing Barcoded Single-Cell Data (10x, Drop-seq)

For data containing single-cell barcodes (e.g., 10x Genomics, Drop-seq), use our integrated SLURM pipeline:

### Introduction: scRNA_VariantCalling_PipelineV2.0.sh

`scRNA_VariantCalling_PipelineV2.0.sh` is an end-to-end pipeline that takes raw FASTQ files through alignment (using Cell Ranger or STARsolo), performs quality control, and runs per-cell variant calling using cellsnp-lite.

It's designed for use on HPC clusters with SLURM and supports 10x Genomics, Drop-seq, and plate-based protocols.

### Key Features

* **detects protocol**: Supports 10x Genomics, Drop-Seq, and SMART-Seq2/plate data.
* **End-to-end**: Handles FastQC, alignment, barcode handling, and per-cell variant calling.
* **No hardcoded paths or user info**: Adaptable for any lab or HPC setup.
* Compatible with modular environments and conda.

### Usage

Submit your job (edit the script to match your paths and reference files):

```bash
sbatch scRNA_VariantCalling_PipelineV2.0.sh <SAMPLE_ID> <PROTOCOL>
```

## Processing sc-RNA without barcode
### 1. Process FASTQ with Cell Ranger or any alignment tools

Generate a coordinate-sorted BAM file with Tag:

```bash
cellranger count \
  --id=YOUR_SAMPLE_ID \
  --transcriptome=/path/to/reference \
  --fastqs=/path/to/fastqs \
  --sample=YOUR_SAMPLE_NAME
```

---

### 2. Run the Split BAM by Cell Barcode & Variant Calling Pipeline

Run the bash pipeline:

```bash
bash bam2snp.sh
```

The first part of the pipeline contains the script that split bam by cell barcode

```bash
python splitbam.py
```
The script will prompt for:
- Sample name (for output labeling)
- Path to input BAM (e.g., Cell Ranger’s `possorted_genome_bam.bam`)
- Barcode tag (`CB` for 10x, `XC` for Drop-seq)
- Minimum reads per barcode to keep

---
The second half picks up the output of bamsplit.py and perform variant calling with cellsnp-lite

The script will prompt for:
- Output directories for per-cell BAMs and SNP results
- Reference SNP VCF (e.g., dbSNP)

---

### 4. Results

Variant calls per cell will be in the specified output directory and ready for downstream analysis.

---

## Notes

- Indexing of per-cell BAMs is automatic.
- Intermediate files (e.g., list files) are cleaned up at the end of the run.
- To save storage, you can optionally delete per-cell BAMs after variant calling.

---

## Citing Tools

Please cite:
- [Cell Ranger](https://support.10xgenomics.com/)
- [samtools](http://www.htslib.org/)
- [cellsnp-lite](https://github.com/single-cell-genetics/cellsnp-lite)
- [pysam](https://github.com/pysam-developers/pysam)
- [STAR](https://github.com/alexdobin/STAR)

---

## License

See [LICENSE](LICENSE) for terms.

---

## Contact

For questions, open a GitHub issue or contact the repository maintainer.
