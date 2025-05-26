import pysam
import os
import collections
import sys

def split_bam_by_cell_barcode(
    input_bam_path: str,
    output_dir: str = "per_cell_bams_samtools/",
    cell_barcode_tag: str = "XC",
    min_num_reads_per_tag: int = 100
):
    """
    Splits a BAM file into multiple per-cell BAM files based on a specified
    cell barcode tag, filtering out cells below a minimum read threshold.

    Args:
        input_bam_path (str): Path to the input BAM file.
        output_dir (str): Directory where the per-cell BAM files will be saved.
                          Defaults to 'per_cell_bams_samtools/'.
        cell_barcode_tag (str): The BAM tag used for cell barcodes (e.g., 'XC', 'CB').
                                Defaults to 'XC'.
        min_num_reads_per_tag (int): Minimum number of reads a cell barcode must have
                                     to be included. Defaults to 100.
    """
    print(f"Starting BAM splitting process for: {input_bam_path}")
    print(f"Output directory: {output_dir}")
    print(f"Cell barcode tag: {cell_barcode_tag}")
    print(f"Minimum reads per cell: {min_num_reads_per_tag}")

    ## Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    print(f"Ensured output directory '{output_dir}' exists.")

    ## Count reads per cell barcode 
    print("\n--- Pass 1: Counting reads per cell barcode ---")
    cell_barcode_counts = collections.Counter()
    total_reads_processed = 0

    try:
        with pysam.AlignmentFile(input_bam_path, "rb") as bamfile:
            for read in bamfile:
                total_reads_processed += 1
                if read.has_tag(cell_barcode_tag):
                    barcode = read.get_tag(cell_barcode_tag)
                    cell_barcode_counts[barcode] += 1
                if total_reads_processed % 1_000_000 == 0:
                    sys.stdout.write(f"\rProcessed {total_reads_processed:,} reads...")
                    sys.stdout.flush()
        sys.stdout.write(f"\rFinished Pass 1. Total reads processed: {total_reads_processed:,}\n")
        sys.stdout.flush()

    except FileNotFoundError:
        print(f"Error: Input BAM file not found at '{input_bam_path}'")
        return
    except Exception as e:
        print(f"Error during Pass 1 (counting barcodes): {e}")
        return

    ## Filter cell barcodes based on minimum read count
    valid_cell_barcodes = {
        barcode for barcode, count in cell_barcode_counts.items()
        if count >= min_num_reads_per_tag
    }

    print(f"\nFound {len(cell_barcode_counts):,} unique cell barcodes.")
    print(f"Identified {len(valid_cell_barcodes):,} valid cell barcodes (>= {min_num_reads_per_tag} reads).")

    if not valid_cell_barcodes:
        print("No valid cell barcodes found meeting the minimum read threshold. Exiting.")
        return

    ## Write reads to per-cell BAM files ---
    print("\n--- Pass 2: Writing reads to per-cell BAM files ---")
    output_bam_handles = {}
    reads_written_count = 0

    try:
        with pysam.AlignmentFile(input_bam_path, "rb") as bamfile:
            ## Get the header from the input BAM file
            header = bamfile.header

            for read in bamfile:
                if read.has_tag(cell_barcode_tag):
                    barcode = read.get_tag(cell_barcode_tag)

                    if barcode in valid_cell_barcodes:
                        if barcode not in output_bam_handles:
                            ## Open a new BAM file for this barcode if not already open
                            output_filename = os.path.join(output_dir, f"{barcode}.bam")
                            print(f"Opening new output BAM: {output_filename}")
                            output_bam_handles[barcode] = pysam.AlignmentFile(
                                output_filename, "wb", header=header
                            )
                        ## Write the read to the appropriate output BAM file
                        output_bam_handles[barcode].write(read)
                        reads_written_count += 1

                if reads_written_count % 1_000_000 == 0 and reads_written_count > 0:
                    sys.stdout.write(f"\rWritten {reads_written_count:,} reads...")
                    sys.stdout.flush()

        sys.stdout.write(f"\rFinished Pass 2. Total reads written: {reads_written_count:,}\n")
        sys.stdout.flush()

    except Exception as e:
        print(f"Error during Pass 2 (writing reads): {e}")
    finally:
        ## Close all opened output BAM files
        for barcode, handle in output_bam_handles.items():
            handle.close()
        print("\nAll output BAM files closed.")

    print("\nBAM splitting process completed.")
    print("Note: This script does not create separate 'temp files' that need cleanup.")
    print("The output BAM files are directly written to the specified output directory.")


## input with sys.argv 
if __name__ == "__main__":
    ## Expected arguments:
    ## sys.argv[1]: input_bam_path (required)
    ## sys.argv[2]: output_dir (optional, default: "per_cell_bams_samtools/")
    ## sys.argv[3]: cell_barcode_tag (optional, default: "XC")
    ## sys.argv[4]: min_num_reads_per_tag (optional, default: 100)

    if len(sys.argv) < 2:
        print("Usage: python your_script_name.py <input_bam_path> [output_dir] [barcode_tag] [min_reads]")
        print("Example: python your_script_name.py my_data.bam per_cell_output CB 50")
        sys.exit(1)

    input_bam = sys.argv[1]
    output_directory = sys.argv[2] if len(sys.argv) > 2 else "per_cell_bams_samtools/"
    barcode_tag = sys.argv[3] if len(sys.argv) > 3 else "XC"
    min_reads_per_tag = int(sys.argv[4]) if len(sys.argv) > 4 else 100

    split_bam_by_cell_barcode(
        input_bam_path=input_bam,
        output_dir=output_directory,
        cell_barcode_tag=barcode_tag,
        min_num_reads_per_tag=min_reads_per_tag
    )