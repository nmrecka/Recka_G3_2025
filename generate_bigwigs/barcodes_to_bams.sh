#!/bin/bash

# Directory containing the _barcodes.txt files.
BARCODES_DIR="/Users/nmrecka/projects/CTRL_only/Barcodes"

# The input BAM file (assumes it contains a CB tag for cell barcodes).
INPUT_BAM="/Users/nmrecka/projects/2023_10X_multiomics_E15.5_murine_skin/raw_data/dnacore454.healthcare.uiowa.edu/CTRL_count/outs/gex_possorted_bam.bam"

# Loop over every barcode file in the directory.
for barcode_file in "$BARCODES_DIR"/*_barcodes.txt; do
    # Extract the base name (e.g., "epidermis" from "epidermis_barcodes.txt")
    base=$(basename "$barcode_file" _barcodes.txt)
    output_bam="${base}.bam"
    
    echo "Processing ${barcode_file} -> ${output_bam}"
    
    # Filter the BAM file based on the barcodes in the file.
    samtools view -h "$INPUT_BAM" | \
    awk -v barcode_file="$barcode_file" 'BEGIN {
        # Read in all the barcodes from the file into an array.
        while ((getline barcode < barcode_file) > 0) { b[barcode] = 1 }
    }
    # Always print header lines.
    /^@/ { print; next }
    # For each alignment line, extract the CB tag.
    {
        if (match($0, /CB:Z:([^ \t]+)/, a)) {
            if (a[1] in b) {
                print;
            }
        }
    }' | \
    samtools view -b -o "$output_bam"
done