#!/usr/bin/env bash
#
# normalize_bam_to_bw.sh
#
# A script to normalize a BAM file using RPKM plus a custom scale factor,
# and produce a BigWig coverage track.
#
# Usage:
#   ./normalize_bam_to_bw.sh [OPTIONS]
#
# Required arguments:
#   -i, --input-bam        Path to the input BAM file
#   -o, --output-bw        Path to the output BigWig file
#   -s, --scale-factor     Scale factor to apply in addition to RPKM
#
# Optional arguments:
#   -b, --bin-size         Bin size for coverage (default: 25)
#   -d, --ignore-duplicates  Ignore PCR duplicates (on by default, pass '--no-ignore-duplicates' to disable)
#
# Example:
#   ./normalize_bam_to_bw.sh \\
#       --input-bam sample.bam \\
#       --output-bw sample.bw \\
#       --scale-factor 0.5 \\
#       --bin-size 25
#

# Default values
BIN_SIZE=25
IGNORE_DUPLICATES=true

# Parse arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -i|--input-bam)
            INPUT_BAM="$2"
            shift
            shift
            ;;
        -o|--output-bw)
            OUTPUT_BW="$2"
            shift
            shift
            ;;
        -s|--scale-factor)
            SCALE_FACTOR="$2"
            shift
            shift
            ;;
        -b|--bin-size)
            BIN_SIZE="$2"
            shift
            shift
            ;;
        -d|--ignore-duplicates)
            IGNORE_DUPLICATES=true
            shift
            ;;
        --no-ignore-duplicates)
            IGNORE_DUPLICATES=false
            shift
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Check required arguments
if [ -z "$INPUT_BAM" ] || [ -z "$OUTPUT_BW" ] || [ -z "$SCALE_FACTOR" ]; then
    echo "Missing required arguments!"
    echo "Usage example: ./normalize_bam_to_bw.sh -i sample.bam -o sample.bw -s 1"
    exit 1
fi

# Build bamCoverage command
COMMAND="bamCoverage -b ${INPUT_BAM} -o ${OUTPUT_BW} --normalizeUsing RPKM --binSize ${BIN_SIZE} --scaleFactor ${SCALE_FACTOR}"

if [ "$IGNORE_DUPLICATES" = true ]; then
    COMMAND="${COMMAND} --ignoreDuplicates"
fi

# Print and run
echo "Running command:"
echo "${COMMAND}"
eval "${COMMAND}"
