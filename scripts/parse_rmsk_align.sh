#!/bin/bash

# Configuration
INPUT_FILE="hs1.fa.align.gz"
OUTPUT_FILE="hs1_simplified_alignments.txt"

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: $INPUT_FILE not found."
    exit 1
fi

echo "Processing $INPUT_FILE..."
echo "Extracting headers and summary statistics to $OUTPUT_FILE..."

# We use zcat to stream the compressed file.
# Refactored awk logic:
# 1. We no longer print /^$/ (empty lines) from the source file.
# 2. We track when we are printing a header to insert a single spacer for readability.

zcat "$INPUT_FILE" | awk '
    # Match the header line (starts with the SW score)
    /^[0-9]+/ && NF >= 10 {
        # Add a single empty line spacer before every entry except the very first one
        if (count > 0) {
            print ""
        }
        print $0
        count++
        next
    }
    
    # Match specific summary metadata lines
    /^Matrix =/ || \
    /^Kimura/ || \
    /^CpG sites/ || \
    /^Transitions \/ transversions/ || \
    /^Gap_init rate/ {
        print $0
    }
' > "$OUTPUT_FILE"

echo "Done! Simplified file saved as $OUTPUT_FILE"
# Calculate size reduction
if [ -f "$OUTPUT_FILE" ]; then
    ORIG_SIZE=$(du -h "$INPUT_FILE" | cut -f1)
    NEW_SIZE=$(du -h "$OUTPUT_FILE" | cut -f1)
    echo "Original size (compressed): $ORIG_SIZE"
    echo "New size (uncompressed): $NEW_SIZE"
fi
