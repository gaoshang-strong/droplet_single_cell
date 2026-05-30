#!/bin/bash

set -euo pipefail

DATA_DIR="/ShangGaoAIProjects/ZhangJW/data"
OUTPUT_DIR="/ShangGaoAIProjects/ZhangJW/fastQC"
THREADS=12

mkdir -p "$OUTPUT_DIR"

mapfile -t FILES < <(find "$DATA_DIR" -name "*.gz")

echo "Found ${#FILES[@]} files to process:"
printf '  %s\n' "${FILES[@]}"
echo ""

fastqc "${FILES[@]}" \
    --outdir "$OUTPUT_DIR" \
    --threads "$THREADS"

echo ""
echo "Done. Results saved to $OUTPUT_DIR"
