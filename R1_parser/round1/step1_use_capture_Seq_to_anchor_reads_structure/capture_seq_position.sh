#!/bin/bash

set -euo pipefail

DATA_DIR="/ShangGaoAIProjects/ZhangJW/data"
OUTPUT_DIR="/ShangGaoAIProjects/ZhangJW/R1_parser/round1/step1_use_capture_Seq_to_anchor_reads_structure"

mkdir -p "$OUTPUT_DIR"

mapfile -t R1_FILES < <(find "$DATA_DIR" -name "*_R1.fq.gz")

echo "Found ${#R1_FILES[@]} R1 files:"
printf '  %s\n' "${R1_FILES[@]}"
echo ""

for R1 in "${R1_FILES[@]}"; do
    SAMPLE=$(basename "$R1" _R1.fq.gz)
    echo "Processing: $SAMPLE"

    SUMMARY="$OUTPUT_DIR/${SAMPLE}_anchor_summary.tsv"
    POS_DIST="$OUTPUT_DIR/${SAMPLE}_anchor_position_distribution.tsv"
    POS_SORTED="$OUTPUT_DIR/${SAMPLE}_anchor_position_distribution.sorted.tsv"

    pigz -dc -p 8 "$R1" \
    | awk -v summary="$SUMMARY" -v pos_dist="$POS_DIST" '
    BEGIN {
        anchor = "CACCGTCTCCGCCTC"
        total = 0
        hit = 0
    }
    NR % 4 == 2 {
        total++
        pos = index($0, anchor)
        if (pos > 0) {
            hit++
            pos_count[pos]++
        }
    }
    END {
        printf "anchor\t%s\n", anchor > summary
        printf "total_read1\t%d\n", total >> summary
        printf "anchor_exact_hit\t%d\n", hit >> summary
        printf "anchor_exact_hit_rate\t%.6f\n", hit / total >> summary

        print "pos\tcount\trate_among_hits" > pos_dist
        for (p in pos_count) {
            printf "%d\t%d\t%.6f\n", p, pos_count[p], pos_count[p] / hit >> pos_dist
        }
    }
    '

    sort -k1,1n "$POS_DIST" > "$POS_SORTED"
    echo "  -> $SUMMARY"
    echo "  -> $POS_SORTED"
done

echo ""
echo "Done. All results saved to $OUTPUT_DIR"
