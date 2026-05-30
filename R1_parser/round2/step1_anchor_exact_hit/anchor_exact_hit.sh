#!/bin/bash
# Round2 Step1: exact hit scan for the combined 23 bp anchor
# Anchor = common_fixed (TAAGGCGA) + capture (CACCGTCTCCGCCTC)
# Expected position: starts at pos 28 (1-based)

set -euo pipefail

DATA_DIR="/ShangGaoAIProjects/ZhangJW/data"
OUTPUT_DIR="/ShangGaoAIProjects/ZhangJW/R1_parser/round2/step1_anchor_exact_hit"
ANCHOR="TAAGGCGACACCGTCTCCGCCTC"

mkdir -p "$OUTPUT_DIR"

mapfile -t R1_FILES < <(find "$DATA_DIR" -name "*_R1.fq.gz")

echo "Anchor: $ANCHOR (${#ANCHOR} bp)"
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
    | awk -v anchor="$ANCHOR" -v summary="$SUMMARY" -v pos_dist="$POS_DIST" '
    BEGIN {
        total = 0
        hit   = 0
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
        printf "anchor\t%s\n",              anchor  > summary
        printf "anchor_len\t%d\n",          length(anchor) >> summary
        printf "total_read1\t%d\n",         total   >> summary
        printf "anchor_exact_hit\t%d\n",    hit     >> summary
        printf "anchor_exact_hit_rate\t%.6f\n", hit/total >> summary

        print "pos\tcount\trate_among_hits" > pos_dist
        for (p in pos_count) {
            printf "%d\t%d\t%.6f\n", p, pos_count[p], pos_count[p]/hit >> pos_dist
        }
    }
    '

    sort -k1,1n "$POS_DIST" > "$POS_SORTED"
    echo "  -> $SUMMARY"
    echo "  -> $POS_SORTED"
done

echo ""
echo "Done. Results in $OUTPUT_DIR"
