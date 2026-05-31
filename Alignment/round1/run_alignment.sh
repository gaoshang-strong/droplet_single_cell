#!/usr/bin/env bash
set -euo pipefail

THREADS=8
BOWTIE2_INDEX=/ShangGaoAIProjects/tools/reference/GRCh38/Bowtie2/genome
TRIMMED=/ShangGaoAIProjects/ZhangJW/Adapter_trimming/round1/step2_trim_fastq
HD2=/ShangGaoAIProjects/ZhangJW/R1_parser/round3/filter_HD2
BASE=/ShangGaoAIProjects/ZhangJW/Alignment/round1
SCRIPTS=$BASE
PYTHON=/home/sgao30/micromamba/envs/droplet_single_cell/bin/python
PICARD="java -jar /home/sgao30/picard.jar"

SAMPLES=(1PB 2PB 3PB 3PY 6PY 9PY)

# Step 0: extract CB maps from filter_HD2 R1 headers (fast, read-only scan)
echo "=== Step 0: Extract CB maps ==="
for SAMPLE in "${SAMPLES[@]}"; do
    KEY="260430R-S-XY-${SAMPLE}"
    echo "  $SAMPLE"
    $PYTHON $SCRIPTS/extract_cb_map.py \
        --r1 $HD2/${KEY}_R1.fq.gz \
        --out $BASE/step1_add_cb/${KEY}_cb_map.tsv.gz
done
echo "  Done"

for SAMPLE in "${SAMPLES[@]}"; do
    KEY="260430R-S-XY-${SAMPLE}"
    echo "=== Processing $SAMPLE ==="

    # Step 2: align trimmed FASTQs directly with Bowtie2
    echo "  Step 2: align"
    bowtie2 \
        --very-sensitive -X 2000 -p $THREADS \
        -x $BOWTIE2_INDEX \
        -1 $TRIMMED/${KEY}_R1.fq.gz \
        -2 $TRIMMED/${KEY}_R2.fq.gz \
        2> $BASE/step2_align/${KEY}.bowtie2.log \
    | samtools sort -@ $THREADS -O bam -o $BASE/step2_align/${KEY}.sorted.bam
    samtools index $BASE/step2_align/${KEY}.sorted.bam

    # Step 3: add CB/UMI tags from TSV map
    echo "  Step 3: add CB tag to BAM"
    $PYTHON $SCRIPTS/add_cb_tag.py \
        --bam-in  $BASE/step2_align/${KEY}.sorted.bam \
        --bam-out $BASE/step3_tag_bam/${KEY}.tagged.bam \
        --cb-map  $BASE/step1_add_cb/${KEY}_cb_map.tsv.gz
    samtools index $BASE/step3_tag_bam/${KEY}.tagged.bam

    # Step 4: filter (MAPQ>=30, proper pairs, no chrM)
    echo "  Step 4: filter BAM"
    CHROMS=$(samtools view -H $BASE/step3_tag_bam/${KEY}.tagged.bam \
        | grep "^@SQ" | cut -f2 | sed 's/SN://' | grep -v "^chrM$" | tr '\n' ' ')
    samtools view -@ $THREADS -b -f 0x2 -F 0xC -q 30 \
        $BASE/step3_tag_bam/${KEY}.tagged.bam $CHROMS \
        -o $BASE/step4_filter/${KEY}.filtered.bam
    samtools index $BASE/step4_filter/${KEY}.filtered.bam

    # Step 5: deduplicate per cell barcode
    echo "  Step 5: deduplicate"
    $PICARD MarkDuplicates \
        -I  $BASE/step4_filter/${KEY}.filtered.bam \
        -O  $BASE/step5_dedup/${KEY}.dedup.bam \
        -M  $BASE/step5_dedup/${KEY}.dup_metrics.txt \
        --BARCODE_TAG CB \
        --REMOVE_DUPLICATES true
    samtools index $BASE/step5_dedup/${KEY}.dedup.bam

    echo "  Done: $SAMPLE"
done

echo "=== All samples complete ==="
