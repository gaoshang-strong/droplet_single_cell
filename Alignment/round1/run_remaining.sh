#!/usr/bin/env bash
set -euo pipefail

# Run steps 2-5 for samples 2PB through 9PY (1PB already done)

THREADS=8
BOWTIE2_INDEX=/ShangGaoAIProjects/tools/reference/GRCh38/Bowtie2/genome
TRIMMED=/ShangGaoAIProjects/ZhangJW/Adapter_trimming/round1/step2_trim_fastq
BASE=/ShangGaoAIProjects/ZhangJW/Alignment/round1
PYTHON=/home/sgao30/micromamba/envs/droplet_single_cell/bin/python
PICARD="java -jar /home/sgao30/picard.jar"

SAMPLES=(2PB 3PB 3PY 6PY 9PY)

for SAMPLE in "${SAMPLES[@]}"; do
    KEY="260430R-S-XY-${SAMPLE}"
    echo "=== Processing $SAMPLE ==="

    echo "  Step 2: align"
    bowtie2 \
        --very-sensitive -X 2000 -p $THREADS \
        -x $BOWTIE2_INDEX \
        -1 $TRIMMED/${KEY}_R1.fq.gz \
        -2 $TRIMMED/${KEY}_R2.fq.gz \
        2> $BASE/step2_align/${KEY}.bowtie2.log \
    | samtools sort -@ $THREADS -O bam -o $BASE/step2_align/${KEY}.sorted.bam
    samtools index $BASE/step2_align/${KEY}.sorted.bam

    echo "  Step 3: add CB tag"
    $PYTHON $BASE/add_cb_tag.py \
        --bam-in  $BASE/step2_align/${KEY}.sorted.bam \
        --bam-out $BASE/step3_tag_bam/${KEY}.tagged.bam \
        --cb-map  $BASE/step1_add_cb/${KEY}_cb_map.tsv.gz
    samtools index $BASE/step3_tag_bam/${KEY}.tagged.bam

    echo "  Step 4: filter"
    CHROMS=$(samtools view -H $BASE/step3_tag_bam/${KEY}.tagged.bam \
        | grep "^@SQ" | cut -f2 | sed 's/SN://' | grep -v "^chrM$" | tr '\n' ' ')
    samtools view -@ $THREADS -b -f 0x2 -F 0xC -q 30 \
        $BASE/step3_tag_bam/${KEY}.tagged.bam $CHROMS \
        -o $BASE/step4_filter/${KEY}.filtered.bam
    samtools index $BASE/step4_filter/${KEY}.filtered.bam

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

echo "=== All remaining samples complete ==="
