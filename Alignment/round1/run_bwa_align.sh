#!/usr/bin/env bash
set -euo pipefail

THREADS=8
BWA=/home/sgao30/bwa/bwa
BWA_REF=/ShangGaoAIProjects/tools/reference/GRCh38/BWA/genome.fa
TRIMMED=/ShangGaoAIProjects/ZhangJW/Adapter_trimming/round1/step2_trim_fastq
BASE=/ShangGaoAIProjects/ZhangJW/Alignment/round1

SAMPLES=(3PB 3PY 6PY 9PY)

echo "=== Step 2: BWA MEM paired-end alignment ==="
for SAMPLE in "${SAMPLES[@]}"; do
    KEY="260430R-S-XY-${SAMPLE}"
    echo "  $SAMPLE"
    $BWA mem \
        -t $THREADS \
        $BWA_REF \
        $TRIMMED/${KEY}_R1.fq.gz \
        $TRIMMED/${KEY}_R2.fq.gz \
        2> $BASE/step2_align/${KEY}.bwa.log \
    | samtools sort -@ $THREADS -O bam -o $BASE/step2_align/${KEY}.sorted.bam
    samtools index $BASE/step2_align/${KEY}.sorted.bam
    echo "    done -> ${KEY}.sorted.bam"
done

echo "=== All samples aligned ==="
