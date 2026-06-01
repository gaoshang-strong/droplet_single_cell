#!/usr/bin/env bash
set -euo pipefail

THREADS=8
BOWTIE2=/home/sgao30/bowtie2-2.4.4-source/bowtie2-2.4.4/bowtie2
BT2_REF=/ShangGaoAIProjects/tools/reference/GRCm39/Bowtie2/genome
TRIMMED=/ShangGaoAIProjects/ZhangJW/Adapter_trimming/round1/step2_trim_fastq
OUT=/ShangGaoAIProjects/ZhangJW/Alignment/round1/step2b_mouse_genome

SAMPLES=(3PY 6PY 9PY)

echo "=== Bowtie2 alignment to GRCm39 (PY samples) ==="
for SAMPLE in "${SAMPLES[@]}"; do
    KEY="260430R-S-XY-${SAMPLE}"
    echo "  $SAMPLE"
    $BOWTIE2 -p $THREADS -x ${BT2_REF} \
        -1 ${TRIMMED}/${KEY}_R1.fq.gz \
        -2 ${TRIMMED}/${KEY}_R2.fq.gz \
        2> ${OUT}/${KEY}.bowtie2.log \
    | samtools sort -@ $THREADS -O bam -o ${OUT}/${KEY}.bowtie2.sorted.bam
    samtools index ${OUT}/${KEY}.bowtie2.sorted.bam
    echo "    done -> ${KEY}.bowtie2.sorted.bam"
done

echo "=== All done ==="
