#!/usr/bin/env bash
set -euo pipefail

THREADS=8
BWA=/home/sgao30/bwa/bwa
BOWTIE2=/home/sgao30/bowtie2-2.4.4-source/bowtie2-2.4.4/bowtie2
BWA_REF=/ShangGaoAIProjects/tools/reference/GRCm39/BWA/genome
BT2_REF=/ShangGaoAIProjects/tools/reference/GRCm39/Bowtie2/genome
TRIMMED=/ShangGaoAIProjects/ZhangJW/Adapter_trimming/round1/step2_trim_fastq
OUT=/ShangGaoAIProjects/ZhangJW/Alignment/round1/step2b_mouse_genome

SAMPLES=(1PB 2PB 3PB)

echo "=== BWA MEM alignment to GRCm39 ==="
for SAMPLE in "${SAMPLES[@]}"; do
    KEY="260430R-S-XY-${SAMPLE}"
    echo "  $SAMPLE"
    $BWA mem -t $THREADS \
        ${BWA_REF} \
        ${TRIMMED}/${KEY}_R1.fq.gz \
        ${TRIMMED}/${KEY}_R2.fq.gz \
        2> ${OUT}/${KEY}.bwa.log \
    | samtools sort -@ $THREADS -O bam -o ${OUT}/${KEY}.bwa.sorted.bam
    samtools index ${OUT}/${KEY}.bwa.sorted.bam
    echo "    done -> ${KEY}.bwa.sorted.bam"
done

echo "=== Bowtie2 alignment to GRCm39 ==="
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
