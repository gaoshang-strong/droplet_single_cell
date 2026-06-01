# BWA Alignment Summary — Round 1

**Reference:** GRCh38 (`/ShangGaoAIProjects/tools/reference/GRCh38/BWA/genome.fa`)  
**Aligner:** BWA MEM (`-t 8`)  
**Date:** 2026-05-31  

---

## Alignment Statistics

| Sample | Read Pairs | Total Reads | Mapped Reads | Mapping Rate | Properly Paired | Properly Paired Rate | Supplementary | BAM Size |
|--------|-----------|------------|-------------|-------------|----------------|---------------------|---------------|----------|
| 1PB | 19,273,973 | 38,547,946 | 6,740,288 | 17.35% | 3,300,748 | 8.56% | 296,124 | 2.1G |
| 2PB | 26,444,544 | 52,889,088 | 11,385,553 | 21.38% | 6,791,242 | 12.84% | 375,749 | 2.8G |
| 3PB | 40,035,523 | 80,071,046 | 18,154,994 | 22.51% | 11,251,678 | 14.05% | 579,102 | 4.3G |
| 3PY | 37,422,718 | 74,845,436 | 44,017,561 | 58.03% | 38,584,444 | 51.55% | 1,002,837 | 3.0G |
| 6PY | 10,485,405 | 20,970,810 | 12,640,220 | 59.85% | 11,264,148 | 53.71% | 149,456 | 815M |
| 9PY | 12,780,806 | 25,561,612 | 11,771,607 | 45.45% | 9,601,468 | 37.56% | 337,597 | 1.2G |

> **Mapped Rate** = mapped reads / total reads in BAM (including supplementary).  
> **Properly Paired Rate** = properly paired / paired in sequencing.

---

## Key Observations

### 1. PB vs PY Mapping Rate Gap

PB samples (1PB, 2PB, 3PB) show dramatically lower mapping rates (17–23%) compared to PY samples (45–60%). This is consistent with the findings from the previous Bowtie2 run and the root cause analysis in `alignment_QC_report.md`: PB samples likely contain a high proportion of non-human reads (e.g., bacterial/viral contamination or non-nuclear DNA) or reads derived from non-reference genomic regions.

| Group | Mapping Rate Range | Properly Paired Rate Range |
|-------|--------------------|---------------------------|
| PB (1PB, 2PB, 3PB) | 17.35% – 22.51% | 8.56% – 14.05% |
| PY (3PY, 6PY, 9PY) | 45.45% – 59.85% | 37.56% – 53.71% |

### 2. Properly Paired Rate

Properly paired rates are notably low for PB samples (<15%). This suggests many mapped reads from PB samples either map to different chromosomes or have unexpected insert sizes, which may reflect chimeric fragments or library quality issues.

### 3. Supplementary Alignments

All samples show supplementary reads (chimeric/split-read alignments). PY samples tend to have proportionally more supplementary reads relative to their total mapped reads, possibly reflecting more structural variation or longer insert fragments captured.

### 4. Read Depth

3PB is the largest library (40M pairs), while 6PY is the smallest (10.5M pairs).

---

## Output Files

All sorted BAMs and indexes are in `step2_align/`:

```
step2_align/
├── 260430R-S-XY-1PB.sorted.bam  (.bai)
├── 260430R-S-XY-2PB.sorted.bam  (.bai)
├── 260430R-S-XY-3PB.sorted.bam  (.bai)
├── 260430R-S-XY-3PY.sorted.bam  (.bai)
├── 260430R-S-XY-6PY.sorted.bam  (.bai)
└── 260430R-S-XY-9PY.sorted.bam  (.bai)
```
