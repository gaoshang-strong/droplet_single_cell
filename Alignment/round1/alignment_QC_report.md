# Alignment QC Report

## Summary Table

| Sample | Input reads | Alignment rate | Aligned reads | After filter | After dedup | Dup rate |
|--------|------------|---------------|---------------|-------------|-------------|----------|
| 1PB | 19,273,973 | 1.55% | 192,266 | 156,323 | 15,260 | 90.2% |
| 2PB | 26,444,544 | 3.86% | 1,021,319 | 449,355 | 34,104 | 92.4% |
| 3PB | 40,035,523 | 4.72% | 1,891,676 | 1,553,424 | 115,566 | 92.6% |
| 3PY | 37,422,718 | 34.27% | 12,825,552 | 7,713,769 | 562,406 | 92.7% |
| 6PY | 10,485,405 | 38.69% | 4,056,875 | 1,673,593 | 215,376 | 87.1% |
| 9PY | 12,780,806 | 24.07% | 3,078,320 | 1,244,050 | 179,530 | 85.6% |

> Note: "Aligned reads" = read pairs × 2 (BAM contains both mates). "After filter" = after MAPQ ≥ 30 + proper pairs + no chrM. "Dup rate" = (filtered − dedup) / filtered.

---

## Bowtie2 Alignment Rate

### PB samples — critically low alignment rate

| Sample | Concordant unique | Concordant multi | Overall |
|--------|------------------|-----------------|---------|
| 1PB | 0.91% | 0.09% | **1.55%** |
| 2PB | 1.81% | 0.08% | **3.86%** |
| 3PB | 3.65% | 0.44% | **4.72%** |

**1PB: 99.44% of individual reads fail to align to hg38.** This is far below the expected >60% for ATAC-seq. PB samples appear to have a systematic issue.

### PY samples — moderate alignment rate

| Sample | Concordant unique | Concordant multi | Overall |
|--------|------------------|-----------------|---------|
| 3PY | 22.49% | 0.65% | **34.27%** |
| 6PY | 21.70% | 3.66% | **38.69%** |
| 9PY | 13.90% | 0.85% | **24.07%** |

PY samples align at 24–39%, still lower than expected but far better than PB samples.

---

## Read Retention Through Pipeline

```
Input reads → Aligned → Filter (MAPQ/chrM/proper pair) → Dedup
```

| Sample | Input | Aligned | Retention to filter | Retention to dedup |
|--------|-------|---------|---------------------|-------------------|
| 1PB | 19.3M | 192K | 0.81% | 0.08% |
| 2PB | 26.4M | 1.0M | 1.70% | 0.13% |
| 3PB | 40.0M | 1.9M | 3.88% | 0.29% |
| 3PY | 37.4M | 12.8M | **20.6%** | **1.50%** |
| 6PY | 10.5M | 4.1M | **15.9%** | **2.05%** |
| 9PY | 12.8M | 3.1M | **9.74%** | **1.40%** |

---

## Duplication Rate

All samples show **>85% duplication rate** after filtering, indicating low library complexity. This is common in scATAC-seq (each cell contributes few unique fragments), but very high duplication may also reflect low input DNA or library preparation issues.

| Sample | Filtered reads | Dedup reads | Dup rate |
|--------|---------------|------------|----------|
| 1PB | 156,323 | 15,260 | 90.2% |
| 2PB | 449,355 | 34,104 | 92.4% |
| 3PB | 1,553,424 | 115,566 | 92.6% |
| 3PY | 7,713,769 | 562,406 | 92.7% |
| 6PY | 1,673,593 | 215,376 | 87.1% |
| 9PY | 1,244,050 | 179,530 | 85.6% |

---

## Key Findings & Concerns

### 1. PB samples have extremely low alignment rate (< 5%)
The vast majority of PB reads (>99% for 1PB) fail to align to hg38, even as individual reads. Possible causes:
- Reads contain residual adapter/non-genomic sequence not fully removed by trimming
- Sample contamination (non-human DNA)
- Library preparation failure specific to PB samples
- R1 orientation issue (reads the wrong strand after trimming)

**Action required**: Blast a subset of unaligned 1PB R1 reads to identify their origin.

### 2. PY samples have moderate but acceptable alignment rate
34–39% alignment rate is below the typical ATAC-seq benchmark (>60%), but may reflect the short R1 insert after barcode trimming and high barcode complexity.

### 3. High duplication rate across all samples
90–93% duplication rate is concerning. Suggests low library complexity or sequencing the same molecules many times. Cell calling (Step 6) will determine if enough unique fragments per cell remain for analysis.

---

## Read Length (Trimmed FASTQs)

| Read | 1PB mode | 3PY mode |
|------|---------|---------|
| R1 | 81 bp | ~81 bp |
| R2 | 150 bp | 150 bp |

R1 reads are shorter than R2 because R1 begins within the insert (after barcode/capture structure removal). R2 captures the full 150 bp from the other end of the fragment.
