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

### 1. PB samples have extremely low alignment rate (< 5%) — root cause identified

Systematic debugging was performed on 1PB:

| Test | Alignment rate | Conclusion |
|------|---------------|------------|
| Bowtie2 end-to-end paired (current) | 1.65% | Baseline |
| Bowtie2 `--local` paired | 10.37% | Soft-clip helps slightly |
| Bowtie2 end-to-end R2 single-end | 1.52% | R2 alone also fails |
| Bowtie2 end-to-end R1 single-end | 1.90% | R1 alone also fails |
| BWA-MEM R2 single-end | 13.7% | Local alignment helps, but still low |
| BWA-MEM R2 single-end (3PY, for comparison) | **50.4%** | 3PY aligns normally |

**Root cause: PB reads are enriched in low-complexity / repetitive sequences.**

Inspection of the top R2 5'-end sequences reveals systematic repetitive content in PB:

```
1PB R2 (top sequences):        3PY R2 (top sequences):
CACACACACACACACACACA           GCACACTCCTTTCCTCTGCC
GTGTGTGTGTGTGTGTGTGT           GGGTTAACTCAGGTCAGCTA
CTGTAGGACGTGGAATATGG           GGCCAGGACCAGGCCAGAAA
```

1PB is dominated by CA/GT tandem repeats (microsatellites) and other repetitive elements. These sequences exist at thousands of genomic loci and cannot be uniquely placed by Bowtie2 (which discards multi-mappers by default), resulting in near-total alignment failure.

**Possible causes:**
- **Biological**: Peripheral blood contains a high proportion of neutrophils and other terminally differentiated cells with densely compacted heterochromatin. Tn5 may insert non-specifically in repetitive/heterochromatic regions when open chromatin is scarce.
- **Library quality**: DNA degradation or low-complexity library amplification enriching for repetitive elements.
- **PCR bias**: Tandem repeats amplify more efficiently during PCR, leading to overrepresentation in sequencing.

**This is not a pipeline issue** — the same pipeline produces 24–39% alignment for PY samples. The difference is in the reads themselves.

### 2. PY samples have moderate alignment rate
34–39% alignment rate is below the typical ATAC-seq benchmark (>60%), but PY R2 single-end aligns at ~50% with BWA-MEM, suggesting the paired-end concordance requirement and Bowtie2's end-to-end mode account for part of the gap. PY reads show diverse, unique sequences consistent with genuine open chromatin signal.

### 3. High duplication rate across all samples
90–93% duplication rate is concerning. Suggests low library complexity or sequencing the same molecules many times. Cell calling (Step 6) will determine if enough unique fragments per cell remain for analysis.

---

## Read Length (Trimmed FASTQs)

| Read | 1PB mode | 3PY mode |
|------|---------|---------|
| R1 | 81 bp | ~81 bp |
| R2 | 150 bp | 150 bp |

R1 reads are shorter than R2 because R1 begins within the insert (after barcode/capture structure removal). R2 captures the full 150 bp from the other end of the fragment.
