# Alignment QC Report

## Summary Table (vs GRCh38)

| Sample | Input pairs | Bowtie2 → hg38 | BWA → hg38 | After filter | After dedup | Dup rate |
|--------|------------|----------------|-----------|-------------|-------------|----------|
| 1PB | 19,273,973 | 1.55% | 17.35% | 156,323 | 15,260 | 90.2% |
| 2PB | 26,444,544 | 3.86% | 21.38% | 449,355 | 34,104 | 92.4% |
| 3PB | 40,035,523 | 4.72% | 22.51% | 1,553,424 | 115,566 | 92.6% |
| 3PY | 37,422,718 | 34.27% | 58.03% | 7,713,769 | 562,406 | 92.7% |
| 6PY | 10,485,405 | 38.69% | 59.85% | 1,673,593 | 215,376 | 87.1% |
| 9PY | 12,780,806 | 24.07% | 45.45% | 1,244,050 | 179,530 | 85.6% |

> BWA MEM reports inflated mapping rates due to soft-clipping (local alignment). Bowtie2 end-to-end is more conservative and accurate.  
> "After filter" = MAPQ ≥ 30 + proper pairs + no chrM. "Dup rate" = (filtered − dedup) / filtered.

---

## Mouse Genome Alignment (GRCm39) — Root Cause Analysis

All 6 samples were additionally aligned to the mouse reference genome (GRCm39/mm39) using both BWA MEM and Bowtie2 (end-to-end).

### Results

| Sample | Bowtie2 → hg38 | Bowtie2 → mm39 | BWA → mm39 | BWA properly paired (mm39) |
|--------|----------------|----------------|-----------|--------------------------|
| 1PB | 1.55% | **96.48%** | **99.33%** | 98.32% |
| 2PB | 3.86% | **93.74%** | **99.09%** | 98.09% |
| 3PB | 4.72% | **93.81%** | **98.38%** | 97.32% |
| 3PY | 34.27% | **53.19%** | — | — |
| 6PY | 38.69% | **46.91%** | — | — |
| 9PY | 24.07% | **61.69%** | — | — |

### Interpretation

**PB samples are almost entirely mouse DNA.** With >93% Bowtie2 mapping and >98% BWA mapping to GRCm39, the PB libraries contain essentially no human genomic signal. This explains the near-zero alignment to GRCh38 (1.55–4.72%).

**PY samples also contain substantial mouse DNA** (47–62% Bowtie2 mapping to GRCm39). Given that PY aligns to hg38 at 24–39%, the total mouse + human signal accounts for ~75–90% of reads, with the remainder being ambiguous (reads mapping to both, low-quality, or from other sources).

---

## BWA Soft-Clipping Analysis

BWA MEM's local alignment mode produces misleadingly high mapping rates for PB samples. Analysis of CIGAR strings on mapped reads reveals:

| Sample | Full match (no clip) | Heavy soft-clip (>30%) |
|--------|---------------------|----------------------|
| 1PB (→ hg38) | 8.6% | **84.4%** |
| 3PY (→ hg38) | 58.9% | 38.2% |

For 1PB, 84% of BWA "mapped" reads are only partially aligned — BWA finds weak partial matches in the human genome for what are actually mouse reads. The true human mapping rate for PB is ~1.5%, consistent with Bowtie2.

---

## Key Findings

### 1. PB samples are mouse DNA — not human ATAC-seq signal

**Root cause (revised):** Earlier analysis incorrectly attributed low alignment to repetitive/heterochromatic sequences. Mouse genome alignment conclusively shows PB DNA originates from mouse cells. Likely causes:

- K562 cells were co-cultured on mouse feeder cells at time of PB ATAC-seq preparation
- Mouse feeder cell DNA was co-extracted with K562 DNA
- Tn5 tagmentation captured predominantly mouse chromatin

**Implication:** PB samples cannot serve as a bulk ATAC-seq reference for K562. These data are not usable for human chromatin accessibility analysis.

### 2. PY samples contain a mixed human + mouse signal

PY (new technology) shows 47–62% alignment to mouse and 24–39% alignment to human. The droplet-based encapsulation may partially exclude extracellular mouse DNA, but a significant fraction of PY reads is still from mouse origin. This needs to be addressed before downstream analysis.

### 3. BWA vs Bowtie2 discrepancy

BWA MEM consistently reports 2–15× higher mapping rates than Bowtie2 for all samples due to local alignment soft-clipping. For cross-sample comparisons, Bowtie2 end-to-end rates are more reliable. BWA-mapped BAMs should be used cautiously for any sample with suspected mixed-species content.

### 4. High duplication rate across all samples

All samples show 85–93% duplication after filtering. For PY (single-cell), this is expected — each cell contributes few unique fragments. For PB, duplication analysis is moot given the wrong-species origin of the reads.

---

## Read Retention Through Pipeline

| Sample | Input pairs | Bowtie2-aligned | After filter | After dedup |
|--------|-------------|----------------|-------------|-------------|
| 1PB | 19.3M | 192K (1.0%) | 156K (0.81%) | 15K (0.08%) |
| 2PB | 26.4M | 1.0M (3.9%) | 449K (1.70%) | 34K (0.13%) |
| 3PB | 40.0M | 1.9M (4.7%) | 1.55M (3.88%) | 116K (0.29%) |
| 3PY | 37.4M | 12.8M (34.3%) | 7.7M (20.6%) | 562K (1.50%) |
| 6PY | 10.5M | 4.1M (38.7%) | 1.7M (15.9%) | 215K (2.05%) |
| 9PY | 12.8M | 3.1M (24.1%) | 1.2M (9.74%) | 180K (1.40%) |

---

## Read Structure

| Read | Length after trimming | Notes |
|------|-----------------------|-------|
| R1 | ~81 bp | Barcode + linker removed by R1_parser; Tn5 ME trimmed by cutadapt |
| R2 | ~150 bp | Full genomic read; Tn5 ME trimmed if short insert |

Both PB and PY share identical read structure after preprocessing — the species-of-origin difference is biological, not a pipeline artifact.
