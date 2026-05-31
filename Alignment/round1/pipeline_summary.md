# scATAC-seq Alignment Pipeline Summary

## Input
- Trimmed FASTQ: `Adapter_trimming/round1/step2_trim_fastq/{sample}_R1/R2.fq.gz`
- Filter tier: HD2 (Hamming distance ≤ 2), 6 samples (1PB/2PB/3PB/3PY/6PY/9PY)
- Reference genome: hg38 (`GRCh38/BWA/genome.fa`)

---

## Step 0 — Extract CB Map

**Script**: `extract_cb_map.py` (custom Python)  
**Input**: `R1_parser/round3/filter_HD2/{sample}_R1.fq.gz`  
**Output**: `step1_add_cb/{sample}_cb_map.tsv.gz`

Scans R1 FASTQ header comment field, extracts `bc1:Z:` and `bc2:Z:` tags, concatenates to 20 bp cell barcode (CB), and writes `read_name → CB + UMI` mapping TSV. No FASTQ rewriting — read-only scan.

---

## Step 1 — Build Bowtie2 Index (one-time)

**Software**: `bowtie2-build` v2.4.4  
**Command**:
```bash
bowtie2-build --threads 8 genome.fa /GRCh38/Bowtie2/genome
```
**Output**: `GRCh38/Bowtie2/genome.{1,2,3,4,rev.1,rev.2}.bt2`

---

## Step 2 — Alignment

**Software**: `Bowtie2` v2.4.4 + `samtools` v1.11  
**Output**: `step2_align/{sample}.sorted.bam`

```bash
bowtie2 --very-sensitive -X 2000 -p 8 \
  -x GRCh38/Bowtie2/genome \
  -1 {sample}_R1.fq.gz -2 {sample}_R2.fq.gz \
| samtools sort -@ 8 -O bam -o {sample}.sorted.bam
```

| Parameter | Value | Reason |
|-----------|-------|--------|
| `--very-sensitive` | — | Highest sensitivity mode for ATAC-seq |
| `-X 2000` | 2000 bp | Allow large insert sizes (mono/di/tri-nucleosome fragments) |
| `-p 8` | 8 threads | Parallel alignment |

---

## Step 3 — Add CB/UMI Tags to BAM

**Script**: `add_cb_tag.py` (custom Python, uses `pysam` v0.24.0)  
**Output**: `step3_tag_bam/{sample}.tagged.bam`

Loads `cb_map.tsv.gz` into memory, looks up each read by QNAME, and writes SAM auxiliary tags:
- `CB:Z:` — 20 bp cell barcode (BC1 + BC2)
- `UB:Z:` — 2 bp UMI

---

## Step 4 — Filter BAM

**Software**: `samtools` v1.11  
**Output**: `step4_filter/{sample}.filtered.bam`

```bash
samtools view -@ 8 -b -f 0x2 -F 0xC -q 30 \
  {sample}.tagged.bam {chr1..chrY} \
  -o {sample}.filtered.bam
```

| Filter | Flag/Option | Description |
|--------|-------------|-------------|
| Proper pairs only | `-f 0x2` | Keep only concordantly mapped pairs |
| Remove unmapped | `-F 0x4` | Remove unmapped reads |
| Remove mate unmapped | `-F 0x8` | Remove reads with unmapped mate |
| MAPQ ≥ 30 | `-q 30` | High-confidence alignments only |
| Remove chrM | region list | Exclude mitochondrial reads |

---

## Step 5 — Deduplicate per Cell Barcode

**Software**: `Picard` MarkDuplicates (via `java -jar picard.jar`)  
**Output**: `step5_dedup/{sample}.dedup.bam`, `{sample}.dup_metrics.txt`

```bash
java -jar picard.jar MarkDuplicates \
  -I {sample}.filtered.bam \
  -O {sample}.dedup.bam \
  -M {sample}.dup_metrics.txt \
  --BARCODE_TAG CB \
  --REMOVE_DUPLICATES true
```

`--BARCODE_TAG CB` ensures deduplication is performed **within each cell independently**, not globally. This prevents reads from different cells at the same genomic position from being incorrectly removed.

---

## Step 6 — Cell Calling (TODO)

Count reads per CB → knee plot → whitelist of valid cells.

---

## Step 7 — Generate Fragment File (TODO)

**Script**: `make_fragment_file.py` (custom Python, uses `pysam`)  
**Output**: `step7_fragments/{sample}_fragments.tsv.gz` + `.tbi`

Standard 5-column fragment file (0-based, half-open):
```
chr  start  end  barcode  count
```
Compatible with ArchR, Signac, and SnapATAC2.

---

## Software Versions

| Software | Version | Source |
|----------|---------|--------|
| Bowtie2 | 2.4.4 | `~/bowtie2-2.4.4-source/` |
| samtools | 1.11 | system |
| Picard | — | `~/picard.jar` |
| pysam | 0.24.0 | pip, `droplet_single_cell` conda env |
| Python | 3.14.5 | `micromamba droplet_single_cell` env |
