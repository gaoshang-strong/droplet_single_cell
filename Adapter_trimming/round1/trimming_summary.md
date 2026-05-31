# Adapter Trimming Summary — Round 1

**Input:** `filter_HD2/` (HD≤2, W1 exact, gap=20)  
**Output:** `step2_trim_fastq/` (trimmed R1 + R2, ready for alignment)

---

## Pipeline

### Step 1 — 5′ Capture Sequence Trimming (`trim_R1_5p.py`)

Remove everything at the 5′ end of R1 up to and including the capture sequence.

**R1 structure before trimming:**
```
[BC1 10bp][W1 5bp][BC2 10bp][10bp][CAPTURE_SEQ 15bp] | INSERT ...
 ◄────────────────── removed ──────────────────────►   ◄─ kept ─►
```

- Trim point derived from `cs:i` tag in each read header (1-based capture seq start)
- Cut position (0-based): `cs_val - 1 + 15`
- All reads trimmed, skipped = 0 across all samples

| Sample | Input reads |
|--------|------------|
| 1PB | 19,304,120 |
| 2PB | 26,498,256 |
| 3PB | 40,101,821 |
| 3PY | 37,682,637 |
| 6PY | 10,561,326 |
| 9PY | 12,862,720 |

---

### Step 2 — Nextera ME Adapter Trimming (`trim_nextera.py`)

Remove Nextera Mosaic End (ME) sequences using cutadapt 5.2 (paired-end mode).

**Adapter sequences:**

| Position | Sequence |
|----------|----------|
| R1 5′ | `AGATGTGTATAAGAGACAG` |
| R1 3′ | `CTGTCTCTTATACACATCT` (ME RC, read-through) |
| R2 3′ | `CTGTCTCTTATACACATCT` (ME RC, read-through) |

**cutadapt parameters:**
```
-e 0.1  -O 10  --minimum-length 20  --pair-filter=any  -j 8
```

**Results:**

| Sample | Total pairs | R1 w/ adapter | R2 w/ adapter | Too short | Pairs kept |
|--------|------------|--------------|--------------|-----------|------------|
| 1PB | 19,304,120 | 17,859,542 (92.5%) | 755,161 (3.9%) | 30,147 (0.2%) | 19,273,973 (99.8%) |
| 2PB | 26,498,256 | 26,254,837 (99.1%) | 1,837,937 (6.9%) | 53,712 (0.2%) | 26,444,544 (99.8%) |
| 3PB | 40,101,821 | 39,501,037 (98.5%) | 1,847,426 (4.6%) | 66,298 (0.2%) | 40,035,523 (99.8%) |
| 3PY | 37,682,637 | 34,300,791 (91.0%) | 15,509,640 (41.2%) | 259,919 (0.7%) | 37,422,718 (99.3%) |
| 6PY | 10,561,326 | 7,709,170 (73.0%) | 4,053,092 (38.4%) | 75,921 (0.7%) | 10,485,405 (99.3%) |
| 9PY | 12,862,720 | 9,615,357 (74.8%) | 3,985,748 (31.0%) | 81,914 (0.6%) | 12,780,806 (99.4%) |

---

## Observations

- **R1 adapter rate:** All samples show high R1 ME adapter detection (73–99%), confirming the capture-seq trimming in Step 1 correctly exposed the insert-proximal ME sequence.
- **R2 adapter rate — PB vs PY:**
  - PB samples (1PB/2PB/3PB): R2 adapter rate 4–7%, indicating predominantly **long inserts** where R2 does not read through to the ME.
  - PY samples (3PY/6PY/9PY): R2 adapter rate 31–41%, indicating a larger proportion of **short inserts** with read-through.
- **Pairs too short:** <1% across all samples; minimal data loss.
- **Overall yield:** ≥99.3% of read pairs pass filters in all samples.

---

## Output Files

```
step1_trim_fastq/   {key}_R1.fq.gz                  (capture-seq trimmed R1)
step2_trim_fastq/   {key}_R1.fq.gz                  (final trimmed R1)
                    {key}_R2.fq.gz                  (final trimmed R2)
                    {key}_cutadapt.log               (full cutadapt report)
                    {key}_cutadapt.json              (machine-readable stats)
```
