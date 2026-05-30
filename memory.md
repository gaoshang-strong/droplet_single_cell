# Project Memory — ZhangJW Droplet Single Cell

## Samples

6 samples, paired-end R1/R2 fastq.gz, stored under `data/`:

| Sample | Batch | Quality |
|--------|-------|---------|
| 3PB | PB | High (pass rate ~99%) |
| 1PB | PB | High (~95%) |
| 2PB | PB | High (~99%) |
| 3PY | PY | Medium (~93%) |
| 6PY | PY | Low (~82%) |
| 9PY | PY | Low (~82%) |

---

## R1 Read Structure (1-based, closed)

```
BC1 (1-10) | W1 (11-15) | BC2 (16-25) | UMI_2N (26-27) | common_fixed (28-35) | capture (36-50) | insert (51-)
```

| Segment | Pos | Seq |
|---------|-----|-----|
| BC1 | 1–10 | random, 10 bp droplet barcode part 1 |
| W1 | 11–15 | TCGAG (fixed linker) |
| BC2 | 16–25 | random, 10 bp droplet barcode part 2 |
| UMI_2N | 26–27 | random 2N |
| common_fixed | 28–35 | **TAAGGCGA for PB; different sequence for PY** |
| capture | 36–50 | CACCGTCTCCGCCTC (15 bp, consistent across all samples) |
| insert | 51– | DNA insert |

---

## Key Conclusions

### Anchor strategy
- Use **15 bp capture seq `CACCGTCTCCGCCTC`** as the anchor — consistent across all 6 samples
- **Do NOT use common_fixed (TAAGGCGA) in anchor** — PB and PY batches have different 8 bp sequences at pos 28-35, confirmed by round2 23 bp anchor experiment

### Round2 23 bp anchor experiment (abandoned)
- Tested `TAAGGCGACACCGTCTCCGCCTC` (23 bp)
- Hit rates: PB 58–68%, PY 26–37%
- The ~30% gap between batches confirms PY samples have a different common_fixed sequence

### Capture exact hit rates (15 bp anchor, round1 step1)
- 3PB: 95.2%, 1PB: 86.4%, 2PB: 90.0%
- 3PY: 82.7%, 6PY: 60.9%, 9PY: 63.0%
- All samples: dominant peak at **pos 36**, as expected

### W1 (TCGAG) is a reliable structural anchor
- Exact hit rate: 92–98% across all samples
- Peak: pos 11 in all samples
- High W1 hit rate validates that the overall read structure is intact and consistent

### Hamming rescue (round1 step2/step4)
- For reads without exact capture hit, scan pos 33–60 with 15 bp sliding window, Hamming ≤ 3
- Take **minimum Hamming window** — do NOT add CAC prefix filter (would incorrectly skip true anchor windows where error falls in first 3 bases, leading to wrong cs_pos)
- Rescued reads are tagged with `cs:i:{pos} mt:Z:hamming` in the R1 read name

### Filtered read output (round1 step4)
- Output: `*_filtered_R1.fq.gz` and `*_filtered_R2.fq.gz`
- Read name tag: `cs:i:{1-based capture start} mt:Z:{exact|hamming}`
- Pass rates: 3PB/2PB ~99%, 1PB ~95%, 3PY ~93%, 6PY/9PY ~82%
- Files on disk at `R1_parser/round1/step4_filter_reads_with_anchored_capture_seq/` (not in git — too large)

---

## Environment

- micromamba env: `droplet_single_cell`
- micromamba: `/home/sgao30/micromamba/bin/micromamba`
- Python: `/home/sgao30/micromamba/envs/droplet_single_cell/bin/python`
- GitHub: `gaoshang-strong/droplet_single_cell`

---

## Step 6 — W1 exact hit filter

- Input: step4 filtered reads; output: `*_W1_R1/R2.fq.gz` with `w1:i:{pos}` appended to R1 header
- Pass rates: PB 96–98%, PY 92–94%
- R1 read name format after step6: `... cs:i:{n} mt:Z:{exact|hamming} w1:i:{n}`

## Step 7 — Barcode / UMI extraction

Input: step6 `*_W1_R1.fq.gz`  
Output: `*_bc_umi.tsv.gz` per sample (one row per read)

Columns: `read_id, w1_pos, cs_pos, mt, bc1, bc1_len, bc2, umi, gap_len, gap_seq`

0-based Python slice math (w1_pos is 1-based):
- BC1: `seq[max(0, w1_pos-11) : w1_pos-1]`  (10 bp, may be shorter if w1_pos < 11)
- BC2: `seq[w1_pos+4 : w1_pos+14]`           (10 bp)
- UMI: `seq[w1_pos+14 : w1_pos+16]`          (2 bp)
- gap: `seq[w1_pos+16 : cs_pos-1]`           (common_fixed region, normally 8 bp)

Canonical check (w1_pos=11, cs_pos=36): gap = `TAAGGCGA` for PB ✓

## Next Step

Analyze step7 TSV output:
- gap_len distribution (expect 8bp dominant)
- gap_seq exact/Hamming match to common_fixed (TAAGGCGA for PB; empirically derived for PY)
- bc1_len distribution (flag truncated BC1 reads)
