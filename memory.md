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

## Next Step

Barcode and UMI extraction from step4 filtered reads, using `cs:i` tag as reference:
- BC1: pos 1–10 relative to read start
- W1: pos 11–15 (TCGAG, can validate)
- BC2: pos 16–25
- UMI_2N: pos 26–27
- capture start: `cs:i` value
