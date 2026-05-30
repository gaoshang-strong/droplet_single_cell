#!/usr/bin/env python3
"""
Filter paired-end reads based on anchor capture sequence in R1.
Pass criteria (either):
  1. Exact match of anchor anywhere in R1
  2. Sliding window in pos 33-60 (1-based), length 15: min Hamming <= 3

Capture start position (1-based) and match type are written into the R1
read name comment field:
  @readname cs:i:{pos} mt:Z:{exact|hamming}

Samples are processed in parallel (one process per sample).
Outputs filtered R1 + R2 fastq.gz and a filter_summary.tsv.
"""
import os
import glob
import subprocess
import gzip
import multiprocessing

DATA_DIR    = "/ShangGaoAIProjects/ZhangJW/data"
OUTPUT_DIR  = "/ShangGaoAIProjects/ZhangJW/R1_parser/round1/step4_filter_reads_with_anchored_capture_seq"
ANCHOR      = "CACCGTCTCCGCCTC"
ANCHOR_LEN  = len(ANCHOR)   # 15
MAX_HAMMING = 3
WIN_STARTS  = list(range(32, 46))  # 0-based; 1-based pos 33-46, window end ≤ 60

SAMPLE_LABELS = {
    "260430R-S-XY-3PY": "3PY",
    "260430R-S-XY-1PB": "1PB",
    "260430R-S-XY-2PB": "2PB",
    "260430R-S-XY-6PY": "6PY",
    "260430R-S-XY-9PY": "9PY",
    "260430R-S-XY-3PB": "3PB",
}

def extract_key(path):
    for key in SAMPLE_LABELS:
        if key in path:
            return key
    return None

def hamming(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def annotate_header(header_bytes, cs_pos, mt):
    """Append cs and mt tags to the read name comment field (1-based pos)."""
    header = header_bytes.decode().rstrip()
    return (f"{header} cs:i:{cs_pos} mt:Z:{mt}\n").encode()

def process_sample(r1_path):
    key   = extract_key(r1_path)
    label = SAMPLE_LABELS[key]
    r2_path = r1_path.replace("_R1.fq.gz", "_R2.fq.gz")
    out_r1  = os.path.join(OUTPUT_DIR, f"{key}_filtered_R1.fq.gz")
    out_r2  = os.path.join(OUTPUT_DIR, f"{key}_filtered_R2.fq.gz")

    print(f"[{label}] start", flush=True)

    total = passed = exact = hamming_rescued = 0

    p_r1 = subprocess.Popen(["pigz", "-dc", "-p", "2", r1_path], stdout=subprocess.PIPE)
    p_r2 = subprocess.Popen(["pigz", "-dc", "-p", "2", r2_path], stdout=subprocess.PIPE)

    with gzip.open(out_r1, "wb", compresslevel=1) as f1, \
         gzip.open(out_r2, "wb", compresslevel=1) as f2:
        while True:
            r1_block = [p_r1.stdout.readline() for _ in range(4)]
            r2_block = [p_r2.stdout.readline() for _ in range(4)]
            if not r1_block[0]:
                break
            total += 1
            seq = r1_block[1].decode().rstrip()

            # --- exact hit ---
            if ANCHOR in seq:
                cs_pos = seq.index(ANCHOR) + 1  # 1-based
                exact += 1
                passed += 1
                r1_block[0] = annotate_header(r1_block[0], cs_pos, "exact")
                for line in r1_block: f1.write(line)
                for line in r2_block: f2.write(line)
                continue

            # --- hamming rescue ---
            best_pos  = None
            best_dist = MAX_HAMMING + 1
            for start in WIN_STARTS:
                window = seq[start:start + ANCHOR_LEN]
                if len(window) < ANCHOR_LEN:
                    continue
                d = hamming(window, ANCHOR)
                if d < best_dist:
                    best_dist = d
                    best_pos  = start + 1  # 1-based
                    if d == 0:
                        break

            if best_dist <= MAX_HAMMING:
                hamming_rescued += 1
                passed += 1
                r1_block[0] = annotate_header(r1_block[0], best_pos, "hamming")
                for line in r1_block: f1.write(line)
                for line in r2_block: f2.write(line)

    p_r1.wait()
    p_r2.wait()

    failed = total - passed
    print(f"[{label}] done  total={total:,}  passed={passed:,} ({passed/total*100:.2f}%)"
          f"  exact={exact:,}  hamming_rescued={hamming_rescued:,}  failed={failed:,}", flush=True)

    return {
        "label": label, "key": key,
        "total": total, "passed": passed,
        "exact": exact, "hamming_rescued": hamming_rescued,
        "failed": failed,
    }


if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    r1_files = sorted(glob.glob(os.path.join(DATA_DIR, "**/*_R1.fq.gz"), recursive=True))
    r1_files = [f for f in r1_files if extract_key(f) is not None]

    n_workers = min(len(r1_files), multiprocessing.cpu_count())
    print(f"Processing {len(r1_files)} samples with {n_workers} parallel workers ...\n")

    with multiprocessing.Pool(processes=n_workers) as pool:
        results = pool.map(process_sample, r1_files)

    # Write summary TSV
    summary_path = os.path.join(OUTPUT_DIR, "filter_summary.tsv")
    with open(summary_path, "w") as f:
        f.write("sample\ttotal\tpassed\tpassed_pct\texact_hit\thamming_rescued\tfailed\tfailed_pct\n")
        for r in results:
            t = r["total"]
            f.write(
                f"{r['label']}\t{t}\t{r['passed']}\t{r['passed']/t*100:.4f}\t"
                f"{r['exact']}\t{r['hamming_rescued']}\t{r['failed']}\t{r['failed']/t*100:.4f}\n"
            )

    print(f"\nAll done. Summary -> {summary_path}")
