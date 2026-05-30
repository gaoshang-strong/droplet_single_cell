#!/usr/bin/env python3
"""
Filter step6 W1-filtered reads to keep only reads with gap_len == 8.

gap_len is defined as cs_pos - w1_pos - 17 (number of nt between UMI
end and capture start).  gap_len == 8 guarantees:
  - capture start is correctly identified (no Hamming rescue false positives)
  - common_fixed region is intact (no indels)

Strategy: stream step7 TSV (gap_len column) and step6 R1/R2 fastq in
lockstep — same read order, no memory-intensive set lookup needed.

Output: *_gap8_R1.fq.gz / *_gap8_R2.fq.gz
"""
import os
import glob
import subprocess
import gzip
import multiprocessing

STEP6_DIR  = "/ShangGaoAIProjects/ZhangJW/R1_parser/round1/step6_filter_reads_with_W1"
STEP7_DIR  = "/ShangGaoAIProjects/ZhangJW/R1_parser/round1/step7_extract_barcode_UMI"
OUTPUT_DIR = "/ShangGaoAIProjects/ZhangJW/R1_parser/round1/step8_further_filter_based_on_gap8"

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

def process_sample(tsv_path):
    key   = extract_key(tsv_path)
    label = SAMPLE_LABELS[key]
    r1_path = os.path.join(STEP6_DIR, f"{key}_W1_R1.fq.gz")
    r2_path = os.path.join(STEP6_DIR, f"{key}_W1_R2.fq.gz")
    out_r1  = os.path.join(OUTPUT_DIR, f"{key}_gap8_R1.fq.gz")
    out_r2  = os.path.join(OUTPUT_DIR, f"{key}_gap8_R2.fq.gz")

    print(f"[{label}] start", flush=True)

    total = passed = 0

    p_tsv = subprocess.Popen(["pigz", "-dc", "-p", "2", tsv_path],  stdout=subprocess.PIPE)
    p_r1  = subprocess.Popen(["pigz", "-dc", "-p", "2", r1_path],   stdout=subprocess.PIPE)
    p_r2  = subprocess.Popen(["pigz", "-dc", "-p", "2", r2_path],   stdout=subprocess.PIPE)

    p_tsv.stdout.readline()   # skip header

    with gzip.open(out_r1, "wb", compresslevel=1) as f1, \
         gzip.open(out_r2, "wb", compresslevel=1) as f2:
        while True:
            tsv_line = p_tsv.stdout.readline()
            r1_block = [p_r1.stdout.readline() for _ in range(4)]
            r2_block = [p_r2.stdout.readline() for _ in range(4)]

            if not tsv_line or not r1_block[0]:
                break
            total += 1

            parts    = tsv_line.decode().rstrip("\n").split("\t")
            gap_len  = int(parts[8])

            if gap_len == 8:
                passed += 1
                for line in r1_block: f1.write(line)
                for line in r2_block: f2.write(line)

    p_tsv.wait(); p_r1.wait(); p_r2.wait()

    failed = total - passed
    print(
        f"[{label}] done  total={total:,}  passed={passed:,} ({passed/total*100:.2f}%)"
        f"  failed={failed:,}",
        flush=True,
    )
    return {"label": label, "key": key, "total": total, "passed": passed, "failed": failed}


if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    tsv_files = sorted(glob.glob(os.path.join(STEP7_DIR, "*_bc_umi.tsv.gz")))
    tsv_files = [f for f in tsv_files if extract_key(f) is not None]

    n_workers = min(len(tsv_files), multiprocessing.cpu_count())
    print(f"Processing {len(tsv_files)} samples with {n_workers} workers ...\n")

    with multiprocessing.Pool(processes=n_workers) as pool:
        results = pool.map(process_sample, tsv_files)

    summary_path = os.path.join(OUTPUT_DIR, "gap8_filter_summary.tsv")
    with open(summary_path, "w") as f:
        f.write("sample\ttotal\tpassed\tpassed_pct\tfailed\tfailed_pct\n")
        for r in results:
            t = r["total"]
            f.write(
                f"{r['label']}\t{t}\t{r['passed']}\t{r['passed']/t*100:.4f}\t"
                f"{r['failed']}\t{r['failed']/t*100:.4f}\n"
            )

    print(f"\nAll done. Summary -> {summary_path}")
