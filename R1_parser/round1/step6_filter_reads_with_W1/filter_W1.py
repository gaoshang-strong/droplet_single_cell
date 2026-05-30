#!/usr/bin/env python3
"""
Filter step4-filtered reads by W1 (TCGAG) exact hit in the region before
the capture sequence.  Reads without an exact W1 match are discarded.

Input : step4 *_filtered_R1/R2.fq.gz  (read names carry cs:i and mt:Z tags)
Output: *_W1_R1.fq.gz / *_W1_R2.fq.gz with w1:i:{pos} appended to R1 header

Tag format (1-based, same convention as cs:i):
  w1:i:{start_pos}
"""
import os
import glob
import re
import subprocess
import gzip
import multiprocessing

STEP4_DIR  = "/ShangGaoAIProjects/ZhangJW/R1_parser/round1/step4_filter_reads_with_anchored_capture_seq"
OUTPUT_DIR = "/ShangGaoAIProjects/ZhangJW/R1_parser/round1/step6_filter_reads_with_W1"

W1     = "TCGAG"
CS_RE  = re.compile(rb'cs:i:(\d+)')

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

def annotate_header(header_bytes, w1_pos):
    """Append w1:i tag to the existing R1 header (1-based pos)."""
    header = header_bytes.rstrip(b'\n')
    return header + f" w1:i:{w1_pos}\n".encode()

def process_sample(r1_path):
    key    = extract_key(r1_path)
    label  = SAMPLE_LABELS[key]
    r2_path = r1_path.replace("_filtered_R1.fq.gz", "_filtered_R2.fq.gz")
    out_r1  = os.path.join(OUTPUT_DIR, f"{key}_W1_R1.fq.gz")
    out_r2  = os.path.join(OUTPUT_DIR, f"{key}_W1_R2.fq.gz")

    print(f"[{label}] start", flush=True)

    total = passed = 0

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

            # parse cs:i from header
            m = CS_RE.search(r1_block[0])
            if m is None:
                continue
            cs_pos = int(m.group(1))   # 1-based capture start

            # search W1 only before capture start (0-based: 0 .. cs_pos-2)
            seq    = r1_block[1].decode().rstrip()
            region = seq[:cs_pos - 1]
            p      = region.find(W1)
            if p < 0:
                continue

            w1_pos = p + 1   # convert to 1-based
            passed += 1
            r1_block[0] = annotate_header(r1_block[0], w1_pos)
            for line in r1_block: f1.write(line)
            for line in r2_block: f2.write(line)

    p_r1.wait()
    p_r2.wait()

    failed = total - passed
    print(
        f"[{label}] done  total={total:,}  passed={passed:,} ({passed/total*100:.2f}%)"
        f"  failed={failed:,}",
        flush=True,
    )
    return {
        "label": label, "key": key,
        "total": total, "passed": passed, "failed": failed,
    }


if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    r1_files = sorted(glob.glob(os.path.join(STEP4_DIR, "*_filtered_R1.fq.gz")))
    r1_files = [f for f in r1_files if extract_key(f) is not None]

    n_workers = min(len(r1_files), multiprocessing.cpu_count())
    print(f"Processing {len(r1_files)} samples with {n_workers} workers ...\n")

    with multiprocessing.Pool(processes=n_workers) as pool:
        results = pool.map(process_sample, r1_files)

    summary_path = os.path.join(OUTPUT_DIR, "W1_filter_summary.tsv")
    with open(summary_path, "w") as f:
        f.write("sample\ttotal\tpassed\tpassed_pct\tfailed\tfailed_pct\n")
        for r in results:
            t = r["total"]
            f.write(
                f"{r['label']}\t{t}\t{r['passed']}\t{r['passed']/t*100:.4f}\t"
                f"{r['failed']}\t{r['failed']/t*100:.4f}\n"
            )

    print(f"\nAll done. Summary -> {summary_path}")
