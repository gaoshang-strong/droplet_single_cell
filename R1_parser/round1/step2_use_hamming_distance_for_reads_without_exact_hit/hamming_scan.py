#!/usr/bin/env python3
"""
For reads without an exact anchor hit, scan positions 33-60 (1-based)
with a sliding window of length 15 and report minimum Hamming distance distribution.
"""
import os
import glob
import subprocess

DATA_DIR = "/ShangGaoAIProjects/ZhangJW/data"
OUTPUT_DIR = "/ShangGaoAIProjects/ZhangJW/R1_parser/round1/step2_use_local_alignment_for_reads_without_exact_hit"
ANCHOR = "CACCGTCTCCGCCTC"
ANCHOR_LEN = len(ANCHOR)  # 15

# 1-based window: start 33..46, end 47..60 (window length 15, all within pos 33-60)
# 0-based: start index 32..45
WIN_STARTS = range(32, 46)

SAMPLE_LABELS = {
    "260430R-S-XY-3PY": "3PY",
    "260430R-S-XY-1PB": "1PB",
    "260430R-S-XY-2PB": "2PB",
    "260430R-S-XY-6PY": "6PY",
    "260430R-S-XY-9PY": "9PY",
    "260430R-S-XY-3PB": "3PB",
}

def extract_sample_key(filename):
    for key in SAMPLE_LABELS:
        if key in filename:
            return key
    return None

def hamming(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def min_hamming_in_window(seq):
    min_d = ANCHOR_LEN + 1
    for start in WIN_STARTS:
        window = seq[start:start + ANCHOR_LEN]
        if len(window) < ANCHOR_LEN:
            continue
        d = hamming(window, ANCHOR)
        if d < min_d:
            min_d = d
            if min_d == 0:
                break
    return min_d

os.makedirs(OUTPUT_DIR, exist_ok=True)

r1_files = sorted(glob.glob(os.path.join(DATA_DIR, "**/*_R1.fq.gz"), recursive=True))

all_results = []

for r1 in r1_files:
    key = extract_sample_key(r1)
    if key is None:
        continue
    label = SAMPLE_LABELS[key]
    print(f"Processing {label} ({os.path.basename(r1)}) ...")

    counts = {1: 0, 2: 0, 3: 0, ">3": 0}
    no_hit_total = 0

    proc = subprocess.Popen(["pigz", "-dc", "-p", "4", r1], stdout=subprocess.PIPE)
    line_num = 0
    for raw in proc.stdout:
        line_num += 1
        if line_num % 4 != 2:
            continue
        seq = raw.decode().rstrip()
        if ANCHOR in seq:
            continue
        no_hit_total += 1
        min_d = min_hamming_in_window(seq)
        if min_d == 1:
            counts[1] += 1
        elif min_d == 2:
            counts[2] += 1
        elif min_d == 3:
            counts[3] += 1
        else:
            counts[">3"] += 1
    proc.wait()

    all_results.append({"label": label, "key": key, "no_hit_total": no_hit_total, "counts": counts})

    out_file = os.path.join(OUTPUT_DIR, f"{key}_hamming_summary.tsv")
    with open(out_file, "w") as f:
        f.write("min_hamming\tcount\tproportion\n")
        for cat in [1, 2, 3, ">3"]:
            c = counts[cat]
            prop = c / no_hit_total if no_hit_total > 0 else 0
            f.write(f"{cat}\t{c}\t{prop:.6f}\n")
    print(f"  no_exact_hit_reads={no_hit_total:,}  ->  {os.path.basename(out_file)}")

# Combined summary across all samples
summary_file = os.path.join(OUTPUT_DIR, "all_samples_hamming_summary.tsv")
with open(summary_file, "w") as f:
    f.write("sample\tno_hit_reads\t"
            "hamming_1\thamming_1_pct\t"
            "hamming_2\thamming_2_pct\t"
            "hamming_3\thamming_3_pct\t"
            "hamming_gt3\thamming_gt3_pct\n")
    for r in all_results:
        c = r["counts"]
        n = r["no_hit_total"]
        def pct(x):
            return f"{x / n * 100:.2f}" if n > 0 else "0.00"
        f.write(
            f"{r['label']}\t{n}\t"
            f"{c[1]}\t{pct(c[1])}\t"
            f"{c[2]}\t{pct(c[2])}\t"
            f"{c[3]}\t{pct(c[3])}\t"
            f"{c['>3']}\t{pct(c['>3'])}\n"
        )

print(f"\nDone. Combined summary -> {summary_file}")
