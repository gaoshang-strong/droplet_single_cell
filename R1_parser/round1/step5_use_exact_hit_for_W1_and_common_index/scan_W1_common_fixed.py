#!/usr/bin/env python3
"""
Exact hit search for W1 (TCGAG) and common_fixed (TAAGGCGA) in the region
before the capture sequence in step4-filtered R1 reads.

The capture start position is read from the cs:i tag in each read name.
Search region for each read: pos 1 to (cs_pos - 1), 1-based.
Position distribution is reported per sample.
"""
import os
import glob
import re
import subprocess
import multiprocessing

STEP4_DIR  = "/ShangGaoAIProjects/ZhangJW/R1_parser/round1/step4_filter_reads_with_anchored_capture_seq"
OUTPUT_DIR = "/ShangGaoAIProjects/ZhangJW/R1_parser/round1/step5_use_exact_hit_for_W1_and_common_index"

W1           = "TCGAG"
COMMON_FIXED = "TAAGGCGA"

SAMPLE_LABELS = {
    "260430R-S-XY-3PY": "3PY",
    "260430R-S-XY-1PB": "1PB",
    "260430R-S-XY-2PB": "2PB",
    "260430R-S-XY-6PY": "6PY",
    "260430R-S-XY-9PY": "9PY",
    "260430R-S-XY-3PB": "3PB",
}

CS_RE = re.compile(r'cs:i:(\d+)')

def extract_key(path):
    for key in SAMPLE_LABELS:
        if key in path:
            return key
    return None

def write_dist(path, pos_count, total_hit):
    with open(path, "w") as f:
        f.write("pos\tcount\trate_among_hits\n")
        for pos in sorted(pos_count):
            c = pos_count[pos]
            f.write(f"{pos}\t{c}\t{c/total_hit:.6f}\n")

def write_summary(path, seq, total, hit):
    with open(path, "w") as f:
        f.write(f"sequence\t{seq}\n")
        f.write(f"total_reads\t{total}\n")
        f.write(f"exact_hit\t{hit}\n")
        f.write(f"exact_hit_rate\t{hit/total:.6f}\n")

def process_sample(r1_path):
    key   = extract_key(r1_path)
    label = SAMPLE_LABELS[key]
    print(f"[{label}] start", flush=True)

    total = 0
    w1_hit = cf_hit = 0
    w1_pos_count = {}
    cf_pos_count = {}
    cs_pos = None

    proc = subprocess.Popen(["pigz", "-dc", "-p", "2", r1_path], stdout=subprocess.PIPE)
    line_num = 0
    for raw in proc.stdout:
        line_num += 1
        mod = line_num % 4

        if mod == 1:                          # header
            m = CS_RE.search(raw.decode())
            cs_pos = int(m.group(1)) if m else None

        elif mod == 2:                        # sequence
            if cs_pos is None:
                continue
            seq    = raw.decode().rstrip()
            total += 1
            # 0-based: search only before capture start
            region = seq[:cs_pos - 1]

            # W1
            p = region.find(W1)
            if p >= 0:
                w1_hit += 1
                w1_pos_count[p + 1] = w1_pos_count.get(p + 1, 0) + 1

            # common_fixed
            p = region.find(COMMON_FIXED)
            if p >= 0:
                cf_hit += 1
                cf_pos_count[p + 1] = cf_pos_count.get(p + 1, 0) + 1

    proc.wait()

    write_summary(os.path.join(OUTPUT_DIR, f"{key}_W1_summary.tsv"),           W1,           total, w1_hit)
    write_summary(os.path.join(OUTPUT_DIR, f"{key}_common_fixed_summary.tsv"), COMMON_FIXED, total, cf_hit)
    write_dist(os.path.join(OUTPUT_DIR, f"{key}_W1_position_distribution.tsv"),           w1_pos_count, w1_hit)
    write_dist(os.path.join(OUTPUT_DIR, f"{key}_common_fixed_position_distribution.tsv"), cf_pos_count, cf_hit)

    print(
        f"[{label}] done  total={total:,}"
        f"  W1_hit={w1_hit:,} ({w1_hit/total*100:.2f}%)"
        f"  CF_hit={cf_hit:,} ({cf_hit/total*100:.2f}%)",
        flush=True
    )
    return {"label": label, "key": key, "total": total, "w1_hit": w1_hit, "cf_hit": cf_hit}


if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    r1_files = sorted(glob.glob(os.path.join(STEP4_DIR, "*_filtered_R1.fq.gz")))

    n_workers = min(len(r1_files), multiprocessing.cpu_count())
    print(f"Processing {len(r1_files)} samples with {n_workers} workers ...\n")

    with multiprocessing.Pool(processes=n_workers) as pool:
        results = pool.map(process_sample, r1_files)

    # Combined summary
    summary_path = os.path.join(OUTPUT_DIR, "all_samples_summary.tsv")
    with open(summary_path, "w") as f:
        f.write("sample\ttotal\tW1_hit\tW1_hit_rate\tcommon_fixed_hit\tcommon_fixed_hit_rate\n")
        for r in results:
            t = r["total"]
            f.write(
                f"{r['label']}\t{t}\t{r['w1_hit']}\t{r['w1_hit']/t:.6f}\t"
                f"{r['cf_hit']}\t{r['cf_hit']/t:.6f}\n"
            )

    print(f"\nDone. Summary -> {summary_path}")
