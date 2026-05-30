#!/usr/bin/env python3
"""
Extract BC1, BC2, UMI_2N, and gap sequence from step6 W1-filtered R1 reads.

For each read the following are derived from w1:i and cs:i tags:
  BC1     : 10 nt immediately before W1 (may be < 10 if w1_pos < 11)
  BC2     : 10 nt immediately after W1
  UMI_2N  : 2 nt immediately after BC2
  gap_seq : nt between UMI end and capture start (common_fixed region)

All positions are 1-based in the tags; Python slicing uses 0-based indices.

Canonical structure (w1_pos=11, cs_pos=36):
  BC1  [0:10]   pos 1-10
  W1   [10:15]  pos 11-15
  BC2  [15:25]  pos 16-25
  UMI  [25:27]  pos 26-27
  gap  [27:35]  pos 28-35  (8 bp = common_fixed)
  cap  [35:]    pos 36+

Output: one <key>_bc_umi.tsv.gz per sample with columns:
  read_id, w1_pos, cs_pos, mt, bc1, bc1_len, bc2, umi, gap_len, gap_seq
"""
import os
import glob
import re
import subprocess
import gzip
import multiprocessing

STEP6_DIR  = "/ShangGaoAIProjects/ZhangJW/R1_parser/round1/step6_filter_reads_with_W1"
OUTPUT_DIR = "/ShangGaoAIProjects/ZhangJW/R1_parser/round1/step7_extract_barcode_UMI"

HEADER_RE = re.compile(
    rb'^@(\S+)\s.*cs:i:(\d+).*mt:Z:(\S+).*w1:i:(\d+)'
)

SAMPLE_LABELS = {
    "260430R-S-XY-3PY": "3PY",
    "260430R-S-XY-1PB": "1PB",
    "260430R-S-XY-2PB": "2PB",
    "260430R-S-XY-6PY": "6PY",
    "260430R-S-XY-9PY": "9PY",
    "260430R-S-XY-3PB": "3PB",
}

HEADER = "read_id\tw1_pos\tcs_pos\tmt\tbc1\tbc1_len\tbc2\tumi\tgap_len\tgap_seq\n"

def extract_key(path):
    for key in SAMPLE_LABELS:
        if key in path:
            return key
    return None

def process_sample(r1_path):
    key   = extract_key(r1_path)
    label = SAMPLE_LABELS[key]
    out   = os.path.join(OUTPUT_DIR, f"{key}_bc_umi.tsv.gz")

    print(f"[{label}] start", flush=True)

    total = written = 0

    proc = subprocess.Popen(["pigz", "-dc", "-p", "2", r1_path], stdout=subprocess.PIPE)

    with gzip.open(out, "wt", compresslevel=1) as fout:
        fout.write(HEADER)
        while True:
            header = proc.stdout.readline()
            seq    = proc.stdout.readline()
            proc.stdout.readline()   # +
            proc.stdout.readline()   # qual

            if not header:
                break
            total += 1

            m = HEADER_RE.match(header)
            if m is None:
                continue

            read_id = m.group(1).decode()
            cs_pos  = int(m.group(2))
            mt      = m.group(3).decode()
            w1_pos  = int(m.group(4))

            s = seq.decode().rstrip()

            # 0-based slicing; w1_pos is 1-based
            bc1_start = max(0, w1_pos - 11)   # w1_pos-1 - 10
            bc1_end   = w1_pos - 1
            bc1       = s[bc1_start : bc1_end]

            bc2_start = w1_pos + 4             # w1_pos-1 + 5
            bc2       = s[bc2_start : bc2_start + 10]

            umi_start = bc2_start + 10
            umi       = s[umi_start : umi_start + 2]

            gap_start = umi_start + 2
            gap_end   = cs_pos - 1             # 0-based exclusive = cs_pos-1
            gap_seq   = s[gap_start : gap_end]

            bc1_len  = len(bc1)
            gap_len  = len(gap_seq)

            fout.write(
                f"{read_id}\t{w1_pos}\t{cs_pos}\t{mt}\t"
                f"{bc1}\t{bc1_len}\t{bc2}\t{umi}\t"
                f"{gap_len}\t{gap_seq}\n"
            )
            written += 1

    proc.wait()
    print(f"[{label}] done  total={total:,}  written={written:,}", flush=True)
    return {"label": label, "key": key, "total": total, "written": written}


if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    r1_files = sorted(glob.glob(os.path.join(STEP6_DIR, "*_W1_R1.fq.gz")))
    r1_files = [f for f in r1_files if extract_key(f) is not None]

    n_workers = min(len(r1_files), multiprocessing.cpu_count())
    print(f"Processing {len(r1_files)} samples with {n_workers} workers ...\n")

    with multiprocessing.Pool(processes=n_workers) as pool:
        results = pool.map(process_sample, r1_files)

    summary_path = os.path.join(OUTPUT_DIR, "extract_summary.tsv")
    with open(summary_path, "w") as f:
        f.write("sample\ttotal_reads\twritten\n")
        for r in results:
            f.write(f"{r['label']}\t{r['total']}\t{r['written']}\n")

    print(f"\nAll done. Summary -> {summary_path}")
