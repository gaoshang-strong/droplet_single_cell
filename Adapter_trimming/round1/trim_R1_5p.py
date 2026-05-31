"""
Trim the 5' end of R1 reads at the capture-sequence boundary.

Each filtered R1 read has the structure:
  [BC1 10bp][W1 5bp][BC2 10bp][10bp][CAPTURE_SEQ 15bp][INSERT ...]

The header carries cs:i:<pos> (1-based start of capture seq).
Trim point (0-based) = cs_val - 1 + len(CAPTURE_SEQ) = cs_val + 14

Input : filter_HD2/{key}_R1.fq.gz
Output: step1_trim_fastq/{key}_R1.fq.gz
"""
import os
import re
import subprocess
import multiprocessing

sys_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys; sys.path.insert(0, os.path.join(sys_path, "R1_parser", "round3"))
from config import SAMPLES, OUTPUT_DIR, CAPTURE_SEQ

FILTER_DIR = os.path.join(OUTPUT_DIR, "filter_HD2")
OUT_DIR    = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "step1_trim_fastq"
)

CS_RE    = re.compile(r'cs:i:(\d+)')
CS_LEN   = len(CAPTURE_SEQ)   # 15


def trim_sample(key: str):
    label   = SAMPLES[key]["label"]
    in_r1   = os.path.join(FILTER_DIR, f"{key}_R1.fq.gz")
    out_r1  = os.path.join(OUT_DIR,    f"{key}_R1.fq.gz")

    if not os.path.exists(in_r1):
        print(f"[{label}] missing {in_r1}, skipping", flush=True)
        return

    print(f"[{label}] trimming ...", flush=True)

    reader = subprocess.Popen(
        ["pigz", "-dc", "-p", "2", in_r1],
        stdout=subprocess.PIPE,
    )
    writer = subprocess.Popen(
        ["pigz", "-c", "-p", "2"],
        stdin=subprocess.PIPE,
        stdout=open(out_r1, "wb"),
    )

    total = trimmed = skipped = 0

    while True:
        header = reader.stdout.readline()
        if not header:
            break
        seq  = reader.stdout.readline()
        plus = reader.stdout.readline()
        qual = reader.stdout.readline()

        total += 1
        m = CS_RE.search(header.decode())
        if not m:
            skipped += 1
            writer.stdin.write(header + seq + plus + qual)
            continue

        # 0-based start of first INSERT base = cs_val - 1 + CS_LEN
        cut = int(m.group(1)) - 1 + CS_LEN

        seq_str  = seq.decode()
        qual_str = qual.decode()

        trimmed_seq  = seq_str[cut:].rstrip("\n")
        trimmed_qual = qual_str[cut:].rstrip("\n")

        if not trimmed_seq:
            skipped += 1
            continue

        trimmed += 1
        writer.stdin.write(
            header
            + (trimmed_seq  + "\n").encode()
            + plus
            + (trimmed_qual + "\n").encode()
        )

    reader.stdout.close()
    reader.wait()
    writer.stdin.close()
    writer.wait()

    print(
        f"[{label}] done  total={total:,}  trimmed={trimmed:,}  skipped={skipped:,}",
        flush=True,
    )


if __name__ == "__main__":
    os.makedirs(OUT_DIR, exist_ok=True)

    keys = list(SAMPLES.keys())
    print(f"Trimming {len(keys)} samples in parallel ...\n")

    with multiprocessing.Pool(processes=len(keys)) as pool:
        pool.map(trim_sample, keys)

    print("\nDone.")
