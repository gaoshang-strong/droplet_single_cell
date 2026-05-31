"""
Step 2: Nextera ME adapter trimming with cutadapt (paired-end).

Input:
  R1: step1_trim_fastq/{key}_R1.fq.gz   (5' capture-seq already removed)
  R2: filter_HD2/{key}_R2.fq.gz          (original, untrimmed)

Adapter (Nextera Mosaic End, ME):
  R1 5' : AGATGTGTATAAGAGACAG
  R1 3' : CTGTCTCTTATACACATCT   (ME RC, for read-through)
  R2 3' : CTGTCTCTTATACACATCT

Output:
  step2_trim_fastq/{key}_R1.fq.gz
  step2_trim_fastq/{key}_R2.fq.gz
  step2_trim_fastq/{key}_cutadapt.log
  step2_trim_fastq/{key}_cutadapt.json
"""
import os
import shutil
import subprocess
import multiprocessing

sys_path = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "R1_parser", "round3",
)
import sys; sys.path.insert(0, sys_path)
from config import SAMPLES, OUTPUT_DIR

HERE       = os.path.dirname(os.path.abspath(__file__))
STEP1_DIR  = os.path.join(HERE, "step1_trim_fastq")
FILTER_DIR = os.path.join(OUTPUT_DIR, "filter_HD2")
OUT_DIR    = os.path.join(HERE, "step2_trim_fastq")

CUTADAPT   = (
    shutil.which("cutadapt")
    or "/home/sgao30/micromamba/envs/droplet_single_cell/bin/cutadapt"
)

ME    = "AGATGTGTATAAGAGACAG"
ME_RC = "CTGTCTCTTATACACATCT"


def trim_sample(key: str):
    label   = SAMPLES[key]["label"]
    in_r1   = os.path.join(STEP1_DIR,  f"{key}_R1.fq.gz")
    in_r2   = os.path.join(FILTER_DIR, f"{key}_R2.fq.gz")
    out_r1  = os.path.join(OUT_DIR,    f"{key}_R1.fq.gz")
    out_r2  = os.path.join(OUT_DIR,    f"{key}_R2.fq.gz")
    log     = os.path.join(OUT_DIR,    f"{key}_cutadapt.log")
    json_out = os.path.join(OUT_DIR,   f"{key}_cutadapt.json")

    for path, name in ((in_r1, "step1 R1"), (in_r2, "filter_HD2 R2")):
        if not os.path.exists(path):
            print(f"[{label}] missing {name}: {path}, skipping", flush=True)
            return

    cmd = [
        CUTADAPT,
        "-j", "8",
        "-e", "0.1",
        "-O", "10",
        "-g", ME,
        "-a", ME_RC,
        "-A", ME_RC,
        "--minimum-length", "20",
        "--pair-filter=any",
        "--report=full",
        "--json", json_out,
        "-o", out_r1,
        "-p", out_r2,
        in_r1,
        in_r2,
    ]

    print(f"[{label}] running cutadapt ...", flush=True)
    result = subprocess.run(cmd, capture_output=True, text=True)

    with open(log, "w") as f:
        f.write(result.stdout)
        if result.stderr:
            f.write("\n--- stderr ---\n")
            f.write(result.stderr)

    if result.returncode != 0:
        print(f"[{label}] ERROR (see {os.path.basename(log)})", flush=True)
        print(result.stderr[:500], flush=True)
        return

    # print key summary lines
    for line in result.stdout.splitlines():
        if any(k in line for k in (
            "Total read pairs", "Pairs written", "Reads written",
            "with adapter", "too short",
        )):
            print(f"[{label}]  {line.strip()}", flush=True)

    print(f"[{label}] done", flush=True)


if __name__ == "__main__":
    os.makedirs(OUT_DIR, exist_ok=True)

    keys = list(SAMPLES.keys())
    print(f"Trimming {len(keys)} samples in parallel ...\n")

    with multiprocessing.Pool(processes=len(keys)) as pool:
        pool.map(trim_sample, keys)

    print("\nDone.")
