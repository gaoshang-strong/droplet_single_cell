"""
Stream each sample's raw R1 fastq, apply capture and W1 scanners,
and write one metadata CSV per sample.

Output columns (1-based closed coords):
  read_name, capture_seq_start, capture_seq_end, hamming_distance,
  W1_exact_hit, W1_start, W1_end, gap_length
"""
import os
import gzip
import subprocess
import multiprocessing

from config import SAMPLES, OUTPUT_DIR, N_WORKERS
from scanner import scan_capture, scan_W1

HEADER = (
    "read_name,capture_seq_start,capture_seq_end,hamming_distance,"
    "W1_exact_hit,W1_start,W1_end,gap_length\n"
)


def _gap(capture_start, w1_end) -> str:
    """gap_length = bases between W1 end and capture start (exclusive both ends)."""
    if capture_start is None or w1_end is None:
        return "NA"
    return str(capture_start - w1_end - 1)


def process_sample(key: str) -> dict:
    info  = SAMPLES[key]
    label = info["label"]
    r1    = info["R1"]
    out   = os.path.join(OUTPUT_DIR, f"{key}_metadata.csv.gz")

    print(f"[{label}] start", flush=True)

    total = 0
    proc  = subprocess.Popen(["pigz", "-dc", "-p", "2", r1], stdout=subprocess.PIPE)

    with gzip.open(out, "wt", compresslevel=1) as f:
        f.write(HEADER)
        while True:
            header = proc.stdout.readline()
            seq    = proc.stdout.readline()
            proc.stdout.readline()   # +
            proc.stdout.readline()   # qual
            if not header:
                break
            total += 1

            read_name = header.decode().split()[0].lstrip("@")
            seq_str   = seq.decode().rstrip()

            cs, ce, hd   = scan_capture(seq_str)
            hit, ws, we  = scan_W1(seq_str)

            f.write(
                f"{read_name},"
                f"{cs if cs is not None else 'NA'},"
                f"{ce if ce is not None else 'NA'},"
                f"{hd if hd is not None else 'NA'},"
                f"{'TRUE' if hit else 'FALSE'},"
                f"{ws if ws is not None else 'NA'},"
                f"{we if we is not None else 'NA'},"
                f"{_gap(cs, we)}\n"
            )

    proc.wait()
    print(f"[{label}] done  total={total:,}", flush=True)
    return {"label": label, "key": key, "total": total}


if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    keys = list(SAMPLES.keys())
    print(f"Processing {len(keys)} samples with {N_WORKERS} workers ...\n")

    with multiprocessing.Pool(processes=N_WORKERS) as pool:
        results = pool.map(process_sample, keys)

    summary_path = os.path.join(OUTPUT_DIR, "run_summary.tsv")
    with open(summary_path, "w") as f:
        f.write("sample\ttotal_reads\n")
        for r in results:
            f.write(f"{r['label']}\t{r['total']}\n")

    print(f"\nAll done. Summary -> {summary_path}")
