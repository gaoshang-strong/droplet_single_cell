"""
Filter raw R1/R2 reads using three Hamming-distance tiers.

All three filters share the same criteria:
  - W1_exact_hit == TRUE
  - gap_length == 20
  - hamming_distance <= max_hd  (0 / 2 / 3 depending on filter)

One pass through the metadata CSV + raw R1 + raw R2 per sample;
reads are written to all matching filter folders simultaneously.

Output structure:
  round3/
    filter_HD0/{key}_R1.fq.gz  {key}_R2.fq.gz
    filter_HD2/{key}_R1.fq.gz  {key}_R2.fq.gz
    filter_HD3/{key}_R1.fq.gz  {key}_R2.fq.gz
"""
import os
import gzip
import subprocess
import multiprocessing

from config import SAMPLES, OUTPUT_DIR, FILTERS, N_WORKERS


def process_sample(key: str) -> dict:
    info  = SAMPLES[key]
    label = info["label"]
    csv   = os.path.join(OUTPUT_DIR, f"{key}_metadata.csv.gz")
    r1    = info["R1"]
    r2    = info["R2"]

    print(f"[{label}] start", flush=True)

    # open output files for all filter tiers
    out_handles = []
    for folder, _ in FILTERS:
        out_r1 = os.path.join(OUTPUT_DIR, folder, f"{key}_R1.fq.gz")
        out_r2 = os.path.join(OUTPUT_DIR, folder, f"{key}_R2.fq.gz")
        out_handles.append((
            gzip.open(out_r1, "wb", compresslevel=1),
            gzip.open(out_r2, "wb", compresslevel=1),
        ))

    counts = [0] * len(FILTERS)   # passed reads per filter
    total  = 0

    p_csv = subprocess.Popen(["pigz", "-dc", "-p", "2", csv], stdout=subprocess.PIPE)
    p_r1  = subprocess.Popen(["pigz", "-dc", "-p", "2", r1],  stdout=subprocess.PIPE)
    p_r2  = subprocess.Popen(["pigz", "-dc", "-p", "2", r2],  stdout=subprocess.PIPE)

    p_csv.stdout.readline()   # skip CSV header

    while True:
        meta     = p_csv.stdout.readline()
        r1_block = [p_r1.stdout.readline() for _ in range(4)]
        r2_block = [p_r2.stdout.readline() for _ in range(4)]

        if not meta or not r1_block[0]:
            break
        total += 1

        parts   = meta.decode().rstrip("\n").split(",")
        hd      = int(parts[3])
        hit     = parts[4] == "TRUE"
        gap_raw = parts[7]

        if not hit or gap_raw == "NA" or int(gap_raw) != 20:
            continue

        for i, (_, max_hd) in enumerate(FILTERS):
            if hd <= max_hd:
                f1, f2 = out_handles[i]
                for line in r1_block: f1.write(line)
                for line in r2_block: f2.write(line)
                counts[i] += 1

    for f1, f2 in out_handles:
        f1.close()
        f2.close()

    p_csv.wait(); p_r1.wait(); p_r2.wait()

    result = {"label": label, "key": key, "total": total}
    for i, (folder, max_hd) in enumerate(FILTERS):
        result[folder] = counts[i]

    print(
        f"[{label}] done  total={total:,}  "
        + "  ".join(
            f"{folder}={counts[i]:,} ({counts[i]/total*100:.1f}%)"
            for i, (folder, _) in enumerate(FILTERS)
        ),
        flush=True,
    )
    return result


if __name__ == "__main__":
    # create output directories
    for folder, _ in FILTERS:
        os.makedirs(os.path.join(OUTPUT_DIR, folder), exist_ok=True)

    keys = list(SAMPLES.keys())
    print(f"Processing {len(keys)} samples with {N_WORKERS} workers ...\n")

    with multiprocessing.Pool(processes=N_WORKERS) as pool:
        results = pool.map(process_sample, keys)

    # write summary
    summary_path = os.path.join(OUTPUT_DIR, "filter_summary.tsv")
    header = "sample\ttotal\t" + "\t".join(
        f"{f}_passed\t{f}_pct" for f, _ in FILTERS
    )
    with open(summary_path, "w") as f:
        f.write(header + "\n")
        for r in results:
            t    = r["total"]
            cols = [r["label"], str(t)]
            for folder, _ in FILTERS:
                n = r[folder]
                cols += [str(n), f"{n/t*100:.4f}"]
            f.write("\t".join(cols) + "\n")

    print(f"\nAll done. Summary -> {summary_path}")
