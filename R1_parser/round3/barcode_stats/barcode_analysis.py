"""
Extract BC1, BC2, and combined cell barcode (BC1+BC2) from filter_HD2 R1
read headers, count frequencies, and generate a QC report.

Input : filter_HD2 *_R1.fq.gz  (headers carry bc1:Z and bc2:Z tags)
Output per sample:
  {key}_BC1_counts.tsv.gz       BC1 sequence → read count
  {key}_BC2_counts.tsv.gz       BC2 sequence → read count
  {key}_CB_counts.tsv.gz        CB (BC1+BC2) → read count
  {key}_barcode_stats.json      summary stats (cached)
Across samples:
  barcode_summary.tsv           per-sample summary table
  knee_plot.png                 read-per-CB rank plot (log-log)
  barcode_summary.md            markdown report
"""
import os
import re
import gzip
import json
import subprocess
import multiprocessing
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

sys_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys; sys.path.insert(0, sys_path)
from config import SAMPLES, OUTPUT_DIR

FILTER_DIR  = os.path.join(OUTPUT_DIR, "filter_HD2")
STATS_DIR   = os.path.join(OUTPUT_DIR, "barcode_stats")
SAMPLE_ORDER = ["1PB", "2PB", "3PB", "3PY", "6PY", "9PY"]

BC1_RE = re.compile(r'bc1:Z:([ACGTN]+)')
BC2_RE = re.compile(r'bc2:Z:([ACGTN]+)')

C_PB = "#1565C0"
C_PY = "#E65100"

def sample_color(label):
    return C_PB if label.endswith("PB") else C_PY


# ── per-sample counting ───────────────────────────────────────────────────────

def count_barcodes(key: str) -> dict:
    label     = SAMPLES[key]["label"]
    cache     = os.path.join(STATS_DIR, f"{key}_barcode_stats.json")

    if os.path.exists(cache):
        print(f"[{label}] loaded cached stats", flush=True)
        with open(cache) as f:
            return json.load(f)

    r1_path   = os.path.join(FILTER_DIR, f"{key}_R1.fq.gz")
    print(f"[{label}] counting ...", flush=True)

    bc1_cnt = Counter()
    bc2_cnt = Counter()
    cb_cnt  = Counter()   # combined cell barcode BC1+BC2
    total   = 0

    proc = subprocess.Popen(["pigz", "-dc", "-p", "2", r1_path], stdout=subprocess.PIPE)
    for raw in proc.stdout:
        if not raw.startswith(b"@"):
            continue
        total += 1
        header = raw.decode()
        m1 = BC1_RE.search(header)
        m2 = BC2_RE.search(header)
        if m1 and m2:
            bc1 = m1.group(1)
            bc2 = m2.group(1)
            bc1_cnt[bc1] += 1
            bc2_cnt[bc2] += 1
            cb_cnt[bc1 + bc2] += 1
    proc.wait()

    # write count tables
    def write_counts(counter, path):
        with gzip.open(path, "wt") as f:
            f.write("sequence\tcount\n")
            for seq, cnt in counter.most_common():
                f.write(f"{seq}\t{cnt}\n")

    write_counts(bc1_cnt, os.path.join(STATS_DIR, f"{key}_BC1_counts.tsv.gz"))
    write_counts(bc2_cnt, os.path.join(STATS_DIR, f"{key}_BC2_counts.tsv.gz"))
    write_counts(cb_cnt,  os.path.join(STATS_DIR, f"{key}_CB_counts.tsv.gz"))

    # summary stats
    cb_counts = sorted(cb_cnt.values(), reverse=True)
    stats = {
        "label":        label,
        "key":          key,
        "total_reads":  total,
        "unique_BC1":   len(bc1_cnt),
        "unique_BC2":   len(bc2_cnt),
        "unique_CB":    len(cb_cnt),
        "cb_counts_top1000": cb_counts[:1000],   # store top 1000 for knee plot
        "cb_median_reads":   int(np.median(cb_counts)) if cb_counts else 0,
        "cb_mean_reads":     float(np.mean(cb_counts))  if cb_counts else 0,
        "reads_in_top1k_cb": int(sum(cb_counts[:1000])),
    }

    with open(cache, "w") as f:
        json.dump(stats, f, indent=2)

    print(
        f"[{label}] done  reads={total:,}  unique_CB={len(cb_cnt):,}"
        f"  median_reads/CB={stats['cb_median_reads']}",
        flush=True,
    )
    return stats


# ── plots ─────────────────────────────────────────────────────────────────────

def save(fig, name):
    path = os.path.join(STATS_DIR, name)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {name}")


def plot_knee(all_stats):
    """Log-log knee plot: reads per CB vs CB rank, all samples overlaid."""
    fig, ax = plt.subplots(figsize=(9, 6))

    for s in all_stats:
        counts = s["cb_counts_top1000"]
        if not counts:
            continue
        ranks = np.arange(1, len(counts) + 1)
        ax.plot(ranks, counts, label=s["label"],
                color=sample_color(s["label"]),
                linewidth=1.8,
                linestyle="-" if s["label"].endswith("PB") else "--")

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Cell barcode rank", fontsize=11)
    ax.set_ylabel("Reads per cell barcode", fontsize=11)
    ax.set_title("Knee plot — top 1,000 cell barcodes (BC1+BC2)",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(True, which="both", alpha=0.3)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    save(fig, "knee_plot.png")


def plot_unique_cb(all_stats):
    """Bar: number of unique cell barcodes per sample."""
    labels = [s["label"] for s in all_stats]
    n      = len(labels)
    x      = np.arange(n)
    vals   = [s["unique_CB"] / 1e6 for s in all_stats]
    colors = [sample_color(l) for l in labels]

    fig, ax = plt.subplots(figsize=(8, 4.5))
    bars = ax.bar(x, vals, color=colors, width=0.55)
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                f"{v:.2f}M", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylabel("Unique cell barcodes (millions)", fontsize=10)
    ax.set_title("Unique BC1+BC2 combinations per sample", fontsize=12, fontweight="bold")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    save(fig, "unique_CB.png")


def plot_reads_per_cb(all_stats):
    """Bar: median reads per cell barcode per sample."""
    labels = [s["label"] for s in all_stats]
    n      = len(labels)
    x      = np.arange(n)
    vals   = [s["cb_median_reads"] for s in all_stats]
    colors = [sample_color(l) for l in labels]

    fig, ax = plt.subplots(figsize=(8, 4.5))
    bars = ax.bar(x, vals, color=colors, width=0.55)
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.2,
                str(v), ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylabel("Median reads per cell barcode", fontsize=10)
    ax.set_title("Median reads per cell barcode (BC1+BC2)", fontsize=12, fontweight="bold")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    save(fig, "median_reads_per_CB.png")


# ── markdown report ───────────────────────────────────────────────────────────

def write_report(all_stats):
    by_label = {s["label"]: s for s in all_stats}

    md = []
    md += [
        "# Barcode Statistics Report",
        "",
        "**Input:** filter_HD2 R1 reads (HD≤2 + W1 exact hit + gap=20)  ",
        "**Cell barcode:** BC1 (10 bp) + BC2 (10 bp) = 20 bp combined (CB)",
        "",
        "---",
        "",
        "## 1. Summary Table",
        "",
        "| Sample | Total reads | Unique BC1 | Unique BC2 | Unique CB | Median reads/CB | Mean reads/CB |",
        "|--------|------------|------------|------------|-----------|-----------------|---------------|",
    ]

    for label in SAMPLE_ORDER:
        if label not in by_label:
            continue
        s = by_label[label]
        md.append(
            f"| **{label}** | {s['total_reads']:,} | {s['unique_BC1']:,} |"
            f" {s['unique_BC2']:,} | {s['unique_CB']:,} |"
            f" {s['cb_median_reads']:,} | {s['cb_mean_reads']:.1f} |"
        )

    md += [
        "",
        "---",
        "",
        "## 2. Unique Cell Barcodes (BC1+BC2)",
        "",
        "![Unique CB](unique_CB.png)",
        "",
        "---",
        "",
        "## 3. Median Reads per Cell Barcode",
        "",
        "![Median reads per CB](median_reads_per_CB.png)",
        "",
        "---",
        "",
        "## 4. Knee Plot",
        "",
        "Read count vs rank for the top 1,000 cell barcodes (log-log scale).  ",
        "A steep drop-off indicates a clear separation between real cells and background.",
        "",
        "![Knee plot](knee_plot.png)",
        "",
        "---",
        "",
        "## 5. Observations",
        "",
    ]

    # auto observations
    pb_cbs = [by_label[l]["unique_CB"] for l in ["1PB","2PB","3PB"] if l in by_label]
    py_cbs = [by_label[l]["unique_CB"] for l in ["3PY","6PY","9PY"] if l in by_label]
    if pb_cbs and py_cbs:
        md.append(
            f"- PB samples have {min(pb_cbs):,}–{max(pb_cbs):,} unique cell barcodes; "
            f"PY samples have {min(py_cbs):,}–{max(py_cbs):,}."
        )

    pb_med = [by_label[l]["cb_median_reads"] for l in ["1PB","2PB","3PB"] if l in by_label]
    py_med = [by_label[l]["cb_median_reads"] for l in ["3PY","6PY","9PY"] if l in by_label]
    if pb_med:
        md.append(
            f"- Median reads per CB: PB {min(pb_med)}–{max(pb_med)}, "
            f"PY {min(py_med)}–{max(py_med)}."
        )

    md += [
        "- The knee plot shape indicates whether a clear cell/background boundary exists. "
        "A sharp elbow suggests good library complexity; a flat curve suggests many low-count barcodes.",
        "- High unique CB counts with low median reads per CB typically indicate a large fraction "
        "of background/empty droplets — downstream cell calling (e.g., ArchR, snapATAC2) "
        "will apply a count threshold to separate real cells.",
    ]

    path = os.path.join(STATS_DIR, "barcode_summary.md")
    with open(path, "w") as f:
        f.write("\n".join(md))
    print(f"Saved barcode_summary.md")


# ── main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    os.makedirs(STATS_DIR, exist_ok=True)

    keys = list(SAMPLES.keys())
    print(f"Processing {len(keys)} samples with {len(keys)} workers ...\n")

    with multiprocessing.Pool(processes=len(keys)) as pool:
        all_stats = pool.map(count_barcodes, keys)

    label_map = {s["label"]: s for s in all_stats}
    all_stats  = [label_map[l] for l in SAMPLE_ORDER if l in label_map]

    print("\nGenerating plots ...")
    plot_knee(all_stats)
    plot_unique_cb(all_stats)
    plot_reads_per_cb(all_stats)
    write_report(all_stats)
    print("Done.")
