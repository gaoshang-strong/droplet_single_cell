#!/usr/bin/env python3
"""
Analyze step7 barcode/UMI TSV files.

Per-sample streaming statistics:
  1. cs_pos - w1_pos distribution  (structural integrity)
  2. gap_len distribution
  3. gap_seq Hamming distance to TAAGGCGA  (gap_len==8 only)
  4. TAAGGCGA exact-match rate
  5. bc1_len distribution
  6. mt (exact / hamming) breakdown

Outputs:
  - summary TSVs per sample + one cross-sample summary
  - PNG plots
  - Markdown report
"""
import os
import glob
import re
import subprocess
import multiprocessing
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

STEP7_DIR  = "/ShangGaoAIProjects/ZhangJW/R1_parser/round1/step7_extract_barcode_UMI"
OUTPUT_DIR = STEP7_DIR

COMMON_FIXED_REF = "TAAGGCGA"

SAMPLE_LABELS = {
    "260430R-S-XY-3PY": "3PY",
    "260430R-S-XY-1PB": "1PB",
    "260430R-S-XY-2PB": "2PB",
    "260430R-S-XY-6PY": "6PY",
    "260430R-S-XY-9PY": "9PY",
    "260430R-S-XY-3PB": "3PB",
}
SAMPLE_ORDER = ["1PB", "2PB", "3PB", "3PY", "6PY", "9PY"]

def extract_key(path):
    for key in SAMPLE_LABELS:
        if key in path:
            return key
    return None

def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


# ── streaming analysis ────────────────────────────────────────────────────────

def analyze_sample(tsv_gz_path):
    key   = extract_key(tsv_gz_path)
    label = SAMPLE_LABELS[key]
    print(f"[{label}] start", flush=True)

    total          = 0
    mt_count       = {}          # exact / hamming
    bc1_len_dist   = {}          # 0..10
    gap_len_dist   = {}          # 0..30+
    delta_dist     = {}          # cs_pos - w1_pos
    hamming_dist   = {}          # 0..8  (gap_len==8 only)
    taaggcga_exact = 0           # gap_len==8 and gap_seq==REF
    gap8_total     = 0           # gap_len==8 reads

    proc = subprocess.Popen(["pigz", "-dc", "-p", "2", tsv_gz_path],
                            stdout=subprocess.PIPE)
    proc.stdout.readline()   # skip header

    for raw in proc.stdout:
        parts = raw.decode().rstrip("\n").split("\t")
        if len(parts) < 10:
            continue
        total += 1

        # columns: read_id w1_pos cs_pos mt bc1 bc1_len bc2 umi gap_len gap_seq
        w1_pos  = int(parts[1])
        cs_pos  = int(parts[2])
        mt      = parts[3]
        bc1_len = int(parts[5])
        gap_len = int(parts[8])
        gap_seq = parts[9]

        mt_count[mt]               = mt_count.get(mt, 0) + 1
        bc1_len_dist[bc1_len]      = bc1_len_dist.get(bc1_len, 0) + 1
        gap_len_dist[gap_len]      = gap_len_dist.get(gap_len, 0) + 1
        delta = cs_pos - w1_pos
        delta_dist[delta]          = delta_dist.get(delta, 0) + 1

        if gap_len == 8:
            gap8_total += 1
            if gap_seq == COMMON_FIXED_REF:
                taaggcga_exact += 1
            hd = hamming(gap_seq, COMMON_FIXED_REF)
            hamming_dist[hd] = hamming_dist.get(hd, 0) + 1

    proc.wait()

    stats = {
        "label":          label,
        "key":            key,
        "total":          total,
        "mt_count":       mt_count,
        "bc1_len_dist":   bc1_len_dist,
        "gap_len_dist":   gap_len_dist,
        "delta_dist":     {str(k): v for k, v in delta_dist.items()},
        "hamming_dist":   hamming_dist,
        "taaggcga_exact": taaggcga_exact,
        "gap8_total":     gap8_total,
    }

    json_path = os.path.join(OUTPUT_DIR, f"{key}_stats.json")
    with open(json_path, "w") as f:
        json.dump(stats, f, indent=2)

    print(
        f"[{label}] done  total={total:,}  gap8={gap8_total:,}"
        f"  TAAGGCGA_exact={taaggcga_exact:,} ({taaggcga_exact/total*100:.1f}%)",
        flush=True,
    )
    return stats


# ── plotting helpers ──────────────────────────────────────────────────────────

C_PB     = "#1565C0"
C_PY     = "#E65100"
C_GRAY   = "#90A4AE"
C_GREEN  = "#2E7D32"
C_RED    = "#C62828"

def sample_color(label):
    return C_PB if label.endswith("PB") else C_PY

def save(fig, name):
    p = os.path.join(OUTPUT_DIR, name)
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {name}")


def plot_taaggcga_rate(all_stats):
    """Bar chart: TAAGGCGA exact-match rate per sample (among all reads)."""
    labels = [s["label"] for s in all_stats]
    rates  = [s["taaggcga_exact"] / s["total"] * 100 for s in all_stats]
    colors = [sample_color(l) for l in labels]

    fig, ax = plt.subplots(figsize=(8, 4.5))
    bars = ax.bar(labels, rates, color=colors, width=0.55)
    for bar, v in zip(bars, rates):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.4,
                f"{v:.1f}%", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.set_ylim(0, 80)
    ax.set_ylabel("% of all reads", fontsize=10)
    ax.set_title("TAAGGCGA exact match rate\n(gap_len=8, among all reads)", fontsize=12, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=C_PB, label="PB batch"),
                        Patch(color=C_PY, label="PY batch")], fontsize=9)
    plt.tight_layout()
    save(fig, "taaggcga_rate.png")


def plot_gap_len(all_stats):
    """Stacked bar: gap_len=8 vs other, per sample."""
    labels  = [s["label"] for s in all_stats]
    n       = len(labels)
    x       = np.arange(n)
    gap8    = np.array([s["gap_len_dist"].get(8, 0) / s["total"] * 100 for s in all_stats])
    other   = 100 - gap8

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(x, gap8,  color=C_GREEN,  label="gap_len = 8 (canonical)")
    ax.bar(x, other, bottom=gap8, color=C_RED, label="gap_len ≠ 8")
    for i in range(n):
        ax.text(i, gap8[i] / 2,           f"{gap8[i]:.1f}%",
                ha="center", va="center", fontsize=9, color="white", fontweight="bold")
        ax.text(i, gap8[i] + other[i] / 2, f"{other[i]:.1f}%",
                ha="center", va="center", fontsize=9, color="white", fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylim(0, 108)
    ax.set_ylabel("% of reads", fontsize=10)
    ax.set_title("Gap length (between UMI and capture start)", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    save(fig, "gap_len_breakdown.png")


def plot_hamming(all_stats):
    """Grouped bar: Hamming distance to TAAGGCGA (gap_len=8), per sample."""
    max_hd  = 8
    x       = np.arange(max_hd + 1)
    width   = 0.13
    offsets = np.linspace(-(len(all_stats)-1)/2, (len(all_stats)-1)/2, len(all_stats)) * width

    fig, ax = plt.subplots(figsize=(11, 5))
    for i, s in enumerate(all_stats):
        g8 = s["gap8_total"] or 1
        vals = [s["hamming_dist"].get(hd, 0) / g8 * 100 for hd in range(max_hd + 1)]
        ax.bar(x + offsets[i], vals, width=width,
               label=s["label"], color=sample_color(s["label"]),
               alpha=0.85 if i % 3 == 0 else (0.65 if i % 3 == 1 else 0.45))

    ax.set_xticks(x)
    ax.set_xticklabels([str(h) for h in range(max_hd + 1)], fontsize=11)
    ax.set_xlabel("Hamming distance to TAAGGCGA", fontsize=10)
    ax.set_ylabel("% of gap_len=8 reads", fontsize=10)
    ax.set_title("Hamming distance distribution to TAAGGCGA\n(gap_len=8 reads)", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9, ncol=2)
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    save(fig, "hamming_dist.png")


def plot_delta_dist(all_stats):
    """Bar: % of reads with cs_pos - w1_pos != 25 per sample."""
    labels  = [s["label"] for s in all_stats]
    n       = len(labels)
    x       = np.arange(n)
    ok      = np.array([s["delta_dist"].get("25", 0) / s["total"] * 100 for s in all_stats])
    bad     = 100 - ok

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(x, ok,  color=C_GREEN, label="cs_pos − w1_pos = 25 (canonical)")
    ax.bar(x, bad, bottom=ok, color=C_RED, label="≠ 25 (structural anomaly)")
    for i in range(n):
        ax.text(i, ok[i] / 2,         f"{ok[i]:.1f}%",
                ha="center", va="center", fontsize=9, color="white", fontweight="bold")
        ax.text(i, ok[i] + bad[i] / 2, f"{bad[i]:.1f}%",
                ha="center", va="center", fontsize=9, color="white", fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylim(0, 108)
    ax.set_ylabel("% of reads", fontsize=10)
    ax.set_title("Structural integrity: cs_pos − w1_pos\n(canonical = 25)", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    save(fig, "structural_integrity.png")


def plot_bc1_len(all_stats):
    """Bar: bc1_len < 10 rate per sample."""
    labels  = [s["label"] for s in all_stats]
    n       = len(labels)
    x       = np.arange(n)
    full    = np.array([s["bc1_len_dist"].get(10, 0) / s["total"] * 100 for s in all_stats])
    trunc   = 100 - full

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(x, full,  color=C_GREEN, label="BC1 = 10 bp (full)")
    ax.bar(x, trunc, bottom=full, color=C_RED, label="BC1 < 10 bp (truncated)")
    for i in range(n):
        ax.text(i, full[i] / 2,           f"{full[i]:.1f}%",
                ha="center", va="center", fontsize=9, color="white", fontweight="bold")
        if trunc[i] > 0.05:
            ax.text(i, full[i] + trunc[i] / 2, f"{trunc[i]:.2f}%",
                    ha="center", va="center", fontsize=9, color="white", fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylim(0, 108)
    ax.set_ylabel("% of reads", fontsize=10)
    ax.set_title("BC1 length (expected 10 bp)", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    save(fig, "bc1_len.png")


# ── markdown report ───────────────────────────────────────────────────────────

def write_report(all_stats):
    by_label = {s["label"]: s for s in all_stats}

    md = []
    md.append("# Step 7 Analysis Report — Barcode / UMI Extraction QC")
    md.append("")
    md.append("**Input:** step6 W1-filtered reads (`*_W1_R1.fq.gz`)")
    md.append("")
    md.append("**Reference common_fixed:** `TAAGGCGA` (8 bp, known for PB batch)")
    md.append("")
    md.append("---")
    md.append("")
    md.append("## 1. Structural Integrity — cs_pos − w1_pos")
    md.append("")
    md.append("Expected value: **25** (canonical read structure).  ")
    md.append("Deviation indicates indels or a mis-anchored capture position.")
    md.append("")
    md.append("![Structural integrity](structural_integrity.png)")
    md.append("")

    # delta table
    md.append("| Sample | cs−w1=25 | cs−w1=24 | cs−w1=26 | cs−w1≥27 or ≤23 |")
    md.append("|--------|----------|----------|----------|-----------------|")
    for label in SAMPLE_ORDER:
        s  = by_label[label]
        t  = s["total"]
        dd = {int(k): v for k, v in s["delta_dist"].items()}
        c25 = dd.get(25, 0)
        c24 = dd.get(24, 0)
        c26 = dd.get(26, 0)
        rest = t - c25 - c24 - c26
        md.append(
            f"| **{label}** | {c25/t*100:.1f}% ({c25:,}) "
            f"| {c24/t*100:.1f}% ({c24:,}) "
            f"| {c26/t*100:.1f}% ({c26:,}) "
            f"| {rest/t*100:.1f}% ({rest:,}) |"
        )
    md.append("")
    md.append("---")
    md.append("")
    md.append("## 2. Gap Length Distribution")
    md.append("")
    md.append("Gap = sequence between UMI end and capture start.  ")
    md.append("Canonical length: **8 bp** (common_fixed region).  ")
    md.append("gap_len ≠ 8 implies an indel in common_fixed or a misidentified capture position.")
    md.append("")
    md.append("![Gap length](gap_len_breakdown.png)")
    md.append("")

    # gap_len table
    md.append("| Sample | gap=8 | gap=7 | gap=9 | gap≤6 or ≥10 |")
    md.append("|--------|-------|-------|-------|--------------|")
    for label in SAMPLE_ORDER:
        s  = by_label[label]
        t  = s["total"]
        gd = s["gap_len_dist"]
        g8 = gd.get(8, 0); g7 = gd.get(7, 0); g9 = gd.get(9, 0)
        rest = t - g8 - g7 - g9
        md.append(
            f"| **{label}** | {g8/t*100:.1f}% ({g8:,}) "
            f"| {g7/t*100:.1f}% ({g7:,}) "
            f"| {g9/t*100:.1f}% ({g9:,}) "
            f"| {rest/t*100:.1f}% ({rest:,}) |"
        )
    md.append("")
    md.append("---")
    md.append("")
    md.append("## 3. TAAGGCGA Match Rate")
    md.append("")
    md.append("Exact match of gap_seq to `TAAGGCGA` among all reads.")
    md.append("")
    md.append("![TAAGGCGA rate](taaggcga_rate.png)")
    md.append("")
    md.append("| Sample | Batch | Total | TAAGGCGA exact | Rate |")
    md.append("|--------|-------|-------|----------------|------|")
    for label in SAMPLE_ORDER:
        s     = by_label[label]
        batch = "PB" if label.endswith("PB") else "PY"
        md.append(
            f"| **{label}** | {batch} | {s['total']:,} "
            f"| {s['taaggcga_exact']:,} "
            f"| **{s['taaggcga_exact']/s['total']*100:.1f}%** |"
        )
    md.append("")
    md.append("---")
    md.append("")
    md.append("## 4. Hamming Distance to TAAGGCGA")
    md.append("")
    md.append("Computed for gap_len=8 reads only.")
    md.append("")
    md.append("![Hamming dist](hamming_dist.png)")
    md.append("")
    md.append("| Sample | HD=0 | HD=1 | HD=2 | HD=3 | HD=4 | HD≥5 |")
    md.append("|--------|------|------|------|------|------|------|")
    for label in SAMPLE_ORDER:
        s  = by_label[label]
        g8 = s["gap8_total"] or 1
        hd = s["hamming_dist"]
        h0 = hd.get(0, 0); h1 = hd.get(1, 0); h2 = hd.get(2, 0)
        h3 = hd.get(3, 0); h4 = hd.get(4, 0)
        h5p = g8 - h0 - h1 - h2 - h3 - h4
        md.append(
            f"| **{label}** | {h0/g8*100:.1f}% | {h1/g8*100:.1f}% "
            f"| {h2/g8*100:.1f}% | {h3/g8*100:.1f}% "
            f"| {h4/g8*100:.1f}% | {h5p/g8*100:.1f}% |"
        )
    md.append("")
    md.append("---")
    md.append("")
    md.append("## 5. BC1 Length Distribution")
    md.append("")
    md.append("BC1 is extracted as the 10 nt before W1.  ")
    md.append("BC1 < 10 bp occurs when W1 is within 10 bp of the read start.")
    md.append("")
    md.append("![BC1 length](bc1_len.png)")
    md.append("")
    md.append("| Sample | BC1=10 bp | BC1<10 bp | Min BC1 len |")
    md.append("|--------|-----------|-----------|-------------|")
    for label in SAMPLE_ORDER:
        s   = by_label[label]
        t   = s["total"]
        bd  = s["bc1_len_dist"]
        full = bd.get(10, 0)
        trunc = t - full
        min_len = min(int(k) for k in bd if bd[k] > 0)
        md.append(
            f"| **{label}** | {full/t*100:.2f}% ({full:,}) "
            f"| {trunc/t*100:.2f}% ({trunc:,}) "
            f"| {min_len} |"
        )
    md.append("")
    md.append("---")
    md.append("")
    md.append("## 6. Key Findings")
    md.append("")
    md.append("### 6.1 PY batch has low TAAGGCGA rate and highly diverse gap sequences")
    md.append("")

    pb_rates = [by_label[l]["taaggcga_exact"]/by_label[l]["total"]*100
                for l in ["1PB","2PB","3PB"]]
    py_rates = [by_label[l]["taaggcga_exact"]/by_label[l]["total"]*100
                for l in ["3PY","6PY","9PY"]]
    md.append(
        f"- PB samples: {min(pb_rates):.0f}–{max(pb_rates):.0f}% TAAGGCGA exact match.  "
    )
    md.append(
        f"- PY samples: {min(py_rates):.0f}–{max(py_rates):.0f}% TAAGGCGA exact match — "
        f"the remaining ~{100-max(py_rates):.0f}% carry highly diverse 8-mer sequences.  "
    )
    md.append("- This confirms that **PY's common_fixed is NOT a fixed sequence**: either the region encodes a variable barcode, or a different library preparation was used.")
    md.append("")
    md.append("### 6.2 ~4–9% of reads have gap_len ≠ 8")
    md.append("")
    md.append("- Most anomalies are gap_len=7 (−1 bp): likely a 1-bp deletion in the common_fixed region or a 1-bp shift in the Hamming-rescued capture position.")
    md.append("- A small fraction of PY reads have very large gaps (cs_pos−w1_pos = 44, 56, 63…), indicating that the Hamming rescue identified a false-positive capture position far downstream.")
    md.append("")
    md.append("### 6.3 PY Hamming distance distribution is flat")
    md.append("")
    md.append("- For PB, Hamming distance to TAAGGCGA is bimodal: HD=0 dominates (~60%), then a near-flat tail.")
    md.append("- For PY, HD=0 is ~42% and the distribution is nearly uniform across HD 1–8, consistent with random sequence rather than a single alternative fixed sequence.")
    md.append("")
    md.append("### 6.4 BC1 truncation is rare")
    md.append("")
    md.append("- < 1% of reads have BC1 < 10 bp (W1 appears within 10 bp of the read start).  ")
    md.append("- These reads can be excluded from downstream cell barcode analysis without significant loss.")

    md_path = os.path.join(OUTPUT_DIR, "step7_analysis_report.md")
    with open(md_path, "w") as f:
        f.write("\n".join(md))
    print(f"Saved {md_path}")


# ── main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    tsv_files = sorted(glob.glob(os.path.join(STEP7_DIR, "*_bc_umi.tsv.gz")))
    tsv_files = [f for f in tsv_files if extract_key(f) is not None]

    # check if cached JSON stats exist for all samples
    all_stats = []
    to_process = []
    for f in tsv_files:
        key  = extract_key(f)
        json_path = os.path.join(OUTPUT_DIR, f"{key}_stats.json")
        if os.path.exists(json_path):
            with open(json_path) as jf:
                all_stats.append(json.load(jf))
            print(f"[{SAMPLE_LABELS[key]}] loaded cached stats")
        else:
            to_process.append(f)

    if to_process:
        n_workers = min(len(to_process), multiprocessing.cpu_count())
        print(f"\nProcessing {len(to_process)} samples with {n_workers} workers ...")
        with multiprocessing.Pool(processes=n_workers) as pool:
            new_stats = pool.map(analyze_sample, to_process)
        all_stats.extend(new_stats)

    # sort by SAMPLE_ORDER
    label_to_stats = {s["label"]: s for s in all_stats}
    all_stats = [label_to_stats[l] for l in SAMPLE_ORDER if l in label_to_stats]

    print("\nGenerating plots ...")
    plot_taaggcga_rate(all_stats)
    plot_gap_len(all_stats)
    plot_hamming(all_stats)
    plot_delta_dist(all_stats)
    plot_bc1_len(all_stats)

    write_report(all_stats)
    print("\nAll done.")
