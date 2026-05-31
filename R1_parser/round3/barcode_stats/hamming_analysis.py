"""
Hamming distance analysis for control experiment.

All barcodes are expected to be the reference sequence:
  BC1 = GATCGATCGA   (10 bp)
  BC2 = CTAGCTAGCT   (10 bp)
  CB  = GATCGATCGACTAGCTAGCT  (20 bp, BC1 + BC2 concatenated)

Reads per-sample count TSVs produced by barcode_analysis.py.
Outputs:
  hd_dist_BC1.png / hd_dist_BC2.png / hd_dist_CB.png
  hd0_summary.png
  hamming_report.md
"""
import os
import gzip
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys; sys.path.insert(0, sys_path)
from config import SAMPLES, OUTPUT_DIR

STATS_DIR    = os.path.join(OUTPUT_DIR, "barcode_stats")
SAMPLE_ORDER = ["1PB", "2PB", "3PB", "3PY", "6PY", "9PY"]

BC1_REF = "GATCGATCGA"
BC2_REF = "CTAGCTAGCT"
CB_REF  = BC1_REF + BC2_REF   # GATCGATCGACTAGCTAGCT

C_PB = "#1565C0"
C_PY = "#E65100"


def sample_color(label):
    return C_PB if label.endswith("PB") else C_PY


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def load_hd_dist(key, btype, ref):
    """Return (hd→count dict, total_reads) for one sample / barcode type."""
    path = os.path.join(STATS_DIR, f"{key}_{btype}_counts.tsv.gz")
    hd_dist = {}
    total = 0
    if not os.path.exists(path):
        print(f"  [warn] missing {os.path.basename(path)}", flush=True)
        return hd_dist, total
    with gzip.open(path, "rt") as f:
        next(f)  # skip header
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) < 2:
                continue
            seq, cnt = parts[0], int(parts[1])
            if len(seq) != len(ref):
                continue
            hd = hamming(seq, ref)
            hd_dist[hd] = hd_dist.get(hd, 0) + cnt
            total += cnt
    return hd_dist, total


# ── per-barcode-type figure: 2×3 grid, one subplot per sample ────────────────

def plot_hd_per_sample(btype, ref, key_to_data):
    label_key    = {SAMPLES[k]["label"]: k for k in SAMPLES}
    ordered_keys = [label_key[l] for l in SAMPLE_ORDER if label_key.get(l) in key_to_data]

    n     = len(ordered_keys)
    ncols = 3
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(13, 4.5 * nrows))
    axes = np.array(axes).flatten()

    for i, k in enumerate(ordered_keys):
        label = SAMPLES[k]["label"]
        hd_dist, total = key_to_data[k]
        ax = axes[i]

        if not hd_dist or total == 0:
            ax.set_title(label)
            ax.text(0.5, 0.5, "no data", ha="center", va="center",
                    transform=ax.transAxes)
            continue

        max_hd   = max(hd_dist.keys())
        hd_range = np.arange(0, max_hd + 1)
        counts   = np.array([hd_dist.get(hd, 0) for hd in hd_range])
        fracs    = counts / total * 100

        ax.bar(hd_range, fracs, color=sample_color(label), width=0.7, alpha=0.85)
        hd0_pct = fracs[0] if len(fracs) > 0 else 0
        ax.set_title(f"{label}   HD=0: {hd0_pct:.1f}%", fontsize=11)
        ax.set_xlabel("Hamming distance to reference")
        ax.set_ylabel("% of reads")
        ax.set_xticks(hd_range)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(f"{btype}  Hamming distance to reference  [ {ref} ]",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    out = os.path.join(STATS_DIR, f"hd_dist_{btype}.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved hd_dist_{btype}.png")


# ── summary figure: % HD=0 per sample for BC1 / BC2 / CB ────────────────────

def plot_hd0_summary(all_data):
    label_key = {SAMPLES[k]["label"]: k for k in SAMPLES}
    keys      = [label_key[l] for l in SAMPLE_ORDER if label_key.get(l) in all_data["BC1"]]
    labels    = [SAMPLES[k]["label"] for k in keys]
    x         = np.arange(len(keys))
    width     = 0.25
    colors    = {"BC1": "#2196F3", "BC2": "#FF9800", "CB": "#4CAF50"}

    fig, ax = plt.subplots(figsize=(11, 5))

    for i, btype in enumerate(("BC1", "BC2", "CB")):
        vals = []
        for k in keys:
            hd_dist, total = all_data[btype].get(k, ({}, 0))
            vals.append(hd_dist.get(0, 0) / total * 100 if total else 0)
        offset = (i - 1) * width
        bars = ax.bar(x + offset, vals, width, label=btype,
                      color=colors[btype], alpha=0.85)
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.4,
                    f"{v:.1f}%", ha="center", va="bottom", fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylabel("% reads at HD = 0", fontsize=11)
    ax.set_title("Fraction of reads exactly matching reference barcode  (HD = 0)",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    out = os.path.join(STATS_DIR, "hd0_summary.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved hd0_summary.png")


# ── markdown report ───────────────────────────────────────────────────────────

def write_report(all_data):
    label_key     = {SAMPLES[k]["label"]: k for k in SAMPLES}
    keys_present  = [label_key[l] for l in SAMPLE_ORDER if label_key.get(l) in all_data["BC1"]]
    labels_present = [SAMPLES[k]["label"] for k in keys_present]

    md = [
        "# Hamming Distance Report",
        "",
        "**Experiment:** control test — all barcodes expected to match the reference.",
        "",
        "| Barcode | Reference sequence |",
        "|---------|-------------------|",
        f"| BC1 | `{BC1_REF}` |",
        f"| BC2 | `{BC2_REF}` |",
        f"| CB  | `{CB_REF}` |",
        "",
        "---",
        "",
        "## 1. HD = 0 Summary",
        "",
        "| Sample | BC1 HD=0 | BC1 HD=0 % | BC2 HD=0 | BC2 HD=0 % | CB HD=0 | CB HD=0 % |",
        "|--------|----------|------------|----------|------------|---------|----------|",
    ]

    for k in keys_present:
        label = SAMPLES[k]["label"]
        row = [f"**{label}**"]
        for btype in ("BC1", "BC2", "CB"):
            hd_dist, total = all_data[btype].get(k, ({}, 0))
            hd0 = hd_dist.get(0, 0)
            pct = hd0 / total * 100 if total else 0
            row += [f"{hd0:,}", f"{pct:.2f}%"]
        md.append("| " + " | ".join(row) + " |")

    md += ["", "---", "", "## 2. Full HD Distribution", ""]

    for btype, ref in (("BC1", BC1_REF), ("BC2", BC2_REF), ("CB", CB_REF)):
        md += [f"### {btype}  (ref: `{ref}`)", ""]

        # determine max HD across all samples for this btype
        max_hd = max(
            (max(all_data[btype][k][0].keys())
             for k in keys_present
             if all_data[btype].get(k, ({},))[0]),
            default=0,
        )

        header = ["HD"] + labels_present
        md.append("| " + " | ".join(header) + " |")
        md.append("|" + "|".join(["---"] * len(header)) + "|")

        for hd in range(max_hd + 1):
            row = [str(hd)]
            for k in keys_present:
                hd_dist, total = all_data[btype].get(k, ({}, 0))
                cnt = hd_dist.get(hd, 0)
                pct = cnt / total * 100 if total else 0
                row.append(f"{cnt:,} ({pct:.1f}%)")
            md.append("| " + " | ".join(row) + " |")
        md.append("")

    md += [
        "---",
        "",
        "## 3. Plots",
        "",
        "### HD = 0 across all samples",
        "",
        "![HD=0 summary](hd0_summary.png)",
        "",
        "### BC1 HD distribution per sample",
        "",
        "![BC1 HD dist](hd_dist_BC1.png)",
        "",
        "### BC2 HD distribution per sample",
        "",
        "![BC2 HD dist](hd_dist_BC2.png)",
        "",
        "### CB HD distribution per sample",
        "",
        "![CB HD dist](hd_dist_CB.png)",
    ]

    path = os.path.join(STATS_DIR, "hamming_report.md")
    with open(path, "w") as f:
        f.write("\n".join(md))
    print("Saved hamming_report.md")


# ── main ─────────────────────────────────────────────────────────────────────

# build label→key lookup used in plot helpers
SAMPLES_BY_LABEL = {v["label"]: k for k, v in SAMPLES.items()}

if __name__ == "__main__":
    all_data = {"BC1": {}, "BC2": {}, "CB": {}}

    for k in SAMPLES:
        label = SAMPLES[k]["label"]
        print(f"[{label}] loading ...", flush=True)
        all_data["BC1"][k] = load_hd_dist(k, "BC1", BC1_REF)
        all_data["BC2"][k] = load_hd_dist(k, "BC2", BC2_REF)
        all_data["CB"][k]  = load_hd_dist(k, "CB",  CB_REF)

    print("\nGenerating plots ...")
    for btype, ref in (("BC1", BC1_REF), ("BC2", BC2_REF), ("CB", CB_REF)):
        plot_hd_per_sample(btype, ref, {k: all_data[btype][k] for k in SAMPLES})

    plot_hd0_summary(all_data)
    write_report(all_data)
    print("Done.")
