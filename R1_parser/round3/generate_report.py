"""
Generate QC report for round3 metadata CSVs.
Stats are cached as JSON so re-runs skip re-reading large files.
"""
import os
import json
import subprocess
import multiprocessing

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

from config import SAMPLES, OUTPUT_DIR

SAMPLE_ORDER = ["1PB", "2PB", "3PB", "3PY", "6PY", "9PY"]
KEY_BY_LABEL = {v["label"]: k for k, v in SAMPLES.items()}

C_PB = "#1565C0"
C_PY = "#E65100"

def sample_color(label):
    return C_PB if label.endswith("PB") else C_PY


# ── streaming stats collection ───────────────────────────────────────────────

def collect_stats(key: str) -> dict:
    label    = SAMPLES[key]["label"]
    csv_path = os.path.join(OUTPUT_DIR, f"{key}_metadata.csv.gz")
    cache    = os.path.join(OUTPUT_DIR, f"{key}_report_stats.json")

    if os.path.exists(cache):
        print(f"[{label}] loaded cached stats", flush=True)
        with open(cache) as f:
            return json.load(f)

    print(f"[{label}] scanning ...", flush=True)

    total         = 0
    hd_dist       = {}     # hamming_distance → count
    cs_dist       = {}     # capture_seq_start → count
    w1_hit        = 0
    w1_start_dist = {}     # W1_start → count  (when hit=TRUE)
    gap_dist      = {}     # gap_length → count (when numeric)

    # combined filter counts
    n_hd0          = 0     # hd == 0
    n_hd0_w1       = 0     # hd == 0 AND W1=TRUE
    n_hd0_w1_gap20 = 0     # hd == 0 AND W1=TRUE AND gap == 20

    proc = subprocess.Popen(
        ["pigz", "-dc", "-p", "2", csv_path], stdout=subprocess.PIPE
    )
    proc.stdout.readline()   # skip header

    for raw in proc.stdout:
        parts = raw.decode().rstrip("\n").split(",")
        if len(parts) < 8:
            continue
        total += 1

        # columns: read_name, cs_start, cs_end, hd, W1_hit, W1_start, W1_end, gap
        hd      = int(parts[3])
        hit     = parts[4] == "TRUE"
        ws_raw  = parts[5]
        gap_raw = parts[7]
        cs_s    = int(parts[1])

        hd_dist[hd]   = hd_dist.get(hd, 0) + 1
        cs_dist[cs_s] = cs_dist.get(cs_s, 0) + 1

        if hit:
            w1_hit += 1
            ws = int(ws_raw)
            w1_start_dist[ws] = w1_start_dist.get(ws, 0) + 1

        if gap_raw != "NA":
            g = int(gap_raw)
            gap_dist[g] = gap_dist.get(g, 0) + 1

        if hd == 0:
            n_hd0 += 1
            if hit:
                n_hd0_w1 += 1
                if gap_raw != "NA" and int(gap_raw) == 20:
                    n_hd0_w1_gap20 += 1

    proc.wait()

    stats = {
        "label":          label,
        "key":            key,
        "total":          total,
        "hd_dist":        {str(k): v for k, v in hd_dist.items()},
        "cs_dist":        {str(k): v for k, v in cs_dist.items()},
        "w1_hit":         w1_hit,
        "w1_start_dist":  {str(k): v for k, v in w1_start_dist.items()},
        "gap_dist":       {str(k): v for k, v in gap_dist.items()},
        "n_hd0":          n_hd0,
        "n_hd0_w1":       n_hd0_w1,
        "n_hd0_w1_gap20": n_hd0_w1_gap20,
    }

    with open(cache, "w") as f:
        json.dump(stats, f, indent=2)
    print(f"[{label}] done  total={total:,}", flush=True)
    return stats


# ── plotting helpers ─────────────────────────────────────────────────────────

def save(fig, name):
    path = os.path.join(OUTPUT_DIR, name)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {name}")


def plot_hamming(all_stats):
    """Stacked bar: hamming distance breakdown per sample."""
    labels = [s["label"] for s in all_stats]
    n      = len(labels)
    x      = np.arange(n)
    tiers  = [(0, "#1B5E20"), (1, "#43A047"), (2, "#FDD835"),
              (3, "#FB8C00"), (4, "#E53935")]
    tier_labels = ["HD=0", "HD=1", "HD=2", "HD=3", "HD≥4"]

    pcts = []
    for tier, _ in tiers:
        row = []
        for s in all_stats:
            t  = s["total"]
            if tier < 4:
                c = s["hd_dist"].get(str(tier), 0)
            else:
                c = sum(v for k, v in s["hd_dist"].items() if int(k) >= 4)
            row.append(c / t * 100)
        pcts.append(row)

    fig, ax = plt.subplots(figsize=(9, 5))
    bottom = np.zeros(n)
    for (_, color), lbl, row in zip(tiers, tier_labels, pcts):
        vals = np.array(row)
        bars = ax.bar(x, vals, bottom=bottom, color=color, label=lbl, width=0.55)
        for i, (b, v) in enumerate(zip(bottom, vals)):
            if v > 1.5:
                ax.text(i, b + v / 2, f"{v:.1f}%",
                        ha="center", va="center", fontsize=8,
                        color="white", fontweight="bold")
        bottom += vals

    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylim(0, 108)
    ax.set_ylabel("% of reads", fontsize=10)
    ax.set_title("Capture sequence Hamming distance distribution", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9, loc="upper right")
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    save(fig, "hamming_dist.png")


def plot_capture_start(all_stats):
    """Stacked bar: which window (36/37/38/39) is chosen as best match."""
    labels  = [s["label"] for s in all_stats]
    n       = len(labels)
    x       = np.arange(n)
    windows = [36, 37, 38, 39]
    colors  = ["#1565C0", "#1E88E5", "#90CAF9", "#BBDEFB"]

    fig, ax = plt.subplots(figsize=(9, 5))
    bottom  = np.zeros(n)
    for w, c in zip(windows, colors):
        vals = np.array([
            s["cs_dist"].get(str(w), 0) / s["total"] * 100
            for s in all_stats
        ])
        ax.bar(x, vals, bottom=bottom, color=c, label=f"start={w}", width=0.55)
        for i, (b, v) in enumerate(zip(bottom, vals)):
            if v > 1.5:
                ax.text(i, b + v / 2, f"{v:.1f}%",
                        ha="center", va="center", fontsize=8,
                        color="white", fontweight="bold")
        bottom += vals

    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylim(0, 108)
    ax.set_ylabel("% of reads", fontsize=10)
    ax.set_title("Best-match capture window start position", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    save(fig, "capture_start_dist.png")


def plot_w1_hit(all_stats):
    """Bar: W1 exact hit rate per sample."""
    labels = [s["label"] for s in all_stats]
    n      = len(labels)
    x      = np.arange(n)
    rates  = [s["w1_hit"] / s["total"] * 100 for s in all_stats]
    colors = [sample_color(l) for l in labels]

    fig, ax = plt.subplots(figsize=(8, 4.5))
    bars = ax.bar(x, rates, color=colors, width=0.55)
    for bar, v in zip(bars, rates):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.4,
                f"{v:.1f}%", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.axhline(90, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylim(0, 108)
    ax.set_ylabel("% of reads", fontsize=10)
    ax.set_title("W1 (TCGAG) exact hit rate in pos 1–35", fontsize=12, fontweight="bold")
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    save(fig, "W1_hit_rate.png")


def plot_w1_start(all_stats):
    """Grouped bar: W1_start distribution (11/12/13/14/other) among W1-hit reads."""
    labels   = [s["label"] for s in all_stats]
    n        = len(labels)
    x        = np.arange(n)
    pos_list = [11, 12, 13, 14]
    colors   = ["#1565C0", "#1E88E5", "#42A5F5", "#90CAF9", "#E53935"]
    w        = 0.15
    offsets  = np.arange(-2, 3) * w

    fig, ax = plt.subplots(figsize=(10, 5))
    for pos, c, off in zip(pos_list + ["other"], colors, offsets):
        vals = []
        for s in all_stats:
            denom = s["w1_hit"] or 1
            if pos == "other":
                cnt = denom - sum(
                    s["w1_start_dist"].get(str(p), 0) for p in pos_list
                )
            else:
                cnt = s["w1_start_dist"].get(str(pos), 0)
            vals.append(cnt / denom * 100)
        lbl = f"pos {pos}" if pos != "other" else "other"
        ax.bar(x + off, vals, width=w, color=c, label=lbl)

    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylabel("% of W1-hit reads", fontsize=10)
    ax.set_title("W1 start position distribution (among W1 exact-hit reads)",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    save(fig, "W1_start_dist.png")


def plot_gap(all_stats):
    """Grouped bar: gap_length distribution (19/20/21/other) among W1-hit reads."""
    labels   = [s["label"] for s in all_stats]
    n        = len(labels)
    x        = np.arange(n)
    gaps     = [19, 20, 21]
    colors   = ["#FDD835", "#1B5E20", "#FB8C00", "#E53935"]
    w        = 0.18
    offsets  = np.arange(-1.5, 2.5) * w

    fig, ax = plt.subplots(figsize=(10, 5))
    for g, c, off in zip(gaps + ["other"], colors, offsets):
        vals = []
        for s in all_stats:
            denom = s["w1_hit"] or 1
            if g == "other":
                cnt = denom - sum(
                    s["gap_dist"].get(str(gv), 0) for gv in gaps
                )
                # subtract NA (no W1) which is already excluded from gap_dist
                # gap_dist only contains numeric gap values from W1-hit reads
                # but actually gap_dist includes ALL reads with numeric gap
                # recalculate: use w1_hit as denominator and gap_dist for numerics
                cnt = s["w1_hit"] - sum(s["gap_dist"].get(str(gv), 0) for gv in gaps)
            else:
                cnt = s["gap_dist"].get(str(g), 0)
            vals.append(max(cnt, 0) / denom * 100)
        lbl = f"gap={g}" if g != "other" else "other"
        ax.bar(x + off, vals, width=w, color=c, label=lbl)

    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylabel("% of W1-hit reads", fontsize=10)
    ax.set_title("Gap length between W1 end and capture start\n(among W1 exact-hit reads)",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    save(fig, "gap_dist.png")


def plot_filter_funnel(all_stats):
    """Grouped bar: cumulative pass rates for three filter tiers."""
    labels = [s["label"] for s in all_stats]
    n      = len(labels)
    x      = np.arange(n)
    w      = 0.25
    tiers  = [
        ("HD=0",                  lambda s: s["n_hd0"],          "#43A047"),
        ("HD=0 + W1 hit",         lambda s: s["n_hd0_w1"],       "#1565C0"),
        ("HD=0 + W1 hit + gap=20",lambda s: s["n_hd0_w1_gap20"], "#6A1B9A"),
    ]

    fig, ax = plt.subplots(figsize=(10, 5))
    for i, (lbl, fn, c) in enumerate(tiers):
        vals = [fn(s) / s["total"] * 100 for s in all_stats]
        off  = (i - 1) * w
        bars = ax.bar(x + off, vals, width=w, color=c, label=lbl)
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
                    f"{v:.1f}%", ha="center", va="bottom", fontsize=7.5,
                    fontweight="bold", rotation=90)

    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylim(0, 115)
    ax.set_ylabel("% of raw reads", fontsize=10)
    ax.set_title("Cumulative filter pass rates", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    save(fig, "filter_funnel.png")


# ── markdown report ──────────────────────────────────────────────────────────

def write_report(all_stats):
    by_label = {s["label"]: s for s in all_stats}

    md = []
    md += [
        "# Round3 Metadata QC Report",
        "",
        "**Input:** raw R1 reads (all 6 samples)  ",
        "**Capture scan windows (1-based closed):** [36,50] [37,51] [38,52] [39,53]  ",
        "**W1 search region:** pos 1–35  ",
        "",
        "---",
        "",
        "## 1. Capture Sequence Hamming Distance",
        "",
        "Minimum Hamming distance across the four scan windows.",
        "HD=0 means at least one window has an exact match.",
        "",
        "![Hamming dist](hamming_dist.png)",
        "",
    ]

    md += [
        "| Sample | HD=0 | HD=1 | HD=2 | HD=3 | HD≥4 |",
        "|--------|------|------|------|------|------|",
    ]
    for label in SAMPLE_ORDER:
        s  = by_label[label]
        t  = s["total"]
        hd = s["hd_dist"]
        h  = [hd.get(str(i), 0) for i in range(4)]
        h4 = sum(v for k, v in hd.items() if int(k) >= 4)
        md.append(
            f"| **{label}** | {h[0]/t*100:.1f}% | {h[1]/t*100:.1f}% |"
            f" {h[2]/t*100:.1f}% | {h[3]/t*100:.1f}% | {h4/t*100:.1f}% |"
        )

    md += [
        "",
        "---",
        "",
        "## 2. Best-Match Capture Window Start Position",
        "",
        "Which of the four candidate windows (pos 36/37/38/39) yielded the minimum Hamming distance.",
        "",
        "![Capture start](capture_start_dist.png)",
        "",
    ]

    md += [
        "| Sample | start=36 | start=37 | start=38 | start=39 |",
        "|--------|----------|----------|----------|----------|",
    ]
    for label in SAMPLE_ORDER:
        s  = by_label[label]
        t  = s["total"]
        cs = s["cs_dist"]
        md.append(
            f"| **{label}** |"
            + " | ".join(f" {cs.get(str(w), 0)/t*100:.1f}%" for w in [36,37,38,39])
            + " |"
        )

    md += [
        "",
        "---",
        "",
        "## 3. W1 Exact Hit Rate",
        "",
        "Exact match of `TCGAG` in pos 1–35 of R1.",
        "",
        "![W1 hit rate](W1_hit_rate.png)",
        "",
    ]

    md += [
        "| Sample | W1 hit | W1 miss | Hit rate |",
        "|--------|--------|---------|----------|",
    ]
    for label in SAMPLE_ORDER:
        s    = by_label[label]
        t    = s["total"]
        hit  = s["w1_hit"]
        miss = t - hit
        md.append(
            f"| **{label}** | {hit:,} | {miss:,} | **{hit/t*100:.1f}%** |"
        )

    md += [
        "",
        "---",
        "",
        "## 4. W1 Start Position Distribution",
        "",
        "Among reads with a W1 exact hit, which position the W1 was found at.",
        "Canonical position is 11; positions 12–14 indicate a 1–3 bp prefix before BC1.",
        "",
        "![W1 start dist](W1_start_dist.png)",
        "",
    ]

    md += [
        "| Sample | pos 11 | pos 12 | pos 13 | pos 14 | other |",
        "|--------|--------|--------|--------|--------|-------|",
    ]
    for label in SAMPLE_ORDER:
        s     = by_label[label]
        denom = s["w1_hit"] or 1
        wsd   = s["w1_start_dist"]
        vals  = [wsd.get(str(p), 0) for p in [11,12,13,14]]
        other = denom - sum(vals)
        md.append(
            f"| **{label}** |"
            + "".join(f" {v/denom*100:.1f}% |" for v in vals)
            + f" {other/denom*100:.1f}% |"
        )

    md += [
        "",
        "---",
        "",
        "## 5. Gap Length Distribution",
        "",
        "Number of bases between W1 end and capture start, among W1-hit reads.  ",
        "Expected value: **20** (BC2 10bp + UMI 2bp + common_fixed 8bp).",
        "",
        "![Gap dist](gap_dist.png)",
        "",
    ]

    md += [
        "| Sample | gap=19 | gap=20 | gap=21 | other |",
        "|--------|--------|--------|--------|-------|",
    ]
    for label in SAMPLE_ORDER:
        s     = by_label[label]
        denom = s["w1_hit"] or 1
        gd    = s["gap_dist"]
        vals  = [gd.get(str(g), 0) for g in [19, 20, 21]]
        other = denom - sum(vals)
        md.append(
            f"| **{label}** |"
            + "".join(f" {v/denom*100:.1f}% |" for v in vals)
            + f" {max(other,0)/denom*100:.1f}% |"
        )

    md += [
        "",
        "---",
        "",
        "## 6. Cumulative Filter Pass Rates",
        "",
        "Three progressively stricter filter tiers, each expressed as % of raw reads.",
        "",
        "![Filter funnel](filter_funnel.png)",
        "",
    ]

    md += [
        "| Sample | HD=0 | HD=0 + W1 hit | HD=0 + W1 hit + gap=20 |",
        "|--------|------|---------------|------------------------|",
    ]
    for label in SAMPLE_ORDER:
        s = by_label[label]
        t = s["total"]
        md.append(
            f"| **{label}** | {s['n_hd0']/t*100:.1f}% ({s['n_hd0']:,})"
            f" | {s['n_hd0_w1']/t*100:.1f}% ({s['n_hd0_w1']:,})"
            f" | **{s['n_hd0_w1_gap20']/t*100:.1f}%** ({s['n_hd0_w1_gap20']:,}) |"
        )

    path = os.path.join(OUTPUT_DIR, "qc_report.md")
    with open(path, "w") as f:
        f.write("\n".join(md))
    print(f"Saved qc_report.md")


# ── main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    keys = list(SAMPLES.keys())

    with multiprocessing.Pool(processes=len(keys)) as pool:
        all_stats = pool.map(collect_stats, keys)

    label_map = {s["label"]: s for s in all_stats}
    all_stats  = [label_map[l] for l in SAMPLE_ORDER]

    print("\nGenerating plots ...")
    plot_hamming(all_stats)
    plot_capture_start(all_stats)
    plot_w1_hit(all_stats)
    plot_w1_start(all_stats)
    plot_gap(all_stats)
    plot_filter_funnel(all_stats)

    write_report(all_stats)
    print("Done.")
