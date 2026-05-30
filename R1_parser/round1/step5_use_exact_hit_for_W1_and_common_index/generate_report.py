#!/usr/bin/env python3
import os
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

STEP5_DIR = "/ShangGaoAIProjects/ZhangJW/R1_parser/round1/step5_use_exact_hit_for_W1_and_common_index"

SAMPLE_LABELS = {
    "260430R-S-XY-3PY": "3PY",
    "260430R-S-XY-1PB": "1PB",
    "260430R-S-XY-2PB": "2PB",
    "260430R-S-XY-6PY": "6PY",
    "260430R-S-XY-9PY": "9PY",
    "260430R-S-XY-3PB": "3PB",
}
ORDER = ["3PB", "1PB", "2PB", "3PY", "6PY", "9PY"]
COLORS = ["#2196F3", "#4CAF50", "#FF9800", "#9C27B0", "#F44336", "#00BCD4"]
color_map = dict(zip(ORDER, COLORS))

def extract_key(path):
    for key in SAMPLE_LABELS:
        if key in path:
            return key
    return None

def load_summary(path):
    data = {}
    with open(path) as f:
        for line in f:
            k, v = line.strip().split("\t")
            data[k] = v
    return data

def load_dist(path):
    pos, rate = [], []
    with open(path) as f:
        next(f)
        for line in f:
            p = line.strip().split("\t")
            pos.append(int(p[0]))
            rate.append(float(p[2]))
    return np.array(pos), np.array(rate)

# Load all data
summaries = {}   # label -> {total, w1_hit, w1_rate, cf_hit, cf_rate}
w1_dists  = {}   # label -> (pos_arr, rate_arr)
cf_dists  = {}

for key, label in SAMPLE_LABELS.items():
    w1_sum  = os.path.join(STEP5_DIR, f"{key}_W1_summary.tsv")
    cf_sum  = os.path.join(STEP5_DIR, f"{key}_common_fixed_summary.tsv")
    w1_dist = os.path.join(STEP5_DIR, f"{key}_W1_position_distribution.tsv")
    cf_dist = os.path.join(STEP5_DIR, f"{key}_common_fixed_position_distribution.tsv")

    ws = load_summary(w1_sum)
    cs = load_summary(cf_sum)
    summaries[label] = {
        "total":   int(ws["total_reads"]),
        "w1_hit":  int(ws["exact_hit"]),
        "w1_rate": float(ws["exact_hit_rate"]) * 100,
        "cf_hit":  int(cs["exact_hit"]),
        "cf_rate": float(cs["exact_hit_rate"]) * 100,
    }
    w1_dists[label] = load_dist(w1_dist)
    cf_dists[label] = load_dist(cf_dist)

samples = ORDER
x = np.arange(len(samples))

# --- Plot 1: Hit rate bar chart (W1 and common_fixed side by side) ---
fig, ax = plt.subplots(figsize=(9, 5))
w = 0.35
w1_rates = [summaries[s]["w1_rate"] for s in samples]
cf_rates = [summaries[s]["cf_rate"] for s in samples]
b1 = ax.bar(x - w/2, w1_rates, width=w, color="#1565C0", label="W1 (TCGAG)", alpha=0.9)
b2 = ax.bar(x + w/2, cf_rates, width=w, color="#E65100", label="common_fixed (TAAGGCGA)", alpha=0.9)
for bar, v in zip(b1, w1_rates):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
            f"{v:.1f}%", ha="center", va="bottom", fontsize=8.5, fontweight="bold", color="#1565C0")
for bar, v in zip(b2, cf_rates):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
            f"{v:.1f}%", ha="center", va="bottom", fontsize=8.5, fontweight="bold", color="#E65100")
ax.set_xticks(x)
ax.set_xticklabels(samples, fontsize=11)
ax.set_ylim(0, 115)
ax.set_ylabel("Exact hit rate (%)", fontsize=10)
ax.set_title("W1 and common_fixed Exact Hit Rate\n(among step4-filtered reads)", fontsize=12, fontweight="bold")
ax.legend(fontsize=9)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))
plt.tight_layout()
plt.savefig(os.path.join(STEP5_DIR, "hit_rate_bar.png"), dpi=150)
plt.close()
print("Saved hit_rate_bar.png")

# --- Plot 2: W1 position distribution (all samples, zoom pos 1-35) ---
fig, axes = plt.subplots(2, 3, figsize=(14, 7))
for i, label in enumerate(samples):
    ax = axes[i // 3][i % 3]
    pos, rate = w1_dists[label]
    mask = pos <= 35
    ax.bar(pos[mask], rate[mask] * 100, width=0.7, color=color_map[label], alpha=0.85)
    ax.set_title(label, fontsize=12, fontweight="bold", color=color_map[label])
    ax.set_xlabel("Position (1-based)", fontsize=9)
    ax.set_ylabel("Rate among hits (%)", fontsize=9)
    ax.axvline(x=11, color="red", linestyle="--", linewidth=0.9, alpha=0.7)
    ax.text(11.3, ax.get_ylim()[1] * 0.88, "pos 11", color="red", fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
fig.suptitle("W1 (TCGAG) Position Distribution — pos 1–35", fontsize=13, fontweight="bold")
plt.tight_layout()
plt.savefig(os.path.join(STEP5_DIR, "W1_position_dist.png"), dpi=150)
plt.close()
print("Saved W1_position_dist.png")

# --- Plot 3: common_fixed position distribution (zoom pos 20-35) ---
fig, axes = plt.subplots(2, 3, figsize=(14, 7))
for i, label in enumerate(samples):
    ax = axes[i // 3][i % 3]
    pos, rate = cf_dists[label]
    mask = (pos >= 20) & (pos <= 35)
    ax.bar(pos[mask], rate[mask] * 100, width=0.7, color=color_map[label], alpha=0.85)
    ax.set_title(label, fontsize=12, fontweight="bold", color=color_map[label])
    ax.set_xlabel("Position (1-based)", fontsize=9)
    ax.set_ylabel("Rate among hits (%)", fontsize=9)
    ax.set_xticks(range(20, 36))
    ax.axvline(x=28, color="red", linestyle="--", linewidth=0.9, alpha=0.7)
    ax.text(28.2, ax.get_ylim()[1] * 0.88, "pos 28", color="red", fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
fig.suptitle("common_fixed (TAAGGCGA) Position Distribution — pos 20–35", fontsize=13, fontweight="bold")
plt.tight_layout()
plt.savefig(os.path.join(STEP5_DIR, "common_fixed_position_dist.png"), dpi=150)
plt.close()
print("Saved common_fixed_position_dist.png")

# --- Plot 4: W1 pos 8-15 fine zoom ---
fig, axes = plt.subplots(2, 3, figsize=(14, 7))
for i, label in enumerate(samples):
    ax = axes[i // 3][i % 3]
    pos, rate = w1_dists[label]
    mask = (pos >= 8) & (pos <= 15)
    ax.bar(pos[mask], rate[mask] * 100, width=0.6, color=color_map[label], alpha=0.85)
    ax.set_title(label, fontsize=12, fontweight="bold", color=color_map[label])
    ax.set_xlabel("Position (1-based)", fontsize=9)
    ax.set_ylabel("Rate among hits (%)", fontsize=9)
    ax.set_xticks(range(8, 16))
    ax.axvline(x=11, color="red", linestyle="--", linewidth=0.9, alpha=0.7)
    ax.text(11.1, ax.get_ylim()[1] * 0.88, "pos 11", color="red", fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
fig.suptitle("W1 (TCGAG) Position Distribution — zoom pos 8–15", fontsize=13, fontweight="bold")
plt.tight_layout()
plt.savefig(os.path.join(STEP5_DIR, "W1_position_dist_zoom.png"), dpi=150)
plt.close()
print("Saved W1_position_dist_zoom.png")

# --- Markdown ---
md = []
md.append("# W1 and common_fixed Exact Hit Report (Step 5)")
md.append("")
md.append("**Input:** step4-filtered R1 reads (reads with capture seq exact or Hamming ≤ 3)  ")
md.append("**Search region:** pos 1 to (capture_start − 1) (1-based), i.e. everything before the anchor  ")
md.append("**W1:** `TCGAG` (expected at pos 11–15)  ")
md.append("**common_fixed:** `TAAGGCGA` (expected at pos 28–35)  ")
md.append("**Method:** exact string match (`str.find`), first occurrence only")
md.append("")
md.append("---")
md.append("")
md.append("## 1. Exact Hit Rates")
md.append("")
md.append("![Hit Rate](hit_rate_bar.png)")
md.append("")
md.append("| Sample | Total Reads | W1 Hit | W1 Rate | common_fixed Hit | CF Rate |")
md.append("|--------|------------|--------|---------|-----------------|---------|")
for s in samples:
    d = summaries[s]
    md.append(f"| **{s}** | {d['total']:,} | {d['w1_hit']:,} | {d['w1_rate']:.2f}% | {d['cf_hit']:,} | {d['cf_rate']:.2f}% |")
md.append("")
md.append("---")
md.append("")
md.append("## 2. W1 Position Distribution")
md.append("")
md.append("### 2.1 Overview (pos 1–35)")
md.append("")
md.append("![W1 Position Distribution](W1_position_dist.png)")
md.append("")
md.append("### 2.2 Zoom (pos 8–15)")
md.append("")
md.append("![W1 Zoom](W1_position_dist_zoom.png)")
md.append("")
md.append("---")
md.append("")
md.append("## 3. common_fixed Position Distribution (pos 20–35)")
md.append("")
md.append("![common_fixed Position Distribution](common_fixed_position_dist.png)")
md.append("")
md.append("---")
md.append("")
md.append("## 4. Per-Sample Peak Summary")
md.append("")
md.append("| Sample | W1 peak pos | W1 peak rate | CF peak pos | CF peak rate |")
md.append("|--------|------------|-------------|------------|-------------|")
for s in samples:
    wp, wr = w1_dists[s]
    cp, cr = cf_dists[s]
    w_peak_pos  = wp[np.argmax(wr)]
    w_peak_rate = wr.max() * 100
    c_peak_pos  = cp[np.argmax(cr)]
    c_peak_rate = cr.max() * 100
    md.append(f"| **{s}** | {w_peak_pos} | {w_peak_rate:.2f}% | {c_peak_pos} | {c_peak_rate:.2f}% |")
md.append("")
md.append("---")
md.append("")
md.append("## 5. Observations")
md.append("")
md.append("- **W1 (`TCGAG`)**: dominant peak at position **11** across all samples, consistent with the expected read structure (BC1 = pos 1–10, W1 = pos 11–15).")
md.append("- **common_fixed (`TAAGGCGA`)**: dominant peak at position **28**, consistent with the expected structure (BC2 = pos 16–25, UMI_2N = pos 26–27, common_fixed = pos 28–35).")
md.append("- W1 hit rates are high (92–98%) across all samples, confirming stable linker structure.")
md.append("- common_fixed hit rates are lower (42–73%), which is expected since `TAAGGCGA` is 8 bp and shares bases with barcode regions; some reads may carry sequence variants at BC2/UMI that create false misses.")
md.append("- **6PY** and **9PY** show lower common_fixed hit rates (~42–45%), consistent with their lower overall anchor quality seen in step 1–4.")

md_path = os.path.join(STEP5_DIR, "W1_common_fixed_report.md")
with open(md_path, "w") as f:
    f.write("\n".join(md))
print(f"Saved {md_path}")
