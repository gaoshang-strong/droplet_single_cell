#!/usr/bin/env python3
import os
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

STEP_DIR = "/ShangGaoAIProjects/ZhangJW/R1_parser/round2/step1_anchor_exact_hit"
ANCHOR   = "TAAGGCGACACCGTCTCCGCCTC"

SAMPLE_LABELS = {
    "260430R-S-XY-3PY": "3PY",
    "260430R-S-XY-1PB": "1PB",
    "260430R-S-XY-2PB": "2PB",
    "260430R-S-XY-6PY": "6PY",
    "260430R-S-XY-9PY": "9PY",
    "260430R-S-XY-3PB": "3PB",
}
ORDER  = ["3PB", "1PB", "2PB", "3PY", "6PY", "9PY"]
COLORS = ["#2196F3", "#4CAF50", "#FF9800", "#9C27B0", "#F44336", "#00BCD4"]
cmap   = dict(zip(ORDER, COLORS))

# Round1 15bp exact hit rates for comparison
R1_RATES = {"3PB": 95.2, "1PB": 86.4, "2PB": 90.0, "3PY": 82.7, "6PY": 60.9, "9PY": 63.0}

def extract_key(path):
    for key in SAMPLE_LABELS:
        if key in path:
            return key
    return None

summaries = {}
pos_data  = {}

for f in glob.glob(os.path.join(STEP_DIR, "*_anchor_summary.tsv")):
    key = extract_key(f)
    if not key: continue
    d = {}
    with open(f) as fh:
        for line in fh:
            k, v = line.strip().split("\t")
            d[k] = v
    label = SAMPLE_LABELS[key]
    summaries[label] = {
        "total": int(d["total_read1"]),
        "hit":   int(d["anchor_exact_hit"]),
        "rate":  float(d["anchor_exact_hit_rate"]) * 100,
    }

for f in glob.glob(os.path.join(STEP_DIR, "*_anchor_position_distribution.sorted.tsv")):
    key = extract_key(f)
    if not key: continue
    pos, rate = [], []
    with open(f) as fh:
        next(fh)
        for line in fh:
            p = line.strip().split("\t")
            pos.append(int(p[0])); rate.append(float(p[2]))
    pos_data[SAMPLE_LABELS[key]] = (np.array(pos), np.array(rate))

x = np.arange(len(ORDER))

# --- Plot 1: Hit rate comparison (round1 15bp vs round2 23bp) ---
fig, ax = plt.subplots(figsize=(10, 5))
w = 0.35
r1 = [R1_RATES[s] for s in ORDER]
r2 = [summaries[s]["rate"] for s in ORDER]
b1 = ax.bar(x - w/2, r1, width=w, color="#90CAF9", label="Round1 — 15 bp (capture only)", alpha=0.9)
b2 = ax.bar(x + w/2, r2, width=w, color=[cmap[s] for s in ORDER], label="Round2 — 23 bp (common_fixed + capture)", alpha=0.9)
for bar, v in zip(b1, r1):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
            f"{v:.1f}%", ha="center", va="bottom", fontsize=8, color="#1565C0")
for bar, v in zip(b2, r2):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
            f"{v:.1f}%", ha="center", va="bottom", fontsize=8, color="#333")
ax.set_xticks(x); ax.set_xticklabels(ORDER, fontsize=11)
ax.set_ylim(0, 115)
ax.set_ylabel("Exact hit rate (%)", fontsize=10)
ax.set_title("Exact Hit Rate: Round1 (15 bp) vs Round2 (23 bp)", fontsize=12, fontweight="bold")
ax.legend(fontsize=9)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))
plt.tight_layout()
plt.savefig(os.path.join(STEP_DIR, "hit_rate_comparison.png"), dpi=150)
plt.close(); print("Saved hit_rate_comparison.png")

# --- Plot 2: Position distribution zoom pos 20-35 ---
fig, axes = plt.subplots(2, 3, figsize=(14, 7))
for i, label in enumerate(ORDER):
    ax = axes[i//3][i%3]
    pos, rate = pos_data[label]
    mask = (pos >= 20) & (pos <= 35)
    ax.bar(pos[mask], rate[mask]*100, width=0.7, color=cmap[label], alpha=0.85)
    ax.set_title(label, fontsize=12, fontweight="bold", color=cmap[label])
    ax.set_xlabel("Position (1-based)", fontsize=9)
    ax.set_ylabel("Rate among hits (%)", fontsize=9)
    ax.set_xticks(range(20, 36))
    ax.axvline(x=28, color="red", linestyle="--", linewidth=0.9, alpha=0.7)
    ax.text(28.15, ax.get_ylim()[1]*0.88, "pos 28", color="red", fontsize=8)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
fig.suptitle("23 bp Anchor Position Distribution (zoom pos 20–35)", fontsize=13, fontweight="bold")
plt.tight_layout()
plt.savefig(os.path.join(STEP_DIR, "position_distribution_zoom.png"), dpi=150)
plt.close(); print("Saved position_distribution_zoom.png")

# --- Markdown ---
md = []
md.append("# Round2 Step1 — 23 bp Combined Anchor Exact Hit")
md.append("")
md.append(f"**Anchor:** `{ANCHOR}` (23 bp)  ")
md.append("**Composition:** `TAAGGCGA` (common_fixed, 8 bp) + `CACCGTCTCCGCCTC` (capture, 15 bp)  ")
md.append("**Expected position:** starts at pos 28 (1-based)  ")
md.append("**Method:** exact string match on R1 reads")
md.append("")
md.append("---")
md.append("")
md.append("## 1. Exact Hit Rate — Round1 vs Round2")
md.append("")
md.append("![Hit Rate Comparison](hit_rate_comparison.png)")
md.append("")
md.append("| Sample | Round1 15 bp | Round2 23 bp | Drop |")
md.append("|--------|-------------|-------------|------|")
for s in ORDER:
    r1v = R1_RATES[s]; r2v = summaries[s]["rate"]
    md.append(f"| **{s}** | {r1v:.1f}% | {r2v:.1f}% | -{r1v-r2v:.1f}% |")
md.append("")
md.append("---")
md.append("")
md.append("## 2. 23 bp Anchor Position Distribution")
md.append("")
md.append("![Position Distribution](position_distribution_zoom.png)")
md.append("")
md.append("All samples show a dominant peak at **position 28**, confirming that when the 23 bp anchor is found, it is at the correct expected location.")
md.append("")
md.append("---")
md.append("")
md.append("## 3. Summary Table")
md.append("")
md.append("| Sample | Total R1 | Exact Hit | Hit Rate |")
md.append("|--------|----------|-----------|----------|")
for s in ORDER:
    d = summaries[s]
    md.append(f"| **{s}** | {d['total']:,} | {d['hit']:,} | {d['rate']:.2f}% |")
md.append("")
md.append("---")
md.append("")
md.append("## 4. Conclusion — Why Round2 Was Abandoned")
md.append("")
md.append("Using the 23 bp combined anchor reveals a systematic batch difference:")
md.append("")
md.append("| Batch | Samples | Round2 Hit Rate |")
md.append("|-------|---------|----------------|")
md.append("| PB | 1PB, 2PB, 3PB | 58–69% |")
md.append("| PY | 3PY, 6PY, 9PY | 26–37% |")
md.append("")
md.append("The ~30% gap between PB and PY batches — far exceeding what sequencing error can explain — strongly suggests that **PY and PB samples carry different 8 bp fixed sequences at pos 28–35** (analogous to different Illumina i7 sample indices). The `TAAGGCGA` common_fixed is present in PB reads but likely replaced by a different sequence in PY reads.")
md.append("")
md.append("**Decision:** Revert to the Round1 strategy using only the 15 bp capture sequence `CACCGTCTCCGCCTC` as anchor, which is consistent across all samples. The `common_fixed` region will not be used as an anchor going forward.")

md_path = os.path.join(STEP_DIR, "round2_step1_report.md")
with open(md_path, "w") as f:
    f.write("\n".join(md))
print(f"Saved {md_path}")
