#!/usr/bin/env python3
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

STEP2_DIR = "/ShangGaoAIProjects/ZhangJW/R1_parser/round1/step2_use_local_alignment_for_reads_without_exact_hit"

# Load combined summary
summary_file = os.path.join(STEP2_DIR, "all_samples_hamming_summary.tsv")
samples, no_hit, h1, h2, h3, hgt3 = [], [], [], [], [], []

with open(summary_file) as f:
    next(f)
    for line in f:
        parts = line.strip().split("\t")
        samples.append(parts[0])
        no_hit.append(int(parts[1]))
        h1.append(float(parts[3]))
        h2.append(float(parts[5]))
        h3.append(float(parts[7]))
        hgt3.append(float(parts[9]))

# Also load exact hit rates from step1 summaries
STEP1_DIR = "/ShangGaoAIProjects/ZhangJW/R1_parser/round1/step1_use_capture_Seq_to_anchor_reads_structure"
SAMPLE_KEYS = {
    "3PY": "260430R-S-XY-3PY",
    "1PB": "260430R-S-XY-1PB",
    "2PB": "260430R-S-XY-2PB",
    "6PY": "260430R-S-XY-6PY",
    "9PY": "260430R-S-XY-9PY",
    "3PB": "260430R-S-XY-3PB",
}
import glob
exact_hit_rates = {}
for f in glob.glob(os.path.join(STEP1_DIR, "*_anchor_summary.tsv")):
    for label, key in SAMPLE_KEYS.items():
        if key in f:
            with open(f) as fh:
                data = dict(line.strip().split("\t") for line in fh)
            exact_hit_rates[label] = float(data["anchor_exact_hit_rate"]) * 100

colors = ["#43A047", "#FDD835", "#FB8C00", "#E53935"]
labels_cat = ["Hamming = 1", "Hamming = 2", "Hamming = 3", "Hamming > 3"]

n = len(samples)
x = np.arange(n)

# --- Plot 1: Stacked bar - proportion of hamming categories among no-hit reads ---
fig, ax = plt.subplots(figsize=(9, 5))
b1 = ax.bar(x, h1, color=colors[0], label=labels_cat[0])
b2 = ax.bar(x, h2, bottom=h1, color=colors[1], label=labels_cat[1])
b3 = ax.bar(x, h3, bottom=np.array(h1)+np.array(h2), color=colors[2], label=labels_cat[2])
b4 = ax.bar(x, hgt3, bottom=np.array(h1)+np.array(h2)+np.array(h3), color=colors[3], label=labels_cat[3])

ax.set_xticks(x)
ax.set_xticklabels(samples, fontsize=11)
ax.set_ylabel("Proportion among no-exact-hit reads (%)", fontsize=10)
ax.set_ylim(0, 108)
ax.set_title("Min Hamming Distance Distribution\n(reads without exact anchor hit)", fontsize=12, fontweight="bold")
ax.legend(loc="upper right", fontsize=9)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))

# annotate hamming=1 bar
for i, v in enumerate(h1):
    ax.text(i, v / 2, f"{v:.1f}%", ha="center", va="center", fontsize=8.5,
            color="white", fontweight="bold")

plt.tight_layout()
plt.savefig(os.path.join(STEP2_DIR, "hamming_stacked_bar.png"), dpi=150)
plt.close()
print("Saved hamming_stacked_bar.png")

# --- Plot 2: Grouped bar ---
fig, ax = plt.subplots(figsize=(11, 5))
width = 0.2
for i, (cat_data, color, label) in enumerate(zip([h1, h2, h3, hgt3], colors, labels_cat)):
    ax.bar(x + i * width - 1.5 * width, cat_data, width=width,
           color=color, label=label, alpha=0.9)
ax.set_xticks(x)
ax.set_xticklabels(samples, fontsize=11)
ax.set_ylabel("Proportion (%)", fontsize=10)
ax.set_title("Min Hamming Distance per Category per Sample", fontsize=12, fontweight="bold")
ax.legend(fontsize=9)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))
plt.tight_layout()
plt.savefig(os.path.join(STEP2_DIR, "hamming_grouped_bar.png"), dpi=150)
plt.close()
print("Saved hamming_grouped_bar.png")

# --- Plot 3: Recovery overview - exact hit + hamming 1/2/3 vs lost ---
fig, ax = plt.subplots(figsize=(9, 5))
exact = [exact_hit_rates.get(s, 0) for s in samples]
# no_hit_pct = 100 - exact
# among no_hit: h1/h2/h3 are already percentages of no_hit reads
no_hit_pct = [100 - e for e in exact]
rec1 = [no_hit_pct[i] * h1[i] / 100 for i in range(n)]
rec2 = [no_hit_pct[i] * h2[i] / 100 for i in range(n)]
rec3 = [no_hit_pct[i] * h3[i] / 100 for i in range(n)]
lost = [no_hit_pct[i] * hgt3[i] / 100 for i in range(n)]

c_exact  = "#1565C0"
c_rec    = ["#43A047", "#FDD835", "#FB8C00"]
c_lost   = "#E53935"

ax.bar(x, exact,  color=c_exact,   label="Exact hit")
ax.bar(x, rec1,   bottom=exact,                         color=c_rec[0], label="No-hit, Hamming=1")
ax.bar(x, rec2,   bottom=np.array(exact)+np.array(rec1), color=c_rec[1], label="No-hit, Hamming=2")
ax.bar(x, rec3,   bottom=np.array(exact)+np.array(rec1)+np.array(rec2), color=c_rec[2], label="No-hit, Hamming=3")
ax.bar(x, lost,   bottom=np.array(exact)+np.array(rec1)+np.array(rec2)+np.array(rec3), color=c_lost,  label="No-hit, Hamming>3")

ax.set_xticks(x)
ax.set_xticklabels(samples, fontsize=11)
ax.set_ylabel("% of all R1 reads", fontsize=10)
ax.set_ylim(0, 108)
ax.set_title("Read Recovery Overview\n(relative to all R1 reads)", fontsize=12, fontweight="bold")
ax.legend(loc="lower right", fontsize=8.5)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))
plt.tight_layout()
plt.savefig(os.path.join(STEP2_DIR, "recovery_overview.png"), dpi=150)
plt.close()
print("Saved recovery_overview.png")

# --- Write Markdown ---
md = []
md.append("# Anchor Hamming Distance Analysis (Step 2)")
md.append("")
md.append("**Anchor:** `CACCGTCTCCGCCTC` (15 bp)  ")
md.append("**Scope:** reads without exact anchor hit  ")
md.append("**Method:** sliding window (length 15, pos 33–60, 1-based), minimum Hamming distance  ")
md.append("")
md.append("---")
md.append("")
md.append("## 1. No-exact-hit Read Counts")
md.append("")
md.append("| Sample | No-exact-hit Reads |")
md.append("|--------|--------------------|")
for i, s in enumerate(samples):
    md.append(f"| **{s}** | {no_hit[i]:,} |")
md.append("")
md.append("---")
md.append("")
md.append("## 2. Min Hamming Distance Distribution (among no-hit reads)")
md.append("")
md.append("### 2.1 Stacked bar")
md.append("")
md.append("![Stacked Bar](hamming_stacked_bar.png)")
md.append("")
md.append("### 2.2 Grouped bar")
md.append("")
md.append("![Grouped Bar](hamming_grouped_bar.png)")
md.append("")
md.append("### 2.3 Summary table")
md.append("")
md.append("| Sample | No-hit Reads | Hamming=1 (%) | Hamming=2 (%) | Hamming=3 (%) | Hamming>3 (%) |")
md.append("|--------|-------------|--------------|--------------|--------------|--------------|")
for i, s in enumerate(samples):
    md.append(f"| **{s}** | {no_hit[i]:,} | {h1[i]:.2f}% | {h2[i]:.2f}% | {h3[i]:.2f}% | {hgt3[i]:.2f}% |")
md.append("")
md.append("---")
md.append("")
md.append("## 3. Read Recovery Overview (relative to all R1 reads)")
md.append("")
md.append("![Recovery Overview](recovery_overview.png)")
md.append("")
md.append("| Sample | Exact Hit | +Hamming=1 | +Hamming=2 | +Hamming=3 | Lost (>3) |")
md.append("|--------|-----------|-----------|-----------|-----------|----------|")
for i, s in enumerate(samples):
    md.append(f"| **{s}** | {exact[i]:.1f}% | {rec1[i]:.1f}% | {rec2[i]:.1f}% | {rec3[i]:.1f}% | {lost[i]:.1f}% |")
md.append("")
md.append("---")
md.append("")
md.append("## 4. Observations")
md.append("")
md.append("- **2PB** has the highest Hamming=1 proportion (71.7%) among no-hit reads, suggesting most misses are single-base sequencing errors — recoverable with 1-mismatch tolerance.")
md.append("- **6PY** and **9PY** have the largest Hamming>3 fractions (~46–50%), consistent with their low exact hit rates and suggesting genuine structural issues in a large fraction of reads.")
md.append("- **3PY** also shows elevated Hamming>3 (38%), despite a higher exact hit rate than 6PY/9PY.")
md.append("- Allowing up to 1 mismatch would recover an additional 20–72% of no-hit reads depending on sample, substantially increasing usable read count.")

md_path = os.path.join(STEP2_DIR, "hamming_analysis_report.md")
with open(md_path, "w") as f:
    f.write("\n".join(md))
print(f"Saved {md_path}")
