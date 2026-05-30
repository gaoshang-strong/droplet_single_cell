#!/usr/bin/env python3
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

STEP6_DIR = "/ShangGaoAIProjects/ZhangJW/R1_parser/round1/step6_filter_reads_with_W1"
summary_file = os.path.join(STEP6_DIR, "W1_filter_summary.tsv")

samples, totals, passed, passed_pct, failed, failed_pct = [], [], [], [], [], []

with open(summary_file) as f:
    next(f)
    for line in f:
        p = line.strip().split("\t")
        samples.append(p[0])
        totals.append(int(p[1]))
        passed.append(int(p[2]))
        passed_pct.append(float(p[3]))
        failed.append(int(p[4]))
        failed_pct.append(float(p[5]))

n       = len(samples)
x       = np.arange(n)
totals  = np.array(totals)
passed  = np.array(passed)
failed  = np.array(failed)
passed_pct = np.array(passed_pct)
failed_pct = np.array(failed_pct)

C_PASS   = "#1565C0"
C_FAILED = "#E53935"

# --- Plot 1: Pass rate per sample ---
fig, ax = plt.subplots(figsize=(8, 4.5))
bar_colors = [C_PASS if p >= 95 else C_FAILED for p in passed_pct]
bars = ax.bar(x, passed_pct, color=bar_colors, width=0.55)
for bar, v in zip(bars, passed_pct):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
            f"{v:.1f}%", ha="center", va="bottom", fontsize=10, fontweight="bold")
ax.axhline(y=95, color="gray", linestyle="--", linewidth=0.8, alpha=0.6, label="95% threshold")
ax.set_xticks(x)
ax.set_xticklabels(samples, fontsize=11)
ax.set_ylim(0, 108)
ax.set_ylabel("Pass rate (%)", fontsize=10)
ax.set_title("W1 (TCGAG) Exact Match Pass Rate per Sample", fontsize=12, fontweight="bold")
ax.legend(fontsize=9)
ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(STEP6_DIR, "W1_pass_rate.png"), dpi=150)
plt.close()
print("Saved W1_pass_rate.png")

# --- Plot 2: Stacked bar — proportion ---
fig, ax = plt.subplots(figsize=(9, 5))
ax.bar(x, passed_pct, color=C_PASS,   label="W1 hit (passed)")
ax.bar(x, failed_pct, bottom=passed_pct, color=C_FAILED, label="No W1 hit (failed)")
for i in range(n):
    ax.text(i, passed_pct[i] / 2, f"{passed_pct[i]:.1f}%",
            ha="center", va="center", fontsize=9, color="white", fontweight="bold")
    ax.text(i, passed_pct[i] + failed_pct[i] / 2, f"{failed_pct[i]:.1f}%",
            ha="center", va="center", fontsize=9, color="white", fontweight="bold")
ax.set_xticks(x)
ax.set_xticklabels(samples, fontsize=11)
ax.set_ylabel("% of total reads", fontsize=10)
ax.set_ylim(0, 108)
ax.set_title("W1 Filtering Result per Sample", fontsize=13, fontweight="bold")
ax.legend(loc="upper right", fontsize=9)
ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(STEP6_DIR, "W1_stacked_bar.png"), dpi=150)
plt.close()
print("Saved W1_stacked_bar.png")

# --- Plot 3: Absolute counts stacked bar ---
fig, ax = plt.subplots(figsize=(9, 5))
ax.bar(x, passed, color=C_PASS,   label="W1 hit (passed)")
ax.bar(x, failed, bottom=passed,  color=C_FAILED, label="No W1 hit (failed)")
ax.set_xticks(x)
ax.set_xticklabels(samples, fontsize=11)
ax.set_ylabel("Read count", fontsize=10)
ax.set_title("W1 Filtering — Absolute Counts", fontsize=13, fontweight="bold")
ax.legend(loc="upper right", fontsize=9)
ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v/1e6:.0f}M"))
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(STEP6_DIR, "W1_count_bar.png"), dpi=150)
plt.close()
print("Saved W1_count_bar.png")

# --- Write Markdown ---
md = []
md.append("# W1 Filtering Report (Step 6)")
md.append("")
md.append("**Input:** step4-filtered reads (`*_filtered_R1/R2.fq.gz`)")
md.append("")
md.append("**Filter criterion:** exact match of W1 (`TCGAG`, 5 bp) in the region before the capture sequence start (`cs:i` position)")
md.append("")
md.append("**Output:** `*_W1_R1/R2.fq.gz`; W1 start position written to R1 read name as `w1:i:{pos}` (1-based)")
md.append("")
md.append("---")
md.append("")
md.append("## 1. Pass Rate per Sample")
md.append("")
md.append("![Pass Rate](W1_pass_rate.png)")
md.append("")
md.append("---")
md.append("")
md.append("## 2. Filtering Breakdown")
md.append("")
md.append("### 2.1 Proportion of total reads")
md.append("")
md.append("![Stacked Bar](W1_stacked_bar.png)")
md.append("")
md.append("### 2.2 Absolute read counts")
md.append("")
md.append("![Count Bar](W1_count_bar.png)")
md.append("")
md.append("---")
md.append("")
md.append("## 3. Summary Table")
md.append("")
md.append("| Sample | Total Reads | Passed (W1 hit) | Failed | Pass Rate |")
md.append("|--------|------------|-----------------|--------|-----------|")
for i, s in enumerate(samples):
    md.append(
        f"| **{s}** | {totals[i]:,} | {passed[i]:,} ({passed_pct[i]:.2f}%) "
        f"| {failed[i]:,} ({failed_pct[i]:.2f}%) "
        f"| **{passed_pct[i]:.2f}%** |"
    )
md.append("")
md.append("---")
md.append("")
md.append("## 4. Observations")
md.append("")
md.append("- **PB batch** (1PB, 2PB, 3PB): W1 pass rates 96–98%, consistent with step5 exact hit survey.")
md.append("- **PY batch** (3PY, 6PY, 9PY): slightly lower at 92–94%, but still high enough for reliable barcode extraction.")
md.append("- The ~4–6% gap between PB and PY batches is consistent across all upstream steps and likely reflects library quality differences.")
md.append("- The `w1:i` tag in each passing read's name records the 1-based W1 start position alongside the existing `cs:i` and `mt:Z` tags, enabling downstream barcode extraction without re-scanning.")

md_path = os.path.join(STEP6_DIR, "W1_filter_report.md")
with open(md_path, "w") as f:
    f.write("\n".join(md))
print(f"Saved {md_path}")
