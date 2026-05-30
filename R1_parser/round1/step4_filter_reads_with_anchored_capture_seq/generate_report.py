#!/usr/bin/env python3
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

STEP4_DIR = "/ShangGaoAIProjects/ZhangJW/R1_parser/round1/step4_filter_reads_with_anchored_capture_seq"

summary_file = os.path.join(STEP4_DIR, "filter_summary.tsv")

samples, totals, passed_pct, exact, hamming_r, failed = [], [], [], [], [], []

with open(summary_file) as f:
    next(f)
    for line in f:
        p = line.strip().split("\t")
        samples.append(p[0])
        totals.append(int(p[1]))
        passed_pct.append(float(p[3]))
        exact.append(int(p[4]))
        hamming_r.append(int(p[5]))
        failed.append(int(p[6]))

n = len(samples)
x = np.arange(n)
totals = np.array(totals)
exact_pct    = np.array(exact)    / totals * 100
hamming_pct  = np.array(hamming_r) / totals * 100
failed_pct   = np.array(failed)   / totals * 100

C_EXACT   = "#1565C0"
C_HAMMING = "#43A047"
C_FAILED  = "#E53935"

# --- Plot 1: Stacked bar (% of total reads) ---
fig, ax = plt.subplots(figsize=(9, 5))
b1 = ax.bar(x, exact_pct,   color=C_EXACT,   label="Exact hit")
b2 = ax.bar(x, hamming_pct, bottom=exact_pct, color=C_HAMMING, label="Hamming rescued (≤3)")
b3 = ax.bar(x, failed_pct,  bottom=exact_pct + hamming_pct, color=C_FAILED, label="Failed (Hamming >3)")

for i in range(n):
    ax.text(i, exact_pct[i] / 2, f"{exact_pct[i]:.1f}%",
            ha="center", va="center", fontsize=8, color="white", fontweight="bold")
    if hamming_pct[i] > 1:
        ax.text(i, exact_pct[i] + hamming_pct[i] / 2, f"{hamming_pct[i]:.1f}%",
                ha="center", va="center", fontsize=8, color="white", fontweight="bold")
    ax.text(i, exact_pct[i] + hamming_pct[i] + failed_pct[i] / 2, f"{failed_pct[i]:.1f}%",
            ha="center", va="center", fontsize=8, color="white", fontweight="bold")

ax.set_xticks(x)
ax.set_xticklabels(samples, fontsize=11)
ax.set_ylabel("% of total R1 reads", fontsize=10)
ax.set_ylim(0, 108)
ax.set_title("Read Filtering Result per Sample", fontsize=13, fontweight="bold")
ax.legend(loc="upper right", fontsize=9)
ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(STEP4_DIR, "filter_stacked_bar.png"), dpi=150)
plt.close()
print("Saved filter_stacked_bar.png")

# --- Plot 2: Absolute read counts stacked bar ---
fig, ax = plt.subplots(figsize=(9, 5))
ax.bar(x, exact,    color=C_EXACT,   label="Exact hit")
ax.bar(x, hamming_r, bottom=exact,   color=C_HAMMING, label="Hamming rescued (≤3)")
ax.bar(x, failed,   bottom=np.array(exact) + np.array(hamming_r), color=C_FAILED, label="Failed")

ax.set_xticks(x)
ax.set_xticklabels(samples, fontsize=11)
ax.set_ylabel("Read count", fontsize=10)
ax.set_title("Read Filtering — Absolute Counts", fontsize=13, fontweight="bold")
ax.legend(loc="upper right", fontsize=9)
ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v/1e6:.0f}M"))
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(STEP4_DIR, "filter_count_bar.png"), dpi=150)
plt.close()
print("Saved filter_count_bar.png")

# --- Plot 3: Pass rate comparison ---
fig, ax = plt.subplots(figsize=(8, 4.5))
bar_colors = [C_EXACT if p >= 95 else C_HAMMING if p >= 85 else C_FAILED for p in passed_pct]
bars = ax.bar(x, passed_pct, color=bar_colors, width=0.55)
for bar, v in zip(bars, passed_pct):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
            f"{v:.1f}%", ha="center", va="bottom", fontsize=10, fontweight="bold")
ax.axhline(y=95, color="gray", linestyle="--", linewidth=0.8, alpha=0.6, label="95% threshold")
ax.set_xticks(x)
ax.set_xticklabels(samples, fontsize=11)
ax.set_ylim(0, 108)
ax.set_ylabel("Pass rate (%)", fontsize=10)
ax.set_title("Overall Pass Rate per Sample\n(Exact + Hamming ≤ 3)", fontsize=12, fontweight="bold")
ax.legend(fontsize=9)
ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(STEP4_DIR, "filter_pass_rate.png"), dpi=150)
plt.close()
print("Saved filter_pass_rate.png")

# --- Write Markdown ---
md = []
md.append("# Read Filtering Report (Step 4)")
md.append("")
md.append("**Anchor:** `CACCGTCTCCGCCTC` (15 bp)  ")
md.append("**Filter criteria (pass if either):**")
md.append("1. Exact anchor match anywhere in R1")
md.append("2. Sliding window pos 33–60 (1-based), length 15: min Hamming ≤ 3  ")
md.append("")
md.append("**Output:** filtered R1 + R2 fastq.gz; capture start position written to read name as `cs:i:{pos} mt:Z:{exact|hamming}`")
md.append("")
md.append("---")
md.append("")
md.append("## 1. Pass Rate per Sample")
md.append("")
md.append("![Pass Rate](filter_pass_rate.png)")
md.append("")
md.append("---")
md.append("")
md.append("## 2. Filtering Breakdown")
md.append("")
md.append("### 2.1 Proportion of total reads")
md.append("")
md.append("![Stacked Bar](filter_stacked_bar.png)")
md.append("")
md.append("### 2.2 Absolute read counts")
md.append("")
md.append("![Count Bar](filter_count_bar.png)")
md.append("")
md.append("---")
md.append("")
md.append("## 3. Summary Table")
md.append("")
md.append("| Sample | Total Reads | Exact Hit | Hamming Rescued | Failed | Pass Rate |")
md.append("|--------|------------|-----------|-----------------|--------|-----------|")
for i, s in enumerate(samples):
    md.append(
        f"| **{s}** | {totals[i]:,} | {exact[i]:,} ({exact_pct[i]:.1f}%) "
        f"| {hamming_r[i]:,} ({hamming_pct[i]:.1f}%) "
        f"| {failed[i]:,} ({failed_pct[i]:.1f}%) "
        f"| **{passed_pct[i]:.1f}%** |"
    )
md.append("")
md.append("---")
md.append("")
md.append("## 4. Observations")
md.append("")
md.append("- **2PB** and **3PB** achieve ~99% pass rate — high-quality libraries with consistent anchor placement.")
md.append("- **Hamming rescue** recovers an additional 1.7M–5.2M reads per sample that would otherwise be discarded.")
md.append("- **6PY** and **9PY** have the highest failure rates (~18%), consistent with findings in Step 1 and Step 2.")
md.append("- The `cs:i` tag in each passing read's name records the 1-based anchor start position, enabling downstream barcode/UMI extraction without re-scanning.")

md_path = os.path.join(STEP4_DIR, "filter_report.md")
with open(md_path, "w") as f:
    f.write("\n".join(md))
print(f"Saved {md_path}")
