import os
import glob
import base64
import io
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

def fig_to_base64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")

def md_img(b64):
    return f"![](data:image/png;base64,{b64})"

STEP1_DIR = "/ShangGaoAIProjects/ZhangJW/R1_parser/round1/step1_use_capture_Seq_to_anchor_reads_structure"

SAMPLE_LABELS = {
    "260430R-S-XY-3PY": "3PY",
    "260430R-S-XY-1PB": "1PB",
    "260430R-S-XY-2PB": "2PB",
    "260430R-S-XY-6PY": "6PY",
    "260430R-S-XY-9PY": "9PY",
    "260430R-S-XY-3PB": "3PB",
}

def extract_sample_key(filename):
    for key in SAMPLE_LABELS:
        if key in filename:
            return key
    return None

# Load summary data
summaries = {}
for f in glob.glob(os.path.join(STEP1_DIR, "*_anchor_summary.tsv")):
    key = extract_sample_key(f)
    if key is None:
        continue
    data = {}
    with open(f) as fh:
        for line in fh:
            k, v = line.strip().split("\t")
            data[k] = v
    summaries[key] = {
        "label": SAMPLE_LABELS[key],
        "total": int(data["total_read1"]),
        "hit": int(data["anchor_exact_hit"]),
        "rate": float(data["anchor_exact_hit_rate"]),
    }

# Load position distribution data
pos_data = {}
for f in glob.glob(os.path.join(STEP1_DIR, "*_anchor_position_distribution.sorted.tsv")):
    key = extract_sample_key(f)
    if key is None:
        continue
    positions, rates = [], []
    with open(f) as fh:
        next(fh)
        for line in fh:
            parts = line.strip().split("\t")
            positions.append(int(parts[0]))
            rates.append(float(parts[2]))
    pos_data[key] = {"pos": positions, "rate": rates}

ordered_keys = ["260430R-S-XY-3PB", "260430R-S-XY-3PY",
                "260430R-S-XY-1PB", "260430R-S-XY-2PB",
                "260430R-S-XY-6PY", "260430R-S-XY-9PY"]

colors = ["#2196F3", "#4CAF50", "#FF9800", "#9C27B0", "#F44336", "#00BCD4"]

# --- Plot 1: Hit rate bar chart ---
fig, ax = plt.subplots(figsize=(8, 4.5))
labels = [SAMPLE_LABELS[k] for k in ordered_keys]
rates  = [summaries[k]["rate"] * 100 for k in ordered_keys]
bars = ax.bar(labels, rates, color=colors, width=0.55, edgecolor="white", linewidth=0.8)
for bar, rate in zip(bars, rates):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
            f"{rate:.1f}%", ha="center", va="bottom", fontsize=10, fontweight="bold")
ax.set_ylim(0, 105)
ax.set_ylabel("Anchor Exact Hit Rate (%)", fontsize=11)
ax.set_title("Anchor Exact Hit Rate per Sample", fontsize=13, fontweight="bold")
ax.axhline(y=80, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
b64_hit = fig_to_base64(fig)
plt.savefig(os.path.join(STEP1_DIR, "hit_rate_barplot.png"), dpi=150)
plt.close()

# --- Plot 2: Position distribution - overview (all positions) ---
fig, axes = plt.subplots(2, 3, figsize=(14, 7))
axes = axes.flatten()
for i, key in enumerate(ordered_keys):
    ax = axes[i]
    pd = pos_data[key]
    pos_arr = np.array(pd["pos"])
    rate_arr = np.array(pd["rate"])
    ax.bar(pos_arr, rate_arr * 100, width=0.8, color=colors[i], alpha=0.85)
    ax.set_title(SAMPLE_LABELS[key], fontsize=12, fontweight="bold", color=colors[i])
    ax.set_xlabel("Position", fontsize=9)
    ax.set_ylabel("Rate among hits (%)", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    peak_pos = pos_arr[np.argmax(rate_arr)]
    ax.axvline(x=peak_pos, color="red", linestyle="--", linewidth=0.9, alpha=0.7)
    ax.text(peak_pos + 1, max(rate_arr) * 95, f"pos {peak_pos}",
            color="red", fontsize=8)
fig.suptitle("Anchor Position Distribution (all positions)", fontsize=14, fontweight="bold")
plt.tight_layout()
b64_overview = fig_to_base64(fig)
plt.savefig(os.path.join(STEP1_DIR, "position_distribution_overview.png"), dpi=150)
plt.close()

# --- Plot 3: Position distribution - zoom in pos 30-45 ---
fig, axes = plt.subplots(2, 3, figsize=(14, 7))
axes = axes.flatten()
for i, key in enumerate(ordered_keys):
    ax = axes[i]
    pd = pos_data[key]
    pos_arr = np.array(pd["pos"])
    rate_arr = np.array(pd["rate"])
    mask = (pos_arr >= 30) & (pos_arr <= 45)
    ax.bar(pos_arr[mask], rate_arr[mask] * 100, width=0.6, color=colors[i], alpha=0.85)
    ax.set_title(SAMPLE_LABELS[key], fontsize=12, fontweight="bold", color=colors[i])
    ax.set_xlabel("Position", fontsize=9)
    ax.set_ylabel("Rate among hits (%)", fontsize=9)
    ax.set_xticks(range(30, 46))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
fig.suptitle("Anchor Position Distribution (zoom: pos 30–45)", fontsize=14, fontweight="bold")
plt.tight_layout()
b64_zoom = fig_to_base64(fig)
plt.savefig(os.path.join(STEP1_DIR, "position_distribution_zoom.png"), dpi=150)
plt.close()

# --- Write Markdown ---
md_lines = []
md_lines.append("# Anchor Sequence QC Report")
md_lines.append("")
md_lines.append(f"**Anchor sequence:** `CACCGTCTCCGCCTC`  ")
md_lines.append(f"**Samples:** {len(ordered_keys)}  ")
md_lines.append(f"**Method:** exact string match on R1 reads")
md_lines.append("")

md_lines.append("---")
md_lines.append("")
md_lines.append("## 1. Anchor Exact Hit Rate")
md_lines.append("")
md_lines.append("![Hit Rate](hit_rate_barplot.png)")
md_lines.append("")
md_lines.append("| Sample | Total R1 Reads | Anchor Hits | Hit Rate |")
md_lines.append("|--------|---------------|-------------|----------|")
for key in ordered_keys:
    s = summaries[key]
    flag = " ⚠️" if s["rate"] < 0.75 else ""
    md_lines.append(f"| **{s['label']}** | {s['total']:,} | {s['hit']:,} | {s['rate']*100:.2f}%{flag} |")
md_lines.append("")
md_lines.append("> ⚠️ = hit rate below 75%, may indicate library quality issue or different read structure.")
md_lines.append("")

md_lines.append("---")
md_lines.append("")
md_lines.append("## 2. Anchor Position Distribution")
md_lines.append("")
md_lines.append("### 2.1 Overview (all positions)")
md_lines.append("")
md_lines.append("![Position Distribution Overview](position_distribution_overview.png)")
md_lines.append("")
md_lines.append("### 2.2 Zoom: positions 30–45 (dominant peak region)")
md_lines.append("")
md_lines.append("![Position Distribution Zoom](position_distribution_zoom.png)")
md_lines.append("")

md_lines.append("---")
md_lines.append("")
md_lines.append("## 3. Per-Sample Peak Summary")
md_lines.append("")
md_lines.append("| Sample | Peak Position | Peak Rate (%) | pos36 Rate (%) | pos37 Rate (%) | pos38 Rate (%) | pos39 Rate (%) |")
md_lines.append("|--------|--------------|--------------|---------------|---------------|---------------|---------------|")
for key in ordered_keys:
    pd = pos_data[key]
    pos_arr = np.array(pd["pos"])
    rate_arr = np.array(pd["rate"])
    peak_idx = np.argmax(rate_arr)
    peak_pos = pos_arr[peak_idx]
    peak_rate = rate_arr[peak_idx] * 100
    def get_rate(p):
        idx = np.where(pos_arr == p)[0]
        return f"{rate_arr[idx[0]]*100:.2f}%" if len(idx) > 0 else "0.00%"
    label = SAMPLE_LABELS[key]
    md_lines.append(f"| **{label}** | {peak_pos} | {peak_rate:.2f}% | {get_rate(36)} | {get_rate(37)} | {get_rate(38)} | {get_rate(39)} |")
md_lines.append("")

md_lines.append("---")
md_lines.append("")
md_lines.append("## 4. Observations")
md_lines.append("")
md_lines.append("- All samples show a dominant peak at **position 36**, consistent with a fixed anchor location in the read structure.")
md_lines.append("- Secondary peaks at positions 37–39 represent minor positional shifts.")
md_lines.append("- **6PY** and **9PY** have notably lower hit rates (~61–63%), which may reflect:")
md_lines.append("  - Higher proportion of reads without the expected library structure")
md_lines.append("  - Sequencing or ligation artifacts")
md_lines.append("  - Different sample preparation batches (B vs A)")

md_path = os.path.join(STEP1_DIR, "anchor_QC_report.md")
with open(md_path, "w") as fh:
    fh.write("\n".join(md_lines))
print(f"Saved {md_path}")
