import pysam
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import defaultdict

BAM_DIR = "/ShangGaoAIProjects/ZhangJW/Alignment/round1/step2b_mouse_genome"
OUT_FIG  = "/ShangGaoAIProjects/ZhangJW/Alignment/round1/atac_coverage_manhattan.png"
BIN_SIZE = 50_000

SAMPLES = ["1PB", "2PB", "3PB", "3PY", "6PY", "9PY"]
BAMS = {s: f"{BAM_DIR}/260430R-S-XY-{s}.bowtie2.sorted.bam" for s in SAMPLES}

# keep only main chromosomes (chr1-19, chrX, chrY)
MAIN_CHROMS = [f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY"]
CHROM_COLORS = ["#4878CF", "#6ACC65"] * 11  # alternating blue/green

def get_chrom_sizes(bam_path):
    bam = pysam.AlignmentFile(bam_path, "rb")
    sizes = {s["SN"]: s["LN"] for s in bam.header["SQ"]}
    bam.close()
    return sizes

def count_bins(bam_path, chrom_sizes):
    counts = {}
    for chrom in MAIN_CHROMS:
        if chrom not in chrom_sizes:
            continue
        n_bins = chrom_sizes[chrom] // BIN_SIZE + 1
        counts[chrom] = np.zeros(n_bins, dtype=np.int32)

    bam = pysam.AlignmentFile(bam_path, "rb")
    for read in bam.fetch():
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        chrom = read.reference_name
        if chrom not in counts:
            continue
        b = read.reference_start // BIN_SIZE
        if b < len(counts[chrom]):
            counts[chrom][b] += 1
    bam.close()
    return counts

def build_x_y(counts, chrom_sizes):
    xs, ys, colors, chrom_ticks = [], [], [], []
    offset = 0
    for i, chrom in enumerate(MAIN_CHROMS):
        if chrom not in counts:
            continue
        n = len(counts[chrom])
        x = np.arange(n) + offset
        y = np.log10(counts[chrom].astype(float) + 1)
        xs.append(x)
        ys.append(y)
        colors.extend([CHROM_COLORS[i]] * n)
        chrom_ticks.append((offset + n / 2, chrom.replace("chr", "")))
        offset += n + 5  # small gap between chromosomes
    return (np.concatenate(xs), np.concatenate(ys),
            np.array(colors), chrom_ticks, offset)

# ── main ──────────────────────────────────────────────────────────────
print("Reading chromosome sizes...")
chrom_sizes = get_chrom_sizes(BAMS["1PB"])

all_counts = {}
for sample in SAMPLES:
    print(f"  Counting bins: {sample}")
    all_counts[sample] = count_bins(BAMS[sample], chrom_sizes)

fig, axes = plt.subplots(6, 1, figsize=(18, 18), sharex=False)
fig.suptitle("ATAC-seq Coverage (GRCm39 · 50 kb bins · Bowtie2)",
             fontsize=14, fontweight="bold", y=0.98)

for ax, sample in zip(axes, SAMPLES):
    x, y, colors, chrom_ticks, total_width = build_x_y(all_counts[sample], chrom_sizes)
    ax.scatter(x, y, c=colors, s=0.3, rasterized=True, linewidths=0)
    ax.set_ylabel("log₁₀(reads + 1)", fontsize=8)
    ax.set_title(sample, fontsize=10, fontweight="bold", loc="left", pad=2)
    ax.set_xlim(-10, total_width)
    ax.set_ylim(0, None)
    ax.set_xticks([t for t, _ in chrom_ticks])
    ax.set_xticklabels([l for _, l in chrom_ticks], fontsize=6)
    ax.tick_params(axis="y", labelsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig(OUT_FIG, dpi=150, bbox_inches="tight")
print(f"Saved: {OUT_FIG}")
