"""
Round3 pipeline configuration.
All coordinates are 1-based closed intervals unless noted otherwise.
"""

# ── sequences ────────────────────────────────────────────────────────────────

CAPTURE_SEQ = "CACCGTCTCCGCCTC"   # 15 bp
W1_SEQ      = "TCGAG"             # 5 bp

# ── capture scan windows (1-based closed) ────────────────────────────────────
# Each entry is (start, end); length must equal len(CAPTURE_SEQ).
# Represents the four candidate positions where capture seq may start.

CAPTURE_SCAN_WINDOWS = [
    (36, 50),
    (37, 51),
    (38, 52),
    (39, 53),
]

# ── samples ──────────────────────────────────────────────────────────────────

SAMPLES = {
    "260430R-S-XY-1PB": {
        "label": "1PB",
        "R1": (
            "/ShangGaoAIProjects/ZhangJW/data/260430R-S-XY-1PB"
            "/20260511B-ZM/Lane03"
            "/20260509_4P260322158US293276A2_B_260509_260430R-S-XY-1PB_L03"
            "/260430R-S-XY-1PB"
            "/20260509_4P260322158US293276A2_B_260509_260430R-S-XY-1PB_L03_R1.fq.gz"
        ),
    },
    "260430R-S-XY-2PB": {
        "label": "2PB",
        "R1": (
            "/ShangGaoAIProjects/ZhangJW/data/260430R-S-XY-2PB"
            "/20260511B-ZM/Lane04"
            "/20260509_4P260322158US293276A2_B_260509_260430R-S-XY-2PB_L04"
            "/260430R-S-XY-2PB"
            "/20260509_4P260322158US293276A2_B_260509_260430R-S-XY-2PB_L04_R1.fq.gz"
        ),
    },
    "260430R-S-XY-3PB": {
        "label": "3PB",
        "R1": (
            "/ShangGaoAIProjects/ZhangJW/data/260430R-S-XY-3PB"
            "/20260513A-ZM"
            "/20260511_4P260322136US293258A2_A_260511_260430R-S-XY-3PB_L00"
            "/260430R-S-XY-3PB"
            "/20260511_4P260322136US293258A2_A_260511_260430R-S-XY-3PB_L00_R1.fq.gz"
        ),
    },
    "260430R-S-XY-3PY": {
        "label": "3PY",
        "R1": (
            "/ShangGaoAIProjects/ZhangJW/data/260430R-S-XY-3PY"
            "/20260511A-ZM"
            "/20260509_4P260322081US293269A2_A_260509_260430R-S-XY-3PY_L00"
            "/260430R-S-XY-3PY"
            "/20260509_4P260322081US293269A2_A_260509_260430R-S-XY-3PY_L00_R1.fq.gz"
        ),
    },
    "260430R-S-XY-6PY": {
        "label": "6PY",
        "R1": (
            "/ShangGaoAIProjects/ZhangJW/data/260430R-S-XY-6PY"
            "/20260511B-ZM/Lane01"
            "/20260509_4P260322158US293276A2_B_260509_260430R-S-XY-6PY_L01"
            "/260430R-S-XY-6PY"
            "/20260509_4P260322158US293276A2_B_260509_260430R-S-XY-6PY_L01_R1.fq.gz"
        ),
    },
    "260430R-S-XY-9PY": {
        "label": "9PY",
        "R1": (
            "/ShangGaoAIProjects/ZhangJW/data/260430R-S-XY-9PY"
            "/20260511B-ZM/Lane02"
            "/20260509_4P260322158US293276A2_B_260509_260430R-S-XY-9PY_L02"
            "/260430R-S-XY-9PY"
            "/20260509_4P260322158US293276A2_B_260509_260430R-S-XY-9PY_L02_R1.fq.gz"
        ),
    },
}

# Derive R2 paths from R1 (same directory, _R1 → _R2)
for _key in SAMPLES:
    SAMPLES[_key]["R2"] = SAMPLES[_key]["R1"].replace("_R1.fq.gz", "_R2.fq.gz")

# ── filter definitions ───────────────────────────────────────────────────────
# Each entry: (folder_name, max_hamming_distance)
# All filters also require W1_exact_hit=TRUE and gap_length=20.

FILTERS = [
    ("filter_HD0", 0),
    ("filter_HD2", 2),
    ("filter_HD3", 3),
]

# ── paths ────────────────────────────────────────────────────────────────────

OUTPUT_DIR = "/ShangGaoAIProjects/ZhangJW/R1_parser/round3"

# ── runtime ──────────────────────────────────────────────────────────────────

N_WORKERS = 6   # one per sample
