"""
Pure-compute scanning functions for round3.
No I/O — import and call from run_scan.py or test directly.

All returned coordinates are 1-based closed intervals.
"""

from config import CAPTURE_SEQ, W1_SEQ, CAPTURE_SCAN_WINDOWS

# W1 is searched only within 1-based pos 1–35 (0-based 0:35)
W1_SEARCH_END = 35


def _hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def scan_capture(seq: str) -> tuple[int, int, int]:
    """
    Scan CAPTURE_SCAN_WINDOWS and return the window with minimum Hamming
    distance to CAPTURE_SEQ.

    Returns (best_start, best_end, min_hamming) in 1-based closed coords.
    Ties broken by the first (smallest-start) window.
    """
    best_start = best_end = best_hd = None

    for start, end in CAPTURE_SCAN_WINDOWS:
        window = seq[start - 1 : end]          # convert to 0-based slice
        if len(window) < len(CAPTURE_SEQ):
            continue                            # read too short for this window
        hd = _hamming(window, CAPTURE_SEQ)
        if best_hd is None or hd < best_hd:
            best_hd    = hd
            best_start = start
            best_end   = end

    return best_start, best_end, best_hd


def scan_W1(seq: str) -> tuple[bool, int | None, int | None]:
    """
    Search for an exact match of W1_SEQ within 1-based pos 1–W1_SEARCH_END.

    Returns (exact_hit, start, end) where start/end are 1-based closed,
    or (False, None, None) if not found.

    If multiple matches exist, returns the leftmost one.
    """
    region = seq[:W1_SEARCH_END]
    p = region.find(W1_SEQ)
    if p == -1:
        return False, None, None
    start = p + 1                               # convert to 1-based
    end   = start + len(W1_SEQ) - 1
    return True, start, end
