"""
Generate a standard fragment file from a deduplicated, CB-tagged BAM.

Output format (5 columns, 0-based half-open, sorted by chr/start):
  chr  start  end  barcode  count

- Only proper pairs where both reads are mapped are used.
- Fragment boundaries: min(r1.start, r2.start), max(r1.end, r2.end)
- Only barcodes in the whitelist are written (pass --whitelist to filter).
- Output is bgzip-compressed; call tabix separately to index.
"""

import sys
import argparse
import gzip
import pysam
from collections import defaultdict


def load_whitelist(path):
    with open(path) as fh:
        return set(line.strip() for line in fh if line.strip())


def process(bam_path, out_path, whitelist=None):
    frags = defaultdict(int)

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for read in bam:
            # Only process read1 of each proper pair to avoid double-counting
            if not read.is_proper_pair or read.is_unmapped or read.mate_is_unmapped:
                continue
            if not read.is_read1:
                continue

            try:
                cb = read.get_tag('CB')
            except KeyError:
                continue

            if whitelist and cb not in whitelist:
                continue

            chrom = read.reference_name
            start = min(read.reference_start, read.next_reference_start)
            end   = start + abs(read.template_length)

            frags[(chrom, start, end, cb)] += 1

    # Sort by chr, start
    sorted_frags = sorted(frags.items(), key=lambda x: (x[0][0], x[0][1]))

    with gzip.open(out_path, 'wt') as fh:
        for (chrom, start, end, cb), count in sorted_frags:
            fh.write(f"{chrom}\t{start}\t{end}\t{cb}\t{count}\n")

    print(f"Fragments written: {len(sorted_frags):,}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--bam', required=True)
    ap.add_argument('--out', required=True, help='Output .tsv.gz path')
    ap.add_argument('--whitelist', default=None,
                    help='One CB per line; only these cells are included')
    args = ap.parse_args()

    wl = load_whitelist(args.whitelist) if args.whitelist else None
    process(args.bam, args.out, wl)


if __name__ == '__main__':
    main()
