"""
Add CB and UMI tags to BAM by looking up read name in a pre-built TSV map.

TSV format (no header): read_name  CB  UMI
Writes SAM tags: CB:Z:{CB}  UB:Z:{UMI}
"""

import sys
import gzip
import argparse
import pysam


def load_cb_map(tsv_gz):
    cb_map = {}
    with gzip.open(tsv_gz, 'rt') as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) == 3:
                cb_map[parts[0]] = (parts[1], parts[2])
    print(f"Loaded {len(cb_map):,} barcode mappings", file=sys.stderr)
    return cb_map


def process(bam_in, bam_out, cb_map):
    tagged = untagged = 0
    with pysam.AlignmentFile(bam_in, 'rb') as bam, \
         pysam.AlignmentFile(bam_out, 'wb', header=bam.header) as out:
        for read in bam:
            entry = cb_map.get(read.query_name)
            if entry:
                read.set_tag('CB', entry[0])
                read.set_tag('UB', entry[1])
                tagged += 1
            else:
                untagged += 1
            out.write(read)

    print(f"Tagged: {tagged:,}  No barcode found: {untagged:,}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--bam-in',  required=True)
    ap.add_argument('--bam-out', required=True)
    ap.add_argument('--cb-map',  required=True, help='TSV.gz from extract_cb_map.py')
    args = ap.parse_args()

    cb_map = load_cb_map(args.cb_map)
    process(args.bam_in, args.bam_out, cb_map)


if __name__ == '__main__':
    main()
