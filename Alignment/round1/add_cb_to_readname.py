"""
Embed cell barcode (CB = bc1+bc2) and UMI from R1 FASTQ comment into read name.

R1 comment format:
  1:N:0:... cs:i:36 hd:i:0 w1s:i:11 w1e:i:15 gap:i:20 bc1:Z:GATCGATCGA bc2:Z:CTAGCTAGCT umi:Z:AG

Output QNAME format (both R1 and R2):
  {original_qname}#{CB}#{UMI}
  e.g. AE01-...:0000:0463#GATCGATCGACTAGCTAGCT#AG
"""

import sys
import re
import gzip
import argparse
from pathlib import Path


BC1_RE = re.compile(r'bc1:Z:(\S+)')
BC2_RE = re.compile(r'bc2:Z:(\S+)')
UMI_RE = re.compile(r'umi:Z:(\S+)')


def open_fastq(path):
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def iter_records(fh):
    while True:
        header = fh.readline()
        if not header:
            break
        seq  = fh.readline().rstrip('\n')
        plus = fh.readline()
        qual = fh.readline().rstrip('\n')
        yield header.rstrip('\n'), seq, qual


def parse_barcode(r1_header):
    """Return (cb, umi) from R1 header comment, or (None, None) if missing."""
    m1 = BC1_RE.search(r1_header)
    m2 = BC2_RE.search(r1_header)
    mu = UMI_RE.search(r1_header)
    if m1 and m2 and mu:
        return m1.group(1) + m2.group(1), mu.group(1)
    return None, None


def rename_header(header, cb, umi):
    """Replace QNAME with QNAME#CB#UMI, keep the rest of the comment."""
    at_name, _, comment = header.partition(' ')
    new_name = f"{at_name}#{cb}#{umi}"
    return f"{new_name} {comment}" if comment else new_name


def process(r1_in, r2_in, r1_out, r2_out):
    skipped = 0
    written = 0

    with open_fastq(r1_in) as f1, open_fastq(r2_in) as f2, \
         gzip.open(r1_out, 'wt') as o1, gzip.open(r2_out, 'wt') as o2:

        for (h1, s1, q1), (h2, s2, q2) in zip(iter_records(f1), iter_records(f2)):
            cb, umi = parse_barcode(h1)
            if cb is None:
                skipped += 1
                continue

            new_h1 = rename_header(h1, cb, umi)
            new_h2 = rename_header(h2, cb, umi)

            o1.write(f"{new_h1}\n{s1}\n+\n{q1}\n")
            o2.write(f"{new_h2}\n{s2}\n+\n{q2}\n")
            written += 1

    print(f"Written: {written:,}  Skipped (no barcode): {skipped:,}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--r1', required=True)
    ap.add_argument('--r2', required=True)
    ap.add_argument('--out-r1', required=True)
    ap.add_argument('--out-r2', required=True)
    args = ap.parse_args()

    process(args.r1, args.r2, args.out_r1, args.out_r2)


if __name__ == '__main__':
    main()
