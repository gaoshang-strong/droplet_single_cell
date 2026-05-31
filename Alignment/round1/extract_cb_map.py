"""
Scan R1 FASTQ from filter_HD2 and write read_name -> CB+UMI mapping TSV.

Output (tab-separated, no header):
  read_name  CB(20bp)  UMI(2bp)
"""

import re
import gzip
import argparse
import sys

BC1_RE = re.compile(r'bc1:Z:(\S+)')
BC2_RE = re.compile(r'bc2:Z:(\S+)')
UMI_RE = re.compile(r'umi:Z:(\S+)')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--r1', required=True, help='filter_HD2 R1 FASTQ (gz)')
    ap.add_argument('--out', required=True, help='Output .tsv.gz')
    args = ap.parse_args()

    written = skipped = 0
    with gzip.open(args.r1, 'rt') as fh, gzip.open(args.out, 'wt') as out:
        while True:
            header = fh.readline()
            if not header:
                break
            fh.readline()   # seq
            fh.readline()   # +
            fh.readline()   # qual

            qname = header[1:].split()[0]
            m1 = BC1_RE.search(header)
            m2 = BC2_RE.search(header)
            mu = UMI_RE.search(header)

            if m1 and m2 and mu:
                cb = m1.group(1) + m2.group(1)
                out.write(f"{qname}\t{cb}\t{mu.group(1)}\n")
                written += 1
            else:
                skipped += 1

    print(f"Written: {written:,}  Skipped: {skipped:,}", file=sys.stderr)


if __name__ == '__main__':
    main()
