# Read structure

## Coordinate rule

Python-style 0-based half-open coordinates.

R1[0:10] means bases 0-9.
R1[10:15] means bases 10-14.
R1[50:] means from base 50 to the end.

## Read1

Structure:

BC1 -- W1 -- BC2 -- UMI_2N -- common_fixed -- capture -- insert -- optional_read2_primer

Segments:

| Segment | Coordinate | Length | Expected sequence / note |
|---|---|---:|---|
| BC1 | R1[0:10] | 10 bp | droplet barcode part 1 |
| W1 | R1[10:15] | 5 bp | TCGAG |
| BC2 | R1[15:25] | 10 bp | droplet barcode part 2 |
| UMI_2N | R1[25:27] | 2 bp | random 2N |
| common_fixed | R1[27:35] | 8 bp | TAAGGCGA |
| capture | R1[35:50] | 15 bp | CACCGTCTCCGCCTC |
| insert | R1[50:] | variable | DNA insert; may contain Read2 primer at the end |

## Read2

Structure:

insert -- optional_read1_primer

Segments:

| Segment | Coordinate | Length | Expected sequence / note |
|---|---|---:|---|
| insert | R2[:] | variable | DNA insert |
| optional_read1_primer | unknown | variable | may appear at the end of Read2 |