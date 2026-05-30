# Read structure

## Coordinate rule

1-based closed coordinates.

R1[1:10] means bases 1–10.
R1[11:15] means bases 11–15.
R1[51:] means from base 51 to the end.

## Read1

Structure:

BC1 -- W1 -- BC2 -- UMI_2N -- common_fixed -- capture -- insert -- optional_read2_primer

Segments:

| Segment | Coordinate | Length | Expected sequence / note |
|---|---|---:|---|
| BC1 | R1[1:10] | 10 bp | droplet barcode part 1 |
| W1 | R1[11:15] | 5 bp | TCGAG |
| BC2 | R1[16:25] | 10 bp | droplet barcode part 2 |
| UMI_2N | R1[26:27] | 2 bp | random 2N |
| common_fixed | R1[28:35] | 8 bp | TAAGGCGA |
| capture | R1[36:50] | 15 bp | CACCGTCTCCGCCTC |
| insert | R1[51:] | variable | DNA insert; may contain Read2 primer at the end |

## Read2

Structure:

insert -- optional_read1_primer

Segments:

| Segment | Coordinate | Length | Expected sequence / note |
|---|---|---:|---|
| insert | R2[1:] | variable | DNA insert |
| optional_read1_primer | unknown | variable | may appear at the end of Read2 |