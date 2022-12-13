# PySmooth

### Input File format

The First row is the header.

The genotype file MUST have the following columns:

- Column 1: chromosome name
- Column 2: Genomic Position in the chromosome. For each chromosome,column2 MUST already be sorted in ascending order.
- Column 3: identification number of the location
- Column 4: Reference allele in the reference or can be left blank cell.
- Column 5: Alternate allele if known or blank cell.
- Column 6 and beyond: Genotype code for the individuals. Four codes can be used. A: parent 1 homozygous, B: parent2 homozygous, H: heterozygous, U: missing data.


