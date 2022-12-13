# PySmooth

### Input File format

The First row is the header.

The genotype file MUST have the following columns:

- Column 1: Chromosome name.
- Column 2: Genomic Position of the marker in the chromosome. For each chromosome,column2 MUST already be sorted in ascending order.
- Column 3: identification number of the marker location.
- Column 4: Reference allele in the reference or can be left blank cell.
- Column 5: Alternate allele if known or blank cell.
- Column 6 and beyond: Genotype code for the individuals. Four codes can be used. A: parent 1 homozygous, B: parent 2 homozygous, H: heterozygous, U: missing data.

### Running PySmooth

PySmooth takes the following arguments

- `-i` or `--input`: Name of the input genotype file. This MUST be provided
- '-o' or '--output': Prefix to name of output files to be generated. If not provided, default is `test`
- `-c` or `--chr`: list of chromosome names to perform analysis on. Names should be separated by comma (e.g `chr1,chr2,chr3`). Default is to run through all the chromosomes in the genotype file.
- `-l` or `--lower`: Lowest threshold for identifying singletons. Default is 0.70
- `-u`or `--upper`: Highest threshold for identifying singletons. Default is 0.98.
- `-g` or `--gap`: PySmooth iteratively identifies singletons starting with the highest threshold till the lowest threshold. This parameter is used to decreased the threshold at each iteration. Default is 0.02.

### Outputs

For each chromosome, PySmooth Generates the following outputs.

- One statistics file that contains `%` of homozygous, heterozygous calls for each individual for the raw genoytpe file.
- One statistics file that contains `%` of homozygous, heterozygous calls for each individual after removing singletons from the raw genoytpe file. Singletons are marked missing.
- One summary statistics File that indicate how many singletons were detected. 
- One Heatmap image that displays genotype calls for the raw genoytpe file. 
- One Heatmap image that displays genotype calls for after removing singletons from the raw genoytpe file. 
- new Genotype File with singletons marked as missing. Singletons are marked as 'U'.

### Running PySmooth

`python run_smooth.py -i my_genotype_file.csv`


`python run_smooth.py -i my_genotype_file.csv -o my_output -c chr1,chr2,chr3 -l 0.80 -u 0.98 -g 0.02`



