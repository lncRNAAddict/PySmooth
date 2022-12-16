
# PySmooth

`PySmooth` is a python implementation of SMOOTH algorithm. Original SMOOTH was written in PASCAL with several limitations.

`PySmooth` contains one main functionality 

- `run_smooth.py`: Detect markers which are singletons in the genotype file and marks them as singletons. Correct the missing and singleton markers using k-nearest neighbor.

## Dependencies

- PySmooth has been tested with Python 3.8.12 version. It should work with Python >= 3.0 version.
- `numpy`
- `Pandas`
- `Sklearn`
- `matplotlib`

## Running `run_smooth.py`

### Input Genotype File format

The First row is the header. Each row represents a unique marker.

The genotype file MUST have the following columns:

- Column 1: Chromosome name.
- Column 2: Genomic Position of the marker in the chromosome. For each chromosome,column 2 MUST already be sorted in ascending order.
- Column 3: Identification id of the marker location. 
- Column 4: Reference allele in the reference genome if known or can be left blank cell.
- Column 5: Alternate allele if known or blank cell.
- Column 6 and beyond: Genotype code for the individuals in the marker location. Four codes can be used. A: parent 1 homozygous, B: parent 2 homozygous, H: heterozygous, U: missing data.

### Running PySmooth

PySmooth takes the following arguments

- `-i` or `--input`: Name of the input genotype file. This MUST be provided
- '-o' or '--output': Prefix to name of output files to be generated. If not provided, default is `test`
- `-c` or `--chr`: list of chromosome names to perform analysis on. Names should be separated by comma (e.g `chr1,chr2,chr3`). Default is to run through all the chromosomes in the genotype file.
- `-l` or `--lower`: Lowest threshold for identifying singletons. Default is 0.70
- `-u`or `--upper`: Highest threshold for identifying singletons. Default is 0.98.
- `-g` or `--gap`: PySmooth iteratively identifies singletons starting with the highest threshold till the lowest threshold. This parameter is used to decreased the threshold at each iteration. Default is 0.02.
- `-k` : number of nearest neighbors to be used to assign correct genotype to singleton or missing item. Default value is 30.

Type the following command in the `terminal`, or `command prompt`, or  `anaconda command prompt` depending on your python installation or OS.

Shown below are two examples of running `PySmooth`.

`python run_smooth.py -i my_genotype_file.csv`

`python run_smooth.py -i my_genotype_file.csv -o my_output -c chr1,chr2,chr3 -l 0.80 -u 0.98 -g 0.02`

### Outputs

For each chromosome, PySmooth Generates the following outputs.

- `<output>_<chr>.stats.csv`, `<output>_<chr>.stats.csv`One statistics file that contains `%` of homozygous, heterozygous calls for each individual for the raw genoytpe file.
- One statistics file that contains `%` of homozygous, heterozygous calls for each individual after removing singletons from the raw genoytpe file. Singletons are marked missing.
- One summary statistics File that indicate how many singletons were detected. 
- One Heatmap image that displays genotype calls for the raw genoytpe file. 
- One Heatmap image that displays genotype calls for after removing singletons from the raw genoytpe file. 
- new Genotype File with singletons marked as missing. Singletons are marked as 'U'.

![alt text](https://github.com/lncRNAAddict/PySmooth/blob/main/example/crosspoints_in_chr1.heatmap.png)





