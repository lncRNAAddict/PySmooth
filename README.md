
# PySmooth

`PySmooth` is a python implementation of SMOOTH algorithm. Original SMOOTH was written in PASCAL with several limitations.

`PySmooth` contains one main functionality 

- `run_smooth.py`: Detect markers which are singletons in the genotype file and marks them as singletons. Correct the missing and singleton markers using k-nearest neighbor.


## Installation and Dependencies

PySmooth has been tested with Python 3.8.12 version. It should work with Python >= 3.0 version. We recommend installing the anaconda python distributon. Download anaconda python distribution from https://www.anaconda.com/products/distribution and install following the instructions provided.

PySmooth depends on the following python libraries. These libraries are already included in the anaconda distribution. Therefore, you do not need to install them.

- `numpy`
- `Pandas`
- `Sklearn`
- `matplotlib`

You can simply download the following scripts from `PySmooth` GitHub page and put them in a single folder. `PySmooth` can be executed by running the script `run_smooth.py`

- `utilities.py`
- `smooth.py`
- `ImputeMissingGenotype.py`
- `run_smooth.py`

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

A screeshot of a portion of an example input file is shown below

![Example Input Genotype File](https://github.com/lncRNAAddict/PySmooth/blob/main/example/GenotypeInput.PNG)

### Running PySmooth

PySmooth takes the following arguments

- `-i` or `--input`: Name of the input genotype file. This MUST be provided
- '-o' or '--output': Prefix to name of output files to be generated. If not provided, default is `test`
- `-c` or `--chr`: list of chromosome names to perform analysis on. Names should be separated by comma (e.g `chr1,chr2,chr3`). Default is to run through all the chromosomes in the genotype file.
- `-l` or `--lower`: Lowest threshold for identifying singletons. Default is 0.70
- `-u`or `--upper`: Highest threshold for identifying singletons. Default is 0.98.
- `-g` or `--gap`: PySmooth iteratively identifies singletons starting with the highest threshold till the lowest threshold. This parameter is used to decreased the threshold at each iteration. Default is 0.02.
- `-k` : number of nearest neighbors to be used to assign correct genotype to singleton or missing item. Default value is 30.

First, change working directory to the folder where the `PySmooth` scripts are stored. You can do that by simply typeing the following command in the `terminal`, or `command prompt`, or  `anaconda command prompt` depending on your python installation or OS.

`cd <path to where PySmooth scripts are stored`

Once the working directory is set, shown below are two examples of running `PySmooth`.

`python run_smooth.py -i <path to the genotype file>/my_genotype_file.csv`

The code above will analyze for each chromosome detected and generate all output files with prefix `test` in the folder `<path to the genotype file>`
  
`python run_smooth.py -i <path to the genotype file>/my_genotype_file.csv -o <path to output folder>/my_output -c chr1,chr2,chr3 -l 0.80 -u 0.98 -g 0.02`

The code above will analyze for chromosomes `chr1, chr2, and chr3`and generate all output files with prefix `my_output` in the folder `<path to output folder>`.


### Outputs

For each chromosome, PySmooth Generates the following outputs.

- Three csv files: `<output>_<chr>.stats.csv`, `<output>_<chr>_singletons_stats.csv`, and `<output>_<chr>_imputed_stats.csv` that contain `%` of homozygous, heterozygous calls for each individual for the raw genoytpe file, after singleton detection, and after error correction. Examples are shown below.

- Three bar plot png files: `<output>_<chr>.stats.png`, `<output>_<chr>_singletons_stats.png`, and `<output>_<chr>_imputed_stats.png` bar plot files that contains `%` of homozygous, heterozygous calls for each individual for the raw genoytpe file, after singleton detection, and after error correction, respectively. Example images are shown below.

![alt text](https://github.com/lncRNAAddict/PySmooth/blob/main/example/Slide3.PNG)

- Three heatmap files: `<output>_<chr>.heatmap.png`, `<output>_<chr>_singletons_heatmap.png`, and `<output>_<chr>_imputed_heatmap.png` that visualize a color-coded image of different genotype codes in the original file, after singleton detection, and after error correction, respectively. Example images are shown below.

![alt text](https://github.com/lncRNAAddict/PySmooth/blob/main/example/Slide2.PNG)

- `<output>_<chr>_singletons.csv`: genotype file with singleton detected. Singletons are marked as `S`. 
- `<output>_<chr>_imputed.csv`: genotype file after error correction.






