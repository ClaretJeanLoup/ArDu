# Introducing Ar(chitecture of)Du(plications)
**Please bear with me as this readme -like ArDu- is still a Work In Progress**

ArDu is a tool designed to screen target loci for genomic duplications, provide estimates of their copy number, and identify their genetic architecture (i.e. breakpoints, copy numbers, secondary rearrangements). As it is entirely written in Python, it should be modular and easy enough to modify to accommodate the users specific needs. 


## Table of contents

- [Introduction](#introducing-ararchitecture-ofduplications)
- [Quick start](#quick-start)
- [Output files](#output-files)
- [Installation](#installation)
- [Special note for Uppmax users](#special-note-for-uppmax-users)
- [Milesi's lab addendum](#milesis-lab-addendum)
- [Debugging](#debugging)
- [Gene copy number estimates, word of caution and best practices](#gene-copy-number-estimates-word-of-caution-and-best-practises)
- [Command line examples](#command-line-examples)
- [JunctionFinder](#junctionfinder)
- [Reference](#reference)

## Quick start 
### Example run with mandatory options:
```bash
python ardu.py -b bamlist.txt -r regions.txt -n Reference -o output_prefix
```

### Example run with plotting and breakpoints estimation:
```bash
python ardu.py -b bamlist.txt -r regions.txt -n Reference -o output_prefix \
--plot png --plot-threshold 1.4 --breakpoint ruptures --bkp-model l2 --bkp-pen 10 
```

## Output files
### _output_prefix_ coverage.tsv
The `coverage.tsv` file produced by ArDu stores multiple summary statistics per sample in a single semicolon-separated field :
| Flag        | Description |
|-------------|---------|
| `uDOC`    | Mean depth of coverage | 
| `sdDOC`   | Standard deviation |
| `medDOC`  | Median depth of coverage    |
| `CovBases`| Number of covered bases  | 
| `RCN`     | Relative copy number (normalised) 

To facilitate downstream analyses, a parser script is provided to extract specific fields into a clean tab-delimited file.

#### Usage

Extract selected fields from a coverage.tsv file:
```bash
python ArDu_Parser.py \
  -i output_prefix.coverage.tsv \
  -o parsed_output.tsv \
  --uDOC --RCN
```
At least one field flag must be provided (--uDOC, --sdDOC, --medDOC, --CovBases, --RCN). Missing or malformed entries are converted to NA.
The output is a wide-format matrix:

First column: locus
Each additional column corresponds to a sample-field combination:
sample1_uDOC   sample1_RCN   sample2_uDOC   sample2_RCN   ...

### Plots
If --plot and either --plot-interval or --plot-proportion are set, a graphic representation of the duplicated loci will be produced in the requested format (png, jpeg, jpg, pdf, svg, eps). Target loci position can be plotted with --plot-. If used conjointly with --breakpoint, the predicted breakpoints will be plotted and numbered in the same order as outputed in the _breakpoint.tsv file. 

### _output_prefix_ _breakpoints.tsv
If --breakpoint is set, ArDu will output a _breakpoints.tsv, containing either the predicted breakpoints position if ruptures was picked, numbered in the same order as the plot, or the genomic regions in which a signification depth of coverage shift was registered if rollingaverage was chosen. 

### _output_prefix_ _mutation.tsv
If --mutation is set, ArDu will output a _mutation.tsv file containign all nucleotides and total depth for the given position, format A=;T=;C=;G=;total_depth=. 


## Installation. 
Here's a quick step by step guide on how to install ArDu: <br>
1- Clone the github repository and create the environment with the `.yml` file. 
```
git clone https://github.com/ClaretJeanLoup/ArDu.git
cd ArDu
```
2- Create the environment
```
conda env create -f ArDu_environment.yml
```
3- Check if the environment was successfully created, it should appear under `ardu`.<br>
```
conda env list
```

4- Activate the environment and run the script
```
conda activate ardu
python /path/to/ardu/ArDu_version.py #replace path to match yours
```
Given that all dependencies are satistfied, ArDu can be used like any other python script.
<br>
<br>
### Special note for Uppmax users. 
For my Uppmax folks, here's how to properly set up ArDu:<br>
Navigate to your home directory and run part 1 and 2 of the installation process.
```
cd ~/
git clone https://github.com/ClaretJeanLoup/ArDu.git
cd ArDu
conda env create -f ArDu_environment.yml
```
Use the following command to add ArDu alisases to your `.bashrc` file. 
```
echo -e '\n# Load ArDu environment\nmlardu() {\n   conda deactivate\n    conda activate ardu\n    ml python/3.9.5/ pysam\n}\n# ArDu script alias\nalias ardu='"'"'python ~/ArDu/ArDu_1.0.py'"'"'' >> ~/.bashrc
```
Source your .bashrc file. 
```
source ~/.bashrc
```
Your ArDu environment can now be loaded with the command `mlardu`.
ArDu can then be run with the following command
```ardu```
>⚠️ Heads up: For reasons beyond my current Conda skills, a few packages might not get installed from the .yml file. If that happens, just activate the environment and install them manually with pip.
<br>

#### Milesi's lab addendum 
As of January 2026, a shared ArDu environment is available on the project:
```
module load Miniconda3
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate /crex/proj/snic2020-6-185/software/envs/ardu_shared/
export PYTHONNOUSERSITE=1



python /crex/proj/snic2020-6-185/software/envs/ardu_shared/ArDu_1.0.py/

```


>⚠️ You'll need to run the full source command line each time you connect to the cluster or launch a job via sbatch or interactive session. 
*Now go out there and hunt some duplications!* 
<br>
<br>
## Debugging
Here are some basic errors that you could encounter and how to fix them:
```
Processing TestRun
Processing genes:   0%|                                   | 0/2 [00:00<?, ?it/s]
Error processing BAM file TestRun.bam: no index available for pileup
```
--> No indexes found in the bam directory. Check for their presence and make sure their names match those of the bam. 

```
Processing TestRun
Processing genes:   0%|                                   | 0/3 [00:00<?, ?it/s]
Error processing BAM file TestRun.bam: invalid literal for int() with base 10: 'st'
```
--> Invalid line in the -r target file. Most likely the file contains a header that should be removed. Next version of ArDu will skip those lines. 

<br>
<br>

## Gene Copy Number Estimates, word of caution and best practises 
ArDu uses information from BAM alignment files to estimate a target loci copy number. The target's depth of coverage is normalised by a user-provided reference loci. Depending on the reference used, this normalised depth can be directly used a a copy number proxy (Claret et al. 2023). 
The choice of reference used for normalisation is critical to the quality of the copy number estimate. A wide range of genomic intervals can be used as reference, from whole chromosomes to a single gene. However, we have seen an improvement in estimation precision by using exonic sequences of a few housekeeping genes. 
TArDu is specifically designed around a candidate locus approach, while you can use it on a large number (i.e. thousands) of targets or even the entirety of an assembly annotation, it is not its intended use (expect long run times and little usability in the results).
<br>
| Argument                 | Description                                                                                                                                                                                                                                                                                                               | Default    |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------- |
| **Mandatory Arguments**  |                                                                                                                                                                                                                                                                                                                           |            |
| `-b, --bam`              | A list of BAM files (one per line). Example format: <br> `~/bamfiles/sample1.bam` <br> `~/bamfiles/sample2.bam` <br> NOTE: all BAM files must be indexed (`samtools index bamfile.bam`) and indexes must be in the same directory as their respective BAM files.                                                          |            |
| `-r, --region`           | A 4-column tab-delimited file. Regions with the same name will be pooled. Format: `chromosome start stop name`, e.g. <br> `chr1 12000 12200 Target1` <br> `chr1 12400 12600 Target1` <br> `chr2 26000 26200 Reference` <br> `chr2 26400 26600 Reference`                                                                  |            |
| `-n, --norm`             | Name of the region used for normalisation. This name **must** be present in the region file.                                                                                                                                                                                                                              |            |
| `-o, --outfile`          | Prefix for output file names.                                                                                                                                                                                                                                                                                             |            |
| **Plotting Arguments**   |                                                                                                                                                                                                                                                                                                                           |            |
| `--plot`                 | Produces a plot of raw and normalised depth of coverage. Accepts: png, jpeg, jpg, pdf, svg, eps. If used with `--breakpoint`, putative breakpoints will be plotted.                                                                                                                                                       | `None`     |
| `--plot-threshold`       | Threshold to consider a duplication. Default reflects expected normalised coverage of a heterozygous diploid organism with one duplicated copy.                                                                                                                                                                           | `1.4`      |
| `--plot-interval`        | Tab-delimited file defining plotting intervals: `gene_name\tchromosome:start-stop`.                                                                                                                                                                                                                                       | `None`     |
| `--plot-proportion`      | Sets plotting window to X times the total gene span.                                                                                                                                                                                                                                                                      | `2`        |
| `--plot-auto`            | Automatically expands plot interval by probing coverage. Probes are extended outward until a drop in coverage is detected.                                                                                                                                                                                                | `False`    |
| `--probe-size`           | Probe length (bp) used during interval extension.                                                                                                                                                                                                                                                                         | `500`      |
| `--probe-threshold`      | Coverage ratio threshold of probes relative to target.                                                                                                                                                                                                                                                                    | `0.8`      |
| `--probe-number`         | Number of probes per extension round.                                                                                                                                                                                                                                                                                     | `20`       |
| `--probe-spacing`        | Distance between probe starts (bp).                                                                                                                                                                                                                                                                                       | `1000`     |
| `--probe-drops`          | Number of consecutive low probes required to stop extension.                                                                                                                                                                                                                                                              | `10`       |
| `--probe-use-median`     | Use median probe coverage instead of mean (useful in noisy regions).                                                                                                                                                                                                                                                      | `False`    |
| `--max-extension`        | Maximum interval extension on each side of the target (bp).                                                                                                                                                                                                                                                               | `5000000`  |
| `--plot-slw`             | Sliding window size (bp) for coverage representation. Larger values produce smoother profiles.                                                                                                                                                                                                                            | `1000`     |
| `--plot-ylim`            | Y-axis limits (`min max`).                                                                                                                                                                                                                                                                                                | `None`     |
| `--plot-doclim`          | Limits for normalised depth of coverage used in plotting and breakpoint analysis. Supports `min`/`max`.                                                                                                                                                                                                                   | `None`     |
| `--plot-target`          | Highlights target region on the plot.                                                                                                                                                                                                                                                                                     | `False`    |
| `--plot-force`           | Forces plotting even if coverage thresholds are not met.                                                                                                                                                                                                                                                                  | `False`    |
| `--plot-covar`           | Plots covariance between mean coverage and variance instead of mean coverage alone.                                                                                                                                                                                                                                       | `False`    |
| **Breakpoint Detection** |                                                                                                                                                                                                                                                                                                                           |            |
| `--bkp, --breakpoint`    | Enables breakpoint detection using either `ruptures` or `rollingaverage`. <br> `ruptures`: uses Ruptures library ([https://centre-borelli.github.io/ruptures-docs/](https://centre-borelli.github.io/ruptures-docs/)). <br> `rollingaverage`: detects shifts using rolling means. <br> Outputs a `.breakpoints.tsv` file. | `None`     |
| `--bkp-slw`              | Window size for breakpoint detection (bp).                                                                                                                                                                                                                                                                                | `1000`     |
| `--bkp-nb`               | Expected number of breakpoints (cannot be used with `--bkp-pen`).                                                                                                                                                                                                                                                         | `2`        |
| `--bkp-pen`              | Penalty parameter for rupture detection. Higher values produce fewer breakpoints.                                                                                                                                                                                                                                         | `None`     |
| `--bkp-model`            | Model used by Ruptures.                                                                                                                                                                                                                                                                                                   | `l2`       |
| `--bkp-algo`             | Algorithm used by Ruptures.                                                                                                                                                                                                                                                                                               | `BottomUp` |
| `--bkp-threshold`        | Threshold for detecting coverage shifts.                                                                                                                                                                                                                                                                                  | `0.5`      |
| `--bkp-passes`           | Number of smoothing passes in rolling average mode.                                                                                                                                                                                                                                                                       | `1`        |
| **Genotyping**           |                                                                                                                                                                                                                                                                                                                           |            |
| `--mutation`             | Counts nucleotide support at given positions from BAM files. Input format: `chromosome position mutation_name(optional)`. Outputs `.mutation.tsv`.                                                                                                                                                                        | `None`     |


## Command line examples
### Example run with bare minimum options:
```
python ardu.py -b bamlist.txt -r regions.txt -n Reference -o output_prefix
```

### Example run with plotting and breakpoints:
```
python ardu.py -b bamlist.txt -r regions.txt -n Reference -o output_prefix \
--plot png --plot-threshold 1.4 --breakpoint ruptures --bkp-model l2 --bkp-pen 10 
```

## Output files
### _output_prefix_ _coverage.tsv
Base ArDu (see _run with bare minimum option_ above) produces a _coverage.tsv file containing all the target loci depth of coverage statistics in lines, one column per bam file screened. Format is raw depth of coverage; mean depth of coverage; SD; median depth of coverage; normalised depth of coverage (ie copy number proxy).

### Plots
If --plot and either --plot-interval or --plot-proportion are set, a graphic representation of the duplicated loci will be produced in the requested format (png, jpeg, jpg, pdf, svg, eps). Target loci position can be plotted with --plot-target. If used conjointly with --breakpoint, the predicted breakpoints will be plotted and numbered in the same order as outputed in the _breakpoint.tsv file. 

### _output_prefix_ _breakpoints.tsv
If --breakpoint is set, ArDu will output a _breakpoints.tsv, containing either the predicted breakpoints position if ruptures was picked, numbered in the same order as the plot, or the genomic regions in which a signification depth of coverage shift was registered if rollingaverage was chosen. 

### _output_prefix_ _mutation.tsv
If --mutation is set, ArDu will output a _mutation.tsv file containign all nucleotides and total depth for the given position, format A=;T=;C=;G=;total_depth=. 

## JunctionFinder
We developped a simple module to check for the existence of junction sequences (i.e. sequences overlapping the breakpoints of duplications events) in short reads. It uses the .bkps.tsv file generqted by ArDu's main script as input, alongside the targets bam file. After detecting softcliped sequences around each putative breakpoints, aligned and unmapped reads stored in the targets bamfile are screened for sequences containing softclips originating from opposing breakpoints using an Aho-Corasick automaton. This allows to extract potential junction sequences overlapping the structural variants breakpoints.  

Breakpoints + BAM files
          │
          ▼
Extract soft-clipped reads
          │
          ▼
Soft-clip FASTA + metadata TSV
          │
          ├── Optional: BLAST alignment
          │        │
          │        ▼
          │   Filter hits near breakpoints
          │
          ▼
Build soft-clip automaton
          │
          ▼
Scan reads for multiple clips
          │
          ▼
Junction reads 

| Argument            | Description                                     |
| ------------------- | ----------------------------------------------- |
| `-b --bam_list`     | File listing BAM paths                          |
| `-i --input`        | Breakpoint file from ArDu main script                                 |
| `-o --output`       | Output prefix                                   |
| `-s --size`         | Minimum soft-clip length (default: 30)          |
| `-e --extension`    | Region around breakpoint (default: 30 bp)       |
| `--blast`           | Reference genome for BLAST                      |
| `--pairs`           | Allowed breakpoints                     |
| `--junction`        | Enable junction read detection                  |
| `--junction-region` | Restrict junction search region (chr:start-end) |


## Reference
If you use ArDu, please cite:
Claret Jean-Loup, Mestre Camille, Milesi Pascal and Labbé Pierrick, _in prep_.
**DOI:** [10.5281/zenodo.14922764](https://doi.org/10.5281/zenodo.14922764)

