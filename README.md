# Introducing Ar(chitecture of)Du(plications)
**Please bear with me as this readme -like ArDu- is still a Work In Progress**

ArDu is a tool designed to screen target loci for genomic duplications, provide estimates of their copy number, and identify their genetic architecture (i.e. breakpoints, copy numbers, secondary rearrangements). As it is entirely written in Python, it should run on pretty much any machine and it should be modular and easy to modify to accommodate the users specific needs. 

## Installation and debugging. 
Simply download the main ArDu script and install it's dependencies. They are all listed in the `ArDu_environment.yml` file. Given that all dependencies are satistfied, ArDu can be used like any other python script. Here are some basic errors that you could encounter and how to fix them:
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



## Gene Copy Number Estimates, word of caution and best practises 

ArDu uses information from BAM alignment files to estimate a target loci copy number. The target's depth of coverage is normalised by a user-provided reference loci. Depending on the reference used, this normalised depth can be directly used a a copy number proxy (Claret et al. 2023). 
The choice of reference used for normalisation is critical to the quality of the copy number estimate. A wide range of genomic intervals can be used as reference, from whole chromosomes to a single housekeeping gene. However, we have seen an improvement in estimation precision by using exonic sequences of a housekeeping gene. WIP: An upcoming version will allow a less crude normalisation method through the use of targeted interval normalisation.
There are many copy number variant annotation tools that work on a genome-wide scale, but ArDu is specifically designed around a candidate locus approach, while you can use it on a large number of targets or even the entirety of an assembly annotation, it is not its intended use (expect long run times and little usability in the results).

## Options 

| Argument               | Description                                                                                                            | Default                     |
|------------------------|------------------------------------------------------------------------------------------------------------------------|--------------------------------------|
| **Mandatory Arguments**|
| `-b, --bam`           | A list of BAM files (one per line). Example format: <br> `~/bamfiles/sample1.bam` <br> `~/bamfiles/sample2.bam` <br> NOTE : all bam files must be indexed `samtools index bamfile.bam` and the indexes must be in the same directory as their respective bam files.  |                     |
  | `-r, --region`        | A 4-column tab delimited file. Regions with the same name will be pooled. Example format _chromosome  start  stop  name_, e.g. <br> `chr1  12000  12200  Target1` <br> `chr1  12400  12600  Target1` <br> `chr2  26000  26200  Reference` <br> `chr2  26400  26600  Reference` | |
| `-n, --norm`         | Name of the region to use for normalisation. This name **MUST** be present in the regions file, e.g. Reference in the example provided above.                                  |                               |
| `-o, --outfile`       | Prefix for output file names.                                                                                          |                      |
| **Plotting Arguments**|
| `--plot`              | Optional: produces a plot of the raw and normalised depth of coverage. Accepts the following extensions: png, jpeg, jpg, pdf, svg, eps. If used jointly with --breakpoint, putative breakpoints will be plotted. | `None`                               |
| `--plot-threshold`    | Threshold value to consider a duplication. Default value is 1.5, i.e. the expected normalised coverage of a heterozygous diploid organism possessing a single-copy and a duplicated copy of a given gene. | `1.5`                                |
| `--plot-interval`     | A tab-delimited file containing specific genomic intervals used for plotting, format: <br>`gene_name\tchromosome:start-stop`. | `None`                               |
| `--plot-proportion`   | Set the plotting interval to X times the size of the total gene span. Default = 2.                                      | `2`                                  |
| `--plot-slw`          | Sliding window size (bp) used for coverage representation, higher values result in a smoother depth of coverage.      | `1000`                               |
| `--plot-ylim`| Set Y-axis limits, format as two space separated numbers (e.g. `--plot-ylim -1 10`).                                    | `None`                               |
| `--plot-doclim`| Minimum and maximum normalised depth of coverage to use for plotting **and** breakpoint analysis (e.g. `--plot-doclim 0.5 5`). To keep minimum or maximum values, use "min" or "max" (e.g. `--plot-doclim min 5`, or `--plot-doclim 0.5 max`).              | `None`                               |
| `--plot-gene-pos`     | Shows target loci position on the plot.                                      | `False`                              |
| `--plot-force`        | Forces plotting regardless of depth of coverage value.                                                                  | `False`                              |
| `--plot-covar`        | Plots the covariance between mean depth of coverage and the variance instead of the mean depth of coverage.            | `False`                              |
| **Breakpoints Detection**  |
| `--bkp, --breakpoint`        | Optional: assess putative breakpoints. Two methods are available: 'ruptures' and 'rollingaverage'. <br> 'ruptures' uses the ruptures python package (see https://centre-borelli.github.io/ruptures-docs/). Associated options are --bkp-model, --bkp-pen or --bkp-nb. <br> 'rollingaverage' uses successive rolling average to detect shifts in depth of coverage. Associated options are --bkp-slw, --bkp-threshold, and --bkp-passes. <br> Setting this argument will create a '.breakpoints.tsv' output file containing the positions found by the rupture package, or the regions with depth of coverage shifting found with 'rollingaverage'. | `None`                               |
| `--bkp-slw`           | Window size for ruptures model and rolling average calculation. Default= 1 kb.                                          | `1000`                               |
| `--bkp-nb`            | Number of expected breakpoints in this structure, cannot be used alongside --bkp-pen.                                   | `None`                               |
| `--bkp-pen`           | Penalty parameter for ruptures algo.predict. A higher value will result in a higher penalty in breakpoint creation. Use this option if you don't have a prior on the breakpoint number in the structure you're looking at. | `None`                               |
| `--bkp-model`         | Model used by ruptures package. Default = 'l2'. Please refer to ruptures's documentation to see which options are best fitting to your specific case.                                                                         | `l2`                                 |
| `--bkp-algo`          | Algorithm used by ruptures package. Default = 'Window'.                                                                  | `Window`                             |
| `--bkp-threshold`     | Threshold for detecting shifts in depth of coverage.                                                                    | `1.0`                                |
| `--bkp-passes`        | Number of rolling average passes. Increasing will lessen the overall variation in depth of coverage.                   | `1`                                  |
| **Genotyping**        |                                                                                                                        |                                      |
| `--mutation`          | Returns the number of reads supporting each nucleotide at the given position(s) in a `.mutation.tsv` file. Takes as input a tab-delimited file formated as follows: <br>`chromosome   position  mutation_name(optional)`. | `None`                               |

## Command line examples
### Example run with bare minimum options:
```
python ardu.py -b bamlist.txt -r regions.txt -n Reference -o output_prefix
```

### Example run with plotting and breakpoints:
```
python ardu.py -b bamlist.txt -r regions.txt -n Reference -o output_prefix \
--plot png --plot_threshold 1.4 --breakpoint ruptures --bkp_model l2 --bkp_pen 10 
```

## Output files
### _output_prefix_ _coverage.tsv
Base ArDu (see _run with bare minimum option_ above) produces a _coverage.tsv file containing all the target loci depth of coverage statistics in lines, one column per bam file screened. Format is raw depth of coverage; mean depth of coverage; SD; median depth of coverage; normalised depth of coverage (ie copy number proxy).

### Plots
If --plot and either --plot-interval or --plot-proportion are set, a graphic representation of the duplicated loci will be produced in the requested format (png, jpeg, jpg, pdf, svg, eps). Target loci position can be plotted with --plot-gene-pos. If used conjointly with --breakpoint, the predicted breakpoints will be plotted and numbered in the same order as outputed in the _breakpoint.tsv file. 

### _output_prefix_ _breakpoints.tsv
If --breakpoint is set, ArDu will output a _breakpoints.tsv, containing either the predicted breakpoints position if ruptures was picked, numbered in the same order as the plot, or the genomic regions in which a signification depth of coverage shift was registered if rollingaverage was chosen. 

### _output_prefix_ _mutation.tsv
If --mutation is set, ArDu will output a _mutation.tsv file containign all nucleotides and total depth for the given position, format A=;T=;C=;G=;total_depth=. 

## Reference
If you use ArDu, please cite:
Claret Jean-Loup, Mestre Camille, Milesi Pascal and Labbé Pierrick, _in prep_.
**DOI:** [10.5281/zenodo.14922764](https://doi.org/10.5281/zenodo.14922764)

