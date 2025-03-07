# Introducing Ar(chitecture of)Du(plications)
**Please bear with us as this readme -like ArDu- is still under construction. WIP**

ArDu is a tool designed to screen target loci for genomic duplications, provide estimates of their copy number, and identify their genetic architecture (i.e. breakpoints, copy numbers, secondary rearrangements). As it is entirely written in Python, it should run on any machine and it is made to be modular and easy to modify by its users to accommodate their specific needs. ArDu depends on a few Python packages that can be found in the `.yml` file available under the main branch of the Git repository.

## Gene Copy Number Estimates

ArDu's main script uses the information from depth of coverage BAM alignment files to compute a proxy of gene copy number for target loci. It uses the depth of coverage of a reference genomic region to compute a normalized depth of coverage for the target loci that can be used as a direct proxy of its copy number (Claret et al. 2023). 

### Word of Caution and best practises 

The choice of the reference is paramount in the quality of the copy number estimate. From entire chromosomes to the exonic sequences of housekeeping genes, a wide array of genomic intervals can be used as a reference. However, we have seen an improvement in the estimation with the use of exonic sequences of a housekeeping gene.

While there are plenty of other copy number variants annotation tools that work on a genome-wide scale, ArDu is specifically designed to approach candidate loci.

**DOI:** [10.5281/zenodo.14922764](https://doi.org/10.5281/zenodo.14922764)

| Argument               | Description                                                                                                            | Default / Example                    |
|------------------------|------------------------------------------------------------------------------------------------------------------------|--------------------------------------|
| **Mandatory Arguments** |                                                                                                                        |                                      |
| `-b, --bam`           | A list of BAM files (one per line). Example format: <br> `~/bamfiles/sample1.bam` <br> `~/bamfiles/sample2.bam`   |                     |
| `-r, --region`        | A 4-column tab delimited file. Regions with the same name will be pooled. Example format: <br> `chr1\t12000\t12200\tTarget1` <br> `chromosome\tstart\tstop\tname` <br> `chr1\t12400\t12600\tTarget1` <br> `chr2\t26000\t26200\tReference` <br> `chr2\t26400\t26600\tReference` | |
| `-n, --norm`         | Name of the region to use for normalisation. This name MUST be present in the regions file.                                   |                               |
| `-o, --outfile`       | Prefix for output file names.                                                                                          |                      |
| **Plotting Arguments** |                                                                                                                        |                                      |
| `--plot`              | Optional: produces a plot of the raw and normalised depth of coverage. Accepts the following extensions: png, jpeg, jpg, pdf, svg, eps. If used jointly with --breakpoint, putative breakpoints will be plotted. | `None`                               |
| `--plot_threshold`    | Threshold value to consider a duplication. Default value is 1.5, i.e. the expected normalised coverage of a heterozygous diploid organism possessing a single-copy and a duplicated copy of a given gene. | `1.5`                                |
| `--plot_interval`     | A tab-delimited file containing specific genomic intervals used for plotting, format: <br>`gene_name\tchromosome:start-stop`. | `None`                               |
| `--plot_proportion`   | Set the plotting interval to X times the size of the total gene span. Default = 2.                                      | `2`                                  |
| `--plot_slw`          | Sliding window size (bp) used for coverage representation, higher values result in a smoother depth of coverage.      | `1000`                               |
| `--plot_min_norm_depth`| Minimum normalised depth of coverage to use for plotting.                                                             | `None`                               |
| `--plot_max_norm_depth`| Maximum normalised depth of coverage to use for plotting.                                                             | `None`                               |
| `--plot_gene_pos`     | When set, add vertical lines to indicate the duplicated gene position on the plot.                                      | `False`                              |
| `--plot_force`        | Forces plotting regardless of depth of coverage value.                                                                  | `False`                              |
| `--plot_covar`        | Plots the covariance between mean depth of coverage and the variance instead of the mean depth of coverage.            | `False`                              |
| **Breakpoints Detection** |                                                                                                                      |                                      |
| `--breakpoint`        | Optional: assess putative breakpoints. Two methods are available: 'ruptures' and 'rollingaverage'. <br> 'ruptures' uses the ruptures python package (see https://centre-borelli.github.io/ruptures-docs/). Associated options are --bkp_model, --bkp_pen or --bkp_nb. <br> 'rollingaverage' uses successive rolling average to detect shifts in depth of coverage. Associated options are --bkp_slw, --bkp_threshold, and --bkp_passes. <br> Setting this argument will create a '.breakpoints.tsv' output file containing the positions found by the rupture package, or the regions with depth of coverage shifting found with 'rollingaverage'. | `None`                               |
| `--bkp_slw`           | Window size for ruptures model and rolling average calculation. Default= 1 kb.                                          | `1000`                               |
| `--bkp_nb`            | Number of expected breakpoints in this structure, cannot be used alongside --bkp_pen.                                   | `None`                               |
| `--bkp_pen`           | Penalty parameter for ruptures algo.predict. A higher value will result in a higher penalty in breakpoint creation. Use this option if you don't have a prior on the breakpoint number in the structure you're looking at. | `None`                               |
| `--bkp_model`         | Model used by ruptures package. Default = 'l2'. Please refer to ruptures's documentation to see which options are best fitting to your specific case.                                                                         | `l2`                                 |
| `--bkp_algo`          | Algorithm used by ruptures package. Default = 'Window'.                                                                  | `Window`                             |
| `--bkp_threshold`     | Threshold for detecting shifts in depth of coverage.                                                                    | `1.0`                                |
| `--bkp_passes`        | Number of rolling average passes. Increasing will lessen the overall variation in depth of coverage.                   | `1`                                  |
| **Genotyping**        |                                                                                                                        |                                      |
| `--mutation`          | Optional: If set, this option will return the number of reads supporting each nucleotide at the given position(s) in a `.mutation.tsv` file. Takes as input a tab-delimited file formated as follows: <br>`chromosome\t position\tmutation_name(optional)`. | `None`                               |


## Example run with bare minimum options:
python ardu.py -b bamlist.txt -r regions.txt -n Reference -o output_prefix
## Example run with plotting and breakpoints:
python ardu.py -b bamlist.txt -r regions.txt -n Reference -o output_prefix \
  --plot png --plot_threshold 1.4 --breakpoint ruptures --bkp_model l2 --bkp_pen 10


