# Loading python built-in modules 
import argparse
import csv
import math
import os
import statistics 
import subprocess
import time

# Loading specialised modules 
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import pysam
import ruptures as rpt
from tqdm import tqdm

def add_plot_params_box(fig, args):
    """Add a summary of plot parameters outside the plotting area."""
    lines = []

    if args.plot_threshold is not None:
        lines.append(f"Threshold: {args.plot_threshold}")
    if args.plot_slw is not None:
        lines.append(f"Sliding Window: {args.plot_slw} bp")
    if args.plot_bin is not None:
        lines.append(f"Bin: {args.plot_bin} bp")
    if args.plot_ylim:
        lines.append(f"Y-limits: {args.plot_ylim[0]}–{args.plot_ylim[1]}")
    if args.plot_doclim:
        lines.append(f"DoC limits: {args.plot_doclim[0]}–{args.plot_doclim[1]}")
    if args.plot_force:
        lines.append("Force Plotting: ON")

    param_text = "\n".join(lines)

    fig.text(
        0.05, 0.025, param_text, #x,y format, 0.5 centers, 0.05 plots toward the bottom because figure axis goes from 0 to 1
        fontsize=8,
        va='top', ha='left',
        bbox=dict(boxstyle="square",facecolor="white", edgecolor="gray", alpha=0.9)
    )


# Functions for slw smoothing
def window_average(arr: np.ndarray, w: int):
    numbers_series = pd.Series(arr)
    windows = numbers_series.rolling(w, center=True, min_periods=1)  # Use center=True to align the window 
    moving_averages = windows.median()
    return moving_averages

def window_variance(arr: np.ndarray, w: int):
    numbers_series = pd.Series(arr)
    windows = numbers_series.rolling(w, center=True, min_periods=1)  # Use center=True to align the window
    moving_variances = windows.var()
    return moving_variances

def calculate_coverage_stats(regions, bam_file, quantile):
    """Compute mean, median, and standard deviation from depth of coverage for each gene."""
    data_list = []
    try:
        total_genes = len(regions)
        
        with pysam.AlignmentFile(bam_file, "rb") as samfile:
            for name, region_list in tqdm(regions.items(), total=total_genes, desc="Processing genes"):
                coverage_values = []

                for region in region_list:
                    if not region:
                        continue

                    # Collect all coverage values for the region
                    for pileupcolumn in samfile.pileup(region=region, min_base_quality=13, min_mapping_quality=20):
                        coverage_values.append(pileupcolumn.n)
                        

                # Remove the smallest coverage values
                #if coverage_values and quantile:
                #    quantile = np.quantile(coverage_values, args.quantile)
                #    coverage_values = [val for val in coverage_values if val >= quantile]

                # Compute statistics if we have coverage values
                if coverage_values:
                    mean_coverage = np.mean(coverage_values)
                    median_coverage = np.median(coverage_values)
                    sd_coverage = np.std(coverage_values, ddof=1) if len(coverage_values) > 1 else 0.0
                    covbases = len(coverage_values)
                else:
                    mean_coverage, median_coverage, sd_coverage, covbases = 'NA', 'NA', 'NA', 'NA'

                # Append results
                data_list.append({
                    "gene": name,
                    "mean": mean_coverage,
                    "median": median_coverage,
                    "sd": sd_coverage,
                    "coveredbases": covbases
                })

    except Exception as e:
        print(f"Error processing BAM file {bam_file}: {e}")
        return pd.DataFrame()

    return pd.DataFrame(data_list)

# Function to detect significant shifts between normalised depth and moving average
def detect_shifts(data: pd.DataFrame, threshold: float, window_size: int, passes: int):
    norm_smooth = data['norm'].values.copy()
    moving_average_smooth = data['moving_average'].values.copy()
    
    for _ in range(passes):
        norm_smooth = window_average(norm_smooth, window_size)
        moving_average_smooth = window_average(moving_average_smooth, window_size)
    
    diff = moving_average_smooth - norm_smooth
    significant_shifts = data[np.abs(diff) > threshold]
    return significant_shifts

# Function to calculate the nucleotide counts and total depth of coverage
def calculate_nucleotide_counts(bam_file, chromosome, position):
    position = int(position)  # pos integer check
    samfile = pysam.AlignmentFile(bam_file, "rb")
    counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    total_depth = 0
    for pileupcolumn in samfile.pileup(region=f"{chromosome}:{position}-{position}", stepper='all'):
        if pileupcolumn.pos == position - 1:  # 0based indexing
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    base = pileupread.alignment.query_sequence[pileupread.query_position]
                    counts[base] += 1
                    total_depth += 1
    samfile.close()
    return counts, total_depth

# Function to get regions from a bed file
def get_regions(regions):
    regions_dict = {}
    with open(regions, "r") as f:
        for line_str in f:
            if not line_str.strip() or line_str.startswith("#"):  # Skip empty lines
                continue
            line = line_str.strip().split('\t')
            if len(line) != 4:
                print(f"Ignoring invalid line: {line}")
                continue
            if line[3] in regions_dict.keys():
                regions_dict[line[3]].append(f"{line[0]}:{line[1]}-{line[2]}") # {name: ["chromosome:start-stop", "chromosome:start-stop"]}
            else:
                regions_dict[line[3]] = [f"{line[0]}:{line[1]}-{line[2]}"]
    return regions_dict

# Function to compute the total span of a candidate gene
def get_totalspan(regions, gene_name):
    gene_intervals = []
    
    # Read bed file
    with open(regions, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) != 4:
                continue  # malformed lines
            
            chromosome = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            gene = parts[3]
            
            # Check if the gene matches the specified gene name
            if gene == gene_name:
                gene_intervals.append((start, end, chromosome))
                chrm = chromosome
    
    # If gene_intervals gene was not found in the file
    if not gene_intervals:
        return None
    
    # Calculate the minimum start and max end for the loci intervals
    min_start = min(interval[0] for interval in gene_intervals)
    max_end = max(interval[1] for interval in gene_intervals)
    
    # Return the interval as tuple format chromosome, min_start, max_end
    totalspan = (chrm, min_start, max_end)
    return totalspan

# Parse tab-delimited file for plot intervals
def parse_plot_intervals(file_path):
    plot_intervals = {}
    with open(file_path, 'r') as f:
        for line in f:
            gene, interval = line.strip().split('\t')
            plot_intervals[gene] = interval
    return plot_intervals

def compute_sliding_covariance(data, window_size):
    """Compute covariance of mean depth and variance in a sliding window."""
    mean_series = data['depth'].rolling(window_size, center=True).mean()
    var_series = data['depth'].rolling(window_size, center=True).var()
    cov_series = mean_series.rolling(window_size, center=True).cov(var_series)
    
    # Scale covariance to match normalised depth range
    if not cov_series.isna().all():  
        cov_min, cov_max = cov_series.min(), cov_series.max()
        norm_min, norm_max = data['norm'].min(), data['norm'].max()
        cov_series = (cov_series - cov_min) / (cov_max - cov_min) * (norm_max - norm_min) + norm_min
        # Center around 1
        cov_series = cov_series - cov_series.mean() + 1
        cov_series[cov_series < 1] = float('nan')

    return data['pos'], cov_series


# Main function
def main():
    # Mandatory arguments
    parser = argparse.ArgumentParser(description="ArDu -Architecture of Duplications- Identify duplications structure via depth of coverage analysis.")
    parser.add_argument("-b", "--bam", help="Input a list of BAM files (one per line)")
    parser.add_argument("-r", "--region", required=True, help="Takes as input a 4 columns file (tab delimited, as bed format: 'chromosome\tstart\tstop\tgene_name') containing regions of \
                        interest. Features like genes composed of multiple regions (e.g. exons) will be pooled together if they share the same name (in the 4th column)")
    parser.add_argument("-n","--norm", required=True, help="Name of the region to use for normalisation describe in the regions file.")
    parser.add_argument("-o", "--outfile", required=True, help="Prefix for output file names")

    #Optional arguments
    parser.add_argument("-q","--quantile", type=float, help="Set the quantile of depth of coverage that is excluded from all analyses.DEPRECATED") #deprecated

    ## Plotting
    parser.add_argument("--plot", 
        type=str, 
        default=None, 
        help="Optional: produces a plot of the raw and normalised depth of coverage. Accepts as value the following extensions for the plot file: \
        png, jpeg, jpg, pdf, svg, eps. If used jointly with --breakpoint, the putative breakpoints will be plotted.")
    parser.add_argument("--plot_threshold", 
        type=float, 
        default= 1.4, 
        help="Treshold value to consider a duplication. Default value is 1.4, i.e. the expected normalised coverage of a heterozygous diploid organism possessing a mono-copy and a duplicated copy of a given gene.")
    parser.add_argument("--plot_interval", 
        help="A tab delimited file containing specific genomic interval used for plotting, forma: 'gene_name\tchromosome:start-stop'")
    parser.add_argument("--plot_proportion", 
        type=float, 
        default= 2, 
        help="Set the plotting interval to X time the size of the total gene span. Default = 2.")
    parser.add_argument("--plot_slw", 
        type=int, 
        default= 1000, 
        help="Sliding window size (pb) used for coverage representation, higher values will result in a smoother depth of coverage representation.")
    parser.add_argument("--plot_gene_pos",
        action='store_true', 
        help="When set, add vertical lines to indicate the duplicated gene position on the plot.")
    parser.add_argument("--plot_force",
        action='store_true', help="Forces plotting regardless of depth of coverage value.")
    parser.add_argument("--plot_covar",
        action='store_true',
        help="Plots the covariance between mean depth of coverage and the variance instead of the mean depth of coverage.")
    parser.add_argument('--plot_ylim',
        type=float,
        nargs=2,
        help='Set Y-axis limits as two separate numbers (e.g. a b)')
    parser.add_argument(
        "--plot_doclim",
        type=str,
        nargs=2,
        help="Minimum and maximum normalised depth of coverage to use for plotting and bkp analysis. "
            "Use 'min' or 'max' to refer to the actual min/max of the dataset, e.g. 'min 0.5' or '0.5 max'.")
    parser.add_argument("--plot_bin",
        type=int,
        help="Bins of X bp are added on top of sliding windows.")
    parser.add_argument("--plot_dpi",
        type=int,
        default= 300,
        help="Set the plot dpi, default= 300.")
    parser.add_argument("--plot_param",
        action='store_true',
        help="Prints the parameters on the plot.")

    ## Breakpoints detection
    parser.add_argument("--breakpoint", choices=["ruptures", "rollingaverage"], help="Optional: assess putative breakpoints. Two methods are available: 'ruptures' and 'rollingaverage'. \n\
                        'ruptures' uses the ruptures python package (see https://centre-borelli.github.io/ruptures-docs/). Associated options are --bkp_model, --bkp_pen or --bkp_nb. \n\
                        'rollingaverage' uses successive rolling average to detect shifts in depth of coverage. Associated options are --bkp_slw, --bkp_threshold and --bkp_passes. \n\
                        Setting this argument will create a '.breakpoints.tsv' output file containing the positions found by the rupture package, or the regions of depth of coverage \
                        shifts found with 'rollingaverage'.")
    parser.add_argument("--bkp_slw", type=int, default=1000, help="Window size for ruptures model and rolling average calculation. Default= 1 kb.")
    parser.add_argument("--bkp_signal", type=str, help="Signal on which to run ruptures, possibilities are norm for normalised depth of coverage, variance or mean for the per sliding window variance and mean or covar for their covariance. \n\
                        Default to mean.")

    parser.add_argument("--bkp_nb", type=int, default=2, help="Number of expected breakpoints in this structure, cannot be used alongside --bkp_pen. Default = 2")
    parser.add_argument("--bkp_pen", type=int, help="Penalty parameter for ruptures algo.predict. A higher value will result in a higher penalty in breakpoint creation.\n\
                        use this option if you don't have a prior on the bkp number in the structure you're looking at.")
    parser.add_argument("--bkp_model", type=str, default="l2", help="Model used by ruptures package. Default ='l2'.")
    parser.add_argument("--bkp_algo", type=str, default="BottomUp", help="Algorithm used by ruptures package. Default ='BottomUp'.")
    parser.add_argument("--bkp_threshold", type=float, default=0.5, help="Threshold for detecting shifts in depth of coverage.")
    parser.add_argument("--bkp_passes", type=int, default=1, help="Number of rolling average passes. Increasing will lessen the overal variation in depth of coverage.")

    ## Diagnostic mutation genotyping 
    parser.add_argument("--mutation", help="Optional: if set, this option will return the number of reads supporting each nucleotides at the given position(s) in a .mutation.tsv file. Takes as input a tab-delimited file\
                         with the following columns: chromosome, position and an optional third column containing the name of the muations you're screening for (tab-delimited).")

    args = parser.parse_args()

    start_time = time.time()  # Start
    date_str = time.strftime("%Y-%m-%d")
    ArDuBase = os.path.splitext(os.path.basename(args.bam))[0]
    ArDuBaseRun = f"{ArDuBase}-{date_str}"
    output_dir = f"ArDuRun-{ArDuBaseRun}/"
    os.makedirs(output_dir, exist_ok=True)

    # Check if the provided bam is a valid file path
    with open(args.bam, 'r') as f:
        bam_files = f.read().splitlines()

    # Create dictionnary for the regions to screen and check if normalisation region is present in the regions file
    regions_dict = get_regions(args.region)
    if args.norm not in regions_dict.keys():
        raise ValueError(f"Normalisation region '{args.norm}' not found in the regions file.")

    # per-bam coverage dataf
    coverage_dfs = []
    breakpoints_data = []
    bam_names_without_ext = [os.path.basename(bam_file).removesuffix('.bam') for bam_file in bam_files] #for coverage column names 

    # Process each bam file
    for bam_file in bam_files:
        try:
            # Extract base name of the bam file
            bam_name = os.path.basename(bam_file).removesuffix('.bam')
            print(f"Processing {bam_name}")

            # Compute coverage values
            gene_coverage_df = calculate_coverage_stats(regions_dict, bam_file, args.quantile)

            # Check if the normalisation gene is present
            if args.norm not in gene_coverage_df['gene'].values:
                print(f"Normalisation gene {args.norm} not found in coverage for {bam_file}. Skipping.")
                continue

            # median coverage for the normalisation gene
            ref_median = gene_coverage_df.loc[gene_coverage_df['gene'] == args.norm, 'mean'].values[0]
            if not isinstance(ref_median, (int, float)):
                print(f"The depth of coverage of the reference used for normalisation is invalid in {bam_file}. Skipping.")
                continue

            # Normalisation
            gene_coverage_df['normalised'] = gene_coverage_df['mean'].apply(lambda x: x / ref_median if isinstance(x, (int, float)) else "NA")
            gene_coverage_df['bam_file'] = bam_name
            coverage_dfs.append(gene_coverage_df)

            # -----------------------------
            # Breakpoints and plotting
            if not args.plot:
                continue
            else:
                if args.plot_force:
                    filtered_genes_df = gene_coverage_df[gene_coverage_df['gene'] != args.norm]
                else:
                    filtered_genes_df = gene_coverage_df[gene_coverage_df['normalised'] > args.plot_threshold]
                # Process filtered gene
                for idx, row in filtered_genes_df.iterrows():
                    gene = row['gene']
                    mean_coverage = row['mean']
                    sd_coverage = row['sd']
                    normalised_coverage = row['normalised']
                    # depth file creaton
                    if args.plot_interval and args.plot:
                        plot_intervals = parse_plot_intervals(args.plot_interval)
                        str_incr_span = plot_intervals[gene]
                        samfile = pysam.AlignmentFile(bam_file, "rb")
                        depth_data = [(p.pos, p.n) for p in samfile.pileup(region=str_incr_span, min_base_quality=13, min_mapping_quality=20)]
                    elif args.plot_proportion and args.plot:
                        totalspan = get_totalspan(args.region, gene)
                        if totalspan is None:
                            print(f"Warning: No intervals found for gene {gene} in BED file.")
                            continue
                        chromosome, plotStart, plotStop = totalspan
                        gene_length = plotStop - plotStart
                        extension = gene_length * args.plot_proportion
                        plotStart -= int(extension)
                        plotStop += int(extension)
                        str_incr_span = f"{chromosome}:{plotStart}-{plotStop}"
                        samfile = pysam.AlignmentFile(bam_file, "rb")
                        depth_data = [(p.pos, p.n) for p in samfile.pileup(region=str_incr_span, min_base_quality=13, min_mapping_quality=20)] #stepper='all'
                    else:
                        continue

                    samfile.close()
                    d = pd.DataFrame(depth_data, columns=["pos", "depth"])
                    d['depth'] = pd.to_numeric(d['depth'], errors='coerce')
                    d["norm"] = d.depth / ref_median
                    if args.plot_doclim: #limits for plotted DoC
                        def parse_bound(val, is_min):
                            """Used for --plot_doclim, parses the inputed values to check for min or max."""
                            if val.lower() == "min":
                                return normalised_min
                                print(normalised_min)
                            elif val.lower() == "max":
                                return normalised_max
                                print(normalised_max)
                            else:
                                return float(val) 
                        normalised_min = d["norm"].min()
                        normalised_max = d["norm"].max()
                        a = parse_bound(args.plot_doclim[0], is_min=True)
                        b = parse_bound(args.plot_doclim[1], is_min=False)

                        d = d[(d["norm"] >= a) & (d["norm"] <= b)].reset_index(drop=True)
                    d["moving_average"] = window_average(d["norm"].to_numpy(), args.plot_slw)
                    d["moving_variances"] = window_variance(d["norm"].to_numpy(), args.plot_slw)
                    d['pos'], d['covariance'] = compute_sliding_covariance(d, args.plot_slw)
                    
                    if args.plot_bin: # Bin the df
                        bin_size = args.plot_bin
                        bin_edges = np.arange(d['pos'].min(), d['pos'].max() + bin_size, bin_size)
                        d['genomic_bins'] = pd.cut(d['pos'], bins=bin_edges, include_lowest=True)
                        for col in ['pos', 'moving_average', 'moving_variances', 'covariance']:
                            d[col] = d.groupby('genomic_bins', observed=True)[col].transform('mean')

                    # Plotting
                    fig, ax = plt.subplots(figsize=(11, 7))
                    if args.plot_ylim:
                        ax.set_ylim(args.plot_ylim)

                    ax.grid(axis='y', linestyle='dashdot', linewidth=0.5)
                    ax.tick_params(axis='both', labelsize=12)
                    ax.plot(d.pos, d.norm, label="Normalised Depth of Coverage", color="lightsteelblue", linestyle='none', marker='.')

                    if args.plot_covar:
                        ax.plot(d.pos, d.covariance, label="Sliding Window Covariance", color="blue", linewidth=2)
                        ax.plot(d.pos, d.moving_variances, label="Sliding Window variance", color="red", alpha=0.5, linewidth=2)
                        ax.plot(d.pos, d.moving_average, label="Smoothed Depth of Coverage", color="black", alpha=0.5, linewidth=2)
                    else:
                        ax.plot(d.pos, d.moving_average, label="Smoothed Depth of Coverage", color="black", linewidth=2)

                    # Breakpoint plotting
                    if args.breakpoint:
                        print(f"Breakpoint method: {args.breakpoint}")
                        if args.breakpoint == "ruptures":
                            choose_algo = args.bkp_algo
                            algo_class = getattr(rpt, choose_algo)
                            if args.bkp_signal:
                                signal = d[args.bkp_signal]
                            else:
                                signal = d["moving_average"]
                            algo = algo_class(model=args.bkp_model).fit(signal.to_numpy().reshape(-1, 1))
                            result = algo.predict(pen=args.bkp_pen) if args.bkp_pen else algo.predict(n_bkps=args.bkp_nb)
                            if result:
                                for line_num, b in enumerate(result[:-1], start=1):
                                    bp_pos = int(d.loc[b, 'pos'])
                                    breakpoints_data.append({
                                        'bam_file': bam_name,
                                        'gene': gene,
                                        'chromosome': str_incr_span.split(':')[0],
                                        'position': bp_pos,
                                        'breakpoint_number': line_num,
                                        'method': 'ruptures'
                                    })
                        elif args.breakpoint == "rollingaverage":
                            significant_shifts = detect_shifts(d, args.bkp_threshold, args.bkp_slw, args.bkp_passes)
                            regions = []
                            consecutive_positions = []
                            for idx, row in significant_shifts.iterrows():
                                if not consecutive_positions or row['pos'] - consecutive_positions[-1] <= 1:
                                    consecutive_positions.append(row['pos'])
                                else:
                                    regions.append((consecutive_positions[0], consecutive_positions[-1]))
                                    consecutive_positions = [row['pos']]
                            if consecutive_positions:
                                regions.append((consecutive_positions[0], consecutive_positions[-1]))

                    # Draw breakpoints
                    if args.breakpoint == "ruptures":
                        for i, b in enumerate(result[:-1]):
                            ax.axvline(d.loc[b]["pos"], color="darkred", linestyle="-", linewidth=1.5)
                            ax.text(d.loc[b]["pos"], d.loc[b]["norm"], f"{i+1}", color="darkred", fontsize=10, weight='bold', va='bottom')
                    elif args.breakpoint == "rollingaverage":
                        for i, region in enumerate(regions):
                            region_center = (region[0] + region[1]) / 2
                            region_data = d[(d['pos'] >= region[0]) & (d['pos'] <= region[1])]
                            if not region_data.empty:
                                max_norm_depth_index = region_data['norm'].idxmax()
                                max_norm_depth = region_data.loc[max_norm_depth_index, 'norm']
                                ax.axvspan(region[0], region[1], color='darkred')
                                ax.text(region_center, max_norm_depth + 0.1, f"{i+1}", color="darkred", fontsize=10, weight='bold')

                    if args.plot_gene_pos:
                        totalspan = get_totalspan(args.region, gene)
                        mask = (d.pos >= totalspan[1]) & (d.pos <= totalspan[2])
                        ax.plot(d.pos[mask], d.moving_average[mask], color="mediumpurple", label=f"{gene}", zorder=5, linewidth=2.5)

                    fig.suptitle(f"{gene} locus in {bam_name}", fontsize=14, y=1.02)
                    ax.set_xlabel(f"Genomic position: {str_incr_span}", fontsize=12)
                    ax.set_ylabel("Normalised Depth of Coverage", fontsize=12)

                    if args.breakpoint:
                        custom_line = Line2D([0], [0], color='darkred', linestyle='-', linewidth=1.5, label='Predicted Breakpoints')
                        handles, labels = ax.get_legend_handles_labels()
                        handles.append(custom_line)
                        labels.append('Predicted Breakpoints')
                        ax.legend(handles, labels, loc='upper right')
                    else:
                        ax.legend(loc='upper right')

                    # Add parameter box
                    if args.plot_param:
                        plt.subplots_adjust(bottom=0.25)
                        add_plot_params_box(fig, args)

                    plt.tight_layout()
                    output_dirgene = f"{output_dir}/{gene}-plots/"
                    os.makedirs(output_dirgene, exist_ok=True)
                    plt.savefig(os.path.join(output_dirgene, f"{bam_name}_{gene}_plot.{args.plot}"),
                                bbox_inches='tight', dpi=args.plot_dpi)
                    plt.close()

        except Exception as e:
            print(f"{bam_file}: {e}")

    # -----------------------------
    # Combine the bam dataf and write the final tsv file

    if coverage_dfs:
        # Concat
        combined_df = pd.concat(coverage_dfs, ignore_index=True)
        # formats the values 
        combined_df['values'] = combined_df.apply(
            lambda row: f"{round(row['mean'], 2)};{round(row['sd'], 2)};{round(row['median'], 2)};{row['coveredbases']};{round(row['normalised'], 1)}" 
                        if isinstance(row['normalised'], (int, float)) else "NA;NA;NA;NA;NA",
            axis=1
        )
        # pivot the combined df so each row is a gene and each column apprt from gene is a bam
        pivot_df = combined_df.pivot(index='gene', columns='bam_file', values='values').reset_index()
        cols = ['gene'] + [bam for bam in bam_names_without_ext if bam in pivot_df.columns]
        pivot_df = pivot_df[cols]
        
        # Output file name
        coverage_output_file = f"{output_dir}/{args.outfile}_coverage.tsv"
        
        # Open the file to write the header and data
        with open(coverage_output_file, 'w') as f:
            # Write custom header lines
            f.write("# ArDu coverage outfile\n")
            f.write("# Generated on: " + str(pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')) + "\n")
            f.write("# Format LOCI  uDoC:sdDoC:medDoC:CovBases:RCN\n")
            f.write("# uDOC: mean Depth of Coverage over the entire locus span.\n")
            f.write("# sdDOC:standard deviation of mean Depth of Coverage.\n")
            f.write("# medDoC:median Depth of Coverage over the entire locus span.\n")
            f.write("# RCN: relative copy number.\n")
            # Write the actual data from the DataFrame
            pivot_df.to_csv(f, sep='\t', index=False)
            

        print(f"Coverage data saved to {output_dir}/{coverage_output_file}")
    else:
        print("No coverage data to write.")


    # -----------------------------
    # Combines the breakpoints df into a single one and write to out 
    if breakpoints_data:
        bp_df = pd.DataFrame(breakpoints_data)
        breakpoints_output_file = f"{output_dir}/{args.outfile}_breakpoints.tsv"
        bp_df.to_csv(breakpoints_output_file, sep='\t', index=False)
        print(f"Breakpoints data saved to {output_dir}/{breakpoints_output_file}")



    # -----------------------------
    if args.mutation:
        mutation_data = {}
        with open(args.mutation, 'r') as mutation_file:
            for line in mutation_file:
                parts = line.strip().split('\t')
                chromosome = parts[0]
                position = parts[1]
                mutation_name = parts[2] if len(parts) > 2 else "noname"
                mutation_id = f"{chromosome}:{position}:{mutation_name}"
                mutation_data[mutation_id] = {}
                for bam_file in bam_files:
                    try:
                        nucleotide_counts, total_depth = calculate_nucleotide_counts(bam_file, chromosome, position)
                        mutation_data[mutation_id][bam_file] = {
                            'A': nucleotide_counts['A'],
                            'T': nucleotide_counts['T'],
                            'C': nucleotide_counts['C'],
                            'G': nucleotide_counts['G'],
                            'depth': total_depth
                        }
                    except Exception as e:
                        print(f"Error processing {mutation_id} in {bam_file}: {e}")
                        mutation_data[mutation_id][bam_file] = {
                            'A': 'NA',
                            'T': 'NA',
                            'C': 'NA',
                            'G': 'NA',
                            'depth': 'NA'
                        }
        mutation_output_file = f"{output_dir}/{args.outfile}_mutations.tsv"
        with open(mutation_output_file, 'w') as out_file:
            out_file.write("mutation_id\t" + "\t".join(bam_files) + "\n")
            for mutation_id, bam_counts in mutation_data.items():
                row = [mutation_id]
                for bam_file in bam_files:
                    if bam_file in bam_counts:
                        counts = bam_counts[bam_file]
                        count_str = f"A={counts['A']};T={counts['T']};C={counts['C']};G={counts['G']};depth={counts['depth']}"
                    else:
                        count_str = "A=NA;T=NA;C=NA;G=NA;depth=NA"
                    row.append(count_str)
                out_file.write("\t".join(row) + "\n")
        print(f"Mutation data saved to {output_dir}/{mutation_output_file}")

    elapsed_time = round(time.time() - start_time, 1)
    print(f"Elapsed time: {elapsed_time} seconds")

if __name__ == "__main__":
    main()
