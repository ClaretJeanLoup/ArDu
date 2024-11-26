# Loading required packages 
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import ruptures as rpt
import argparse
import time
import subprocess
import os

# Function to normalise depth of coverage
def normalise_depth(depth_data, ref_cover):
    return depth_data / ref_cover

# Function for rolling average smoothing
def window_average(arr: np.ndarray, w: int):
    numbers_series = pd.Series(arr)
    windows = numbers_series.rolling(w, center=True)  # Use center=True to align the window
    moving_averages = windows.mean()
    return moving_averages

# Function to calculate the sum coverage from a BAM file within a specified region
def calculate_sum_coverage(bam_file, region):
    samfile = pysam.AlignmentFile(bam_file, "rb")
    coverage_depths = []
    for pileupcolumn in samfile.pileup(region=region):
        coverage_depths.append(pileupcolumn.n)
    samfile.close()
    if coverage_depths:
        sum_coverage = sum(coverage_depths) 
        length = len(coverage_depths)
        return sum_coverage, length 
        
    else:
        return 0, 0

# Function to read BED file and calculate mean coverage for each gene
def calculate_mean_coverage_for_genes(regions, bam_file):
    gene_coverage = {}
    for name in regions.keys():
        for region in regions[name]:
            sum_coverage, length = calculate_mean_coverage(bam_file, region)
            if name in gene_coverage:
                gene_coverage[name].append((sum_coverage,length))
            else:
                gene_coverage[name] = [(sum_coverage,length)]
        sum_coverage = 0
        length = 0
        for region in gene_coverage[name]: 
            sum_coverage += region[0]
            length += region[1]
        mean_coverage = sum_coverage / length if length != 0 else 'NA'
        gene_coverage[name] = mean_coverage
    return gene_coverage

# Function to normalise the depth of coverage per regions
def region_normalisation(gene_coverage, ref_coverage):
    norm_coverage = {}
    for name in gene_coverage.keys():
        if ref_coverage: 
            normalised_coverage = gene_coverage[name] / ref_coverage
        else:
            normalised_coverage = "NA"
        # Append normalised coverage to the list for the gene or create a new list
        norm_coverage[name] = normalised_coverage
    return norm_coverage

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
    position = int(position)  # Ensure position is an integer
    samfile = pysam.AlignmentFile(bam_file, "rb")
    counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    total_depth = 0
    for pileupcolumn in samfile.pileup(region=f"{chromosome}:{position}-{position}"):
        if pileupcolumn.pos == position - 1:  # Adjust for 0-based indexing
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
            if not line_str.strip():  # Skip empty lines
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
    # Initialize variables to store the gene's intervals
    gene_intervals = []
    
    # Read the bed file
    with open(regions, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) != 4:
                continue  # Skip malformed lines
            
            chromosome = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            gene = parts[3]
            
            # Check if the gene matches the specified gene name
            if gene == gene_name:
                gene_intervals.append((start, end, chromosome))
                chrm = chromosome
    
    # If gene_intervals is empty, the gene was not found in the BED file
    if not gene_intervals:
        return None
    
    # Calculate the minimum start and maximum end for the gene intervals
    min_start = min(interval[0] for interval in gene_intervals)
    max_end = max(interval[1] for interval in gene_intervals)
    
    # Return the interval as a tuple in the format (chromosome, min_start, max_end)
    totalspan = (chrm, min_start, max_end)
    return totalspan

# Main function
def main():
    # Mandatory arguments
    parser = argparse.ArgumentParser(description="Identify duplication structure via depth of coverage data")
    parser.add_argument("-b", "--bam", help="Input a list of BAM files (one per line)")
    parser.add_argument("-r", "--region", required=True, help="Takes as input a 4 columns file (tab delimited, as bed format: 'chromosome\tstart\tstop\tname') containing regions of \
                        interest. Features like genes composed of multiple regions (e.g. exons) will be pooled together if they share the same name (in the 4th column)")
    parser.add_argument("--norm", required=True, help="Name of the region to use for normalisation describe in the regions file.")
    parser.add_argument("-o", "--outfile", required=True, help="Prefix for output file names")

    #Optional arguments
    ## Plotting
    parser.add_argument("--plot", type=str, default="png",help="Optional: produces a plot of the raw and normalised depth of coverage. Accepts as value the following extensions for the plot file: \
                        png, jpeg, jpg, pdf, svg, eps. If used jointly with --breakpoint, the putative breakpoints will be plotted.")
    parser.add_argument("--plot_treshold", type=float, default= 1.4, help="Treshold value to consider a duplication. Default value is 1.4, i.e. the expected normalised coverage of a heterozygous diploid organism possessing a mono-copy and a duplicated copy of a given gene.")
    parser.add_argument("--plot_extension", type=int, default= 1000, help="A number of bases that will be added on either sides of the gene interval for plotting")
    parser.add_argument("--plot_interval", type=str, help="A specific genomic interval used for plotting, format chromosome:start-stop")
    parser.add_argument("--plot_proportion", type=float, default= 2, help="Set the plotting interval to X time the size of the total gene span. Default = 2.")
    parser.add_argument("--plot_slw", type=int, default= 1000, help="Sliding window size (pb) used for coverage representation, higher values will result in a smoother depth of coverage representation.")
    parser.add_argument("--plot_min_norm_depth", type=float, help="Minimum normalised depth of coverage to use for plotting.")
    parser.add_argument("--plot_max_norm_depth", type=float, help="Maximum normalised depth of coverage to use for plotting.")
    parser.add_argument("--plot_gene_pos",action='store_true', help="When set, add vertical lines to indicate the duplicated gene position on the plot.")

    ## Breakpoints detection
    parser.add_argument("--breakpoint", choices=["ruptures", "rollingaverage"], help="Optional: assess putative breakpoints. Two methods are available: 'rupture' and 'rollingaverage'. \n\
                        'rupture' uses the ruptures python package (see https://centre-borelli.github.io/ruptures-docs/). Associated options are --window-size, --pen and --model. \n\
                        'rollingaverage' uses successive rolling average to detect shifts in depth of coverage. Associated options are --window-size and --threshold and --passes. \n\
                        Setting this argument will create a '.breakpoints.tsv' output file containing the positions found by the rupture package, or the regions of depth of coverage \
                        shifting found with 'rollingaverage'.")
    parser.add_argument("--window-size", type=int, default=1000, help="Window size for ruptures model and rolling average calculation. Default= 1 kb.")
    parser.add_argument("--bkp-nb", type=int, default=2, help="Number of expected breakpoints in this structure, cannot be used alongside --bkp-pen. Default= 2.")
    parser.add_argument("--bkp-pen", type=int, default=2, help="Penalty parameter for ruptures algo.predict. A higher value will result in a higher penalty in breakpoint creation.\n\
                        use this option if you don't have a prior on the bkp number in the structure you're looking at. Default= 2.")
    parser.add_argument("--model", type=str, default="l2", help="Model used by ruptures package. Default ='l2'.")
    parser.add_argument("--threshold", type=float, default=1.0, help="Threshold for detecting shifts in depth of coverage.")

    ## Genotyping 
    parser.add_argument("--passes", type=int, default=1, help="Number of rolling average passes. Increasing will result in a lessening of the overal variation in depth of coverage.")
    parser.add_argument("--mutation", help="Optional: if set, this option will return the count of bases at the given position(s) in a .mutation.tsv file. Takes as input a file with \
                        chromosome and position pairs (tab-delimited).")

    args = parser.parse_args()

    start_time = time.time()  # Start timing the execution

    # Check if the provided argument is a valid file path
    with open(args.bam, 'r') as f:
        bam_files = f.read().splitlines()
    
    # Create dictionnary for the regions to screen 
    regions_dict = get_regions(args.region)

    # Check if normalisation region is present in the regions file
    if args.norm not in regions_dict.keys():
        raise ValueError(f"Normalisation region '{args.norm}' not found in the regions file.")
    
    # Dictionary to store the coverage for all genes across all BAM files
    all_coverage = {gene: [] for gene in regions_dict.keys()}

    # Process each BAM file
    for bam_file in bam_files:
        try:
            # Compute coverage values for each gene
            bam_name = bam_file[:-len('.bam')]
            gene_coverage = calculate_mean_coverage_for_genes(regions_dict, bam_file)
            normalisation_region = regions_dict[args.norm]
            ref_coverage = calculate_mean_coverage_for_genes({args.norm: normalisation_region}, bam_file)
            ref_coverage = list(ref_coverage.values())[0]
            norm_gene_coverage = region_normalisation(gene_coverage, ref_coverage)

            # Append coverage values to the all_coverage dictionary
            for gene, coverage in norm_gene_coverage.items():
                all_coverage[gene].append(coverage)

            # Check the gene_coverage and processing depth data
            for gene, coverage in norm_gene_coverage.items():
                if coverage > args.plot_treshold:
                    print(f"The {gene} locus in {bam_file} is duplicated.")

                    # Depth file creation to produce a plot
                    if args.plot_interval and args.plot:
                        ch, pos_r = args.plot_interval.split(":")
                        plotStart, plotStop = map(int, pos_r.split("-"))
                        totalspan = f"{ch}:{plotStart}-{plotStop}"
                        samfile = pysam.AlignmentFile(bam_file, "rb")
                        depth_data = [(p.pos, p.n) for p in samfile.pileup(region=str_incr_span)]

                    elif args.plot_proportion and args.plot:
                        totalspan = get_totalspan(args.region, gene)
                        if totalspan is None:
                            print(f"Warning: No intervals found for gene {gene} in BED file.")
                            continue                    
                        chromosome, plotStart, plotStop = totalspan
                        # Adjust plotStart and plotStop according to args.plot_extension
                        gene_length = plotStop - plotStart
                        # Calculate extension based on the proportion of total gene length 
                        extension = gene_length * args.plot_proportion
                        plotStart -= int(extension)
                        plotStop += int(extension)
                        # Create depth file
                        str_incr_span = f"{chromosome}:{plotStart}-{plotStop}"
                        samfile = pysam.AlignmentFile(bam_file, "rb")
                        depth_data = [(p.pos, p.n) for p in samfile.pileup(region=str_incr_span)] 

                    elif args.plot_extension and args.plot:
                        totalspan = get_totalspan(args.region, gene)
                        if totalspan is None:
                            print(f"Warning: No intervals found for gene {gene} in BED file.")
                            continue                    
                        chromosome, plotStart, plotStop = totalspan
                        # Adjust plotStart and plotStop according to args.plot_extension
                        plotStart -= args.plot_extension
                        plotStop += args.plot_extension
                        # Create depth file
                        str_incr_span = f"{chromosome}:{plotStart}-{plotStop}"
                        samfile = pysam.AlignmentFile(bam_file, "rb")
                        depth_data = [(p.pos, p.n) for p in samfile.pileup(region=str_incr_span)]                   

                    else:
                        print("If --plot is set, either --plot_interval, --plot_proportion or --plot_extension should be provided.")
                    
                    samfile.close()
                    d = pd.DataFrame(depth_data, columns=["pos", "depth"])  # Convert depth data into a DataFrame
                    
                    # Normalize depth data
                    d["norm"] = normalise_depth(d.depth, ref_coverage)

                    # Smooth the normalised depth of coverage using numpy window average
                    d["moving_average"] = window_average(d["norm"].to_numpy(), args.plot_slw)

                    # Remove under- or over- covered bases 
                    if args.plot_min_norm_depth:
                        d = d[d["norm"] >= args.plot_min_norm_depth]

                    if args.plot_max_norm_depth:
                        d = d[d["norm"] <= args.plot_max_norm_depth]
                        

                    # Breakpoints identification
                    # Run ruptures on the normalised depth of coverage if --breakpoint is set to "ruptures"
                    if args.breakpoint == "ruptures":
                        algo = rpt.Window(model=args.model, width=args.window_size).fit(d["norm"].to_numpy().reshape(-1, 1))
                        result = algo.predict(pen=args.pen)

                        # Write genomic positions of detected breakpoints to the output file
                        with open(f"{bam_name}.{gene}.breakpoints_ruptures.tsv", "w") as f:
                            f.write("Chromosome\tPosition\tNumber\n")
                            for line_num, b in enumerate(result[:-1], start=1):
                                f.write(f"{args.region.split(':')[0]}\t{int(d.loc[b]['pos'])}\t{line_num}\n")

                    # Run detect_shifts function if --breakpoint is set to "rollingaverage"
                    elif args.breakpoint == "rollingaverage":
                        significant_shifts = detect_shifts(d, args.threshold, args.window_size, args.passes)

                        # Extract consecutive positions and define regions for soft clip identification
                        regions = []
                        consecutive_positions = []
                        for idx, row in significant_shifts.iterrows():
                            if not consecutive_positions:
                                consecutive_positions.append(row['pos'])
                            elif row['pos'] - consecutive_positions[-1] <= 1:  # If consecutive
                                consecutive_positions.append(row['pos'])
                            else:
                                regions.append((consecutive_positions[0], consecutive_positions[-1]))
                                consecutive_positions = [row['pos']]
                        if consecutive_positions:  # Append last region if any consecutive positions left
                            regions.append((consecutive_positions[0], consecutive_positions[-1]))

                        # Writing consecutive regions to an output file
                        with open(f"{bam_name}.{gene}.breakpoints_rollingaverage.tsv", 'w') as outfile:
                            outfile.write("Chromosome\tRegion_Start\tRegion_End\tNumber\n")
                            for i, region in enumerate(regions):
                                outfile.write(f"{args.region.split(':')[0]}\t{int(region[0])}\t{int(region[1])}\t{str(i+1)}\n")




                    # Plotting
                    plt.figure(figsize=(10, 6))
                    plt.plot(d.pos, d.norm, label="Normalised Depth of Coverage", color="lightsteelblue")
                    if args.breakpoint == "ruptures":
                        for i, b in enumerate(result[:-1]):
                            plt.axvline(d.loc[b]["pos"], color="red", linestyle="-", linewidth=0.6) 
                            plt.text(d.loc[b]["pos"], d.loc[b]["norm"], f"{i+1}", color="red", fontsize=10)  # Add breakpoint number
                        plt.title("Depth of Coverage with Ruptures detected breakpoints")
                    elif args.breakpoint == "rollingaverage":
                        for i, region in enumerate(regions):
                            region_center = (region[0] + region[1]) / 2
                            region_data = d[(d['pos'] >= region[0]) & (d['pos'] <= region[1])]
                            if not region_data.empty:
                                max_norm_depth_index = region_data['norm'].idxmax()
                                max_norm_depth = region_data.loc[max_norm_depth_index, 'norm']
                                plt.axvspan(region[0], region[1], color='red', alpha=0.2)
                                plt.text(region_center, max_norm_depth + 0.1, f"{i+1}", color="red", fontsize=10)  # Add region number
                        plt.title("Depth of Coverage with Rolling Average detected breakpoints")
                    else: 
                        plt.title("Depth of Coverage")
                    if args.plot_gene_pos:
                        totalspan = get_totalspan(args.region, gene)
                        plt.axvline(totalspan[1], color="purple", linestyle="-", linewidth=0.8)
                        plt.axvline(totalspan[2], color="purple", linestyle="-", linewidth=0.8)
                    plt.plot(d.pos, d.moving_average, label="Smoothed Depth of Coverage", color="black")
                    plt.suptitle(f"{gene} locus in {bam_file}", fontsize=16)
                    plt.xlabel(f"Genomic position: {str_incr_span}")
                    plt.ylabel("Normalised Depth of Coverage")
                    plt.legend()
                    
                    # Save the plot file
                    plt.savefig(f"{bam_name}_{gene}_plot.{args.plot}")
                    plt.close()

        except Exception as e:
            print(f"Error computing normalised depth of coverage in {bam_file}: {e}")

    # Write the results to a single TSV file
    coverage_output_file = f"{args.outfile}_coverage.tsv"
    with open(coverage_output_file, 'w') as f:
        header = "gene_name\t" + "\t".join(bam_files) + "\n"
        f.write(header)
        for gene, coverages in all_coverage.items():
            line = f"{gene}\t" + "\t".join(map(lambda x: str(round(x,1)), coverages)) + "\n"
            f.write(line)
    print(f"Coverage data saved to {coverage_output_file}")

    elapsed_time = time.time() - start_time
    print(f"Elapsed time: {elapsed_time} seconds")

    # Output mutations if provided
    if args.mutation:
        with open(args.mutation, 'r') as f:
            with open("{args.outfile}.mutations.tsv", 'w') as out_file:
                out_file.write("sample\tchromosome:position\tA\tT\tC\tG\tdepth\n")
                for line in f:
                    chromosome, position = line.strip().split('\t')
                    nucleotide_counts, total_depth = calculate_nucleotide_counts(args.bam, chromosome, position)
                    out_file.write(f"{bam_name}\t{chromosome}:{position}\t{nucleotide_counts['A']}\t{nucleotide_counts['T']}\t{nucleotide_counts['C']}\t{nucleotide_counts['G']}\t{total_depth}\n")
        print("Genotyping file was generated successfully.")

if __name__ == "__main__":
    main()