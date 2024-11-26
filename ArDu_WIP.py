# Loading required packages 
import argparse
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pysam
import ruptures as rpt
import statistics 
import subprocess
import time

# Function for rolling average smoothing
def window_average(arr: np.ndarray, w: int):
    numbers_series = pd.Series(arr)
    windows = numbers_series.rolling(w, center=True)  # Use center=True to align the window
    moving_averages = windows.mean()
    return moving_averages


# Function to calculate the sum coverage from a BAM file within a specified region
def calculate_sum_coverage(bam_file, region):
    try:
        # Open BAM file
        samfile = pysam.AlignmentFile(bam_file, "rb")
        coverage_depths = []

        # Calculate coverage depths for the region
        for pileupcolumn in samfile.pileup(region=region, stepper="all"):
            coverage_depths.append(pileupcolumn.n)
        samfile.close()

        # Calculate sum coverage and region length
        if coverage_depths:
            sum_coverage = sum(coverage_depths)
            region_length = len(coverage_depths)
            return sum_coverage, region_length, coverage_depths  # Return individual coverage depths
        else:
            # Return "NA" if no coverage data is available
            return "NA", "NA", []

    except Exception as e:
        print(f"Error calculating sum coverage for region {region}: {e}")
        return "NA", "NA", []  # Return NA for errors

# Function to read BED file and calculate mean and standard deviation of coverage for each gene
def calculate_mean_and_sd_coverage_for_genes(regions, bam_file):
    gene_coverage = {}
    
    for name, region_list in regions.items():
        all_coverages = []  # To store individual coverage values for calculating SD
        
        for region in region_list:
            region_coverage, region_length, coverage_depths = calculate_sum_coverage(bam_file, region)
            
            if region_coverage == 'NA' or region_length == 'NA':
                mean_coverage = 'NA'
                sd_coverage = 'NA'
                break
            
            all_coverages.extend(coverage_depths)  # Collect all coverage values
            
        else:  # This else corresponds to the for loop
            # Calculate mean and standard deviation
            if all_coverages:
                mean_coverage = sum(all_coverages) / len(all_coverages) if len(all_coverages) > 0 else 'NA'
                sd_coverage = statistics.stdev(all_coverages) if len(all_coverages) > 1 else 0.0  # SD is 0 if only one value
            else:
                mean_coverage = 'NA'
                sd_coverage = 'NA'
        
        gene_coverage[name] = {'mean': mean_coverage, 'sd': sd_coverage}

    return gene_coverage

        
# Function to normalise the depth of coverage per region
def region_normalisation(gene_coverage, ref_coverage):
    norm_coverage = {}
    
    for name in gene_coverage.keys():
        # Access the mean coverage from gene_coverage
        mean_coverage = gene_coverage[name]['mean']
        
        # Check the type of the value in ref_coverage corresponding to the key `name`
        if isinstance(ref_coverage, float) and isinstance(mean_coverage, float) and mean_coverage != 0:
            normalised_coverage = mean_coverage / ref_coverage  
        else:
            normalised_coverage = "NA" 
            if mean_coverage == 0:
                print(f"Mean coverage for gene '{name}' is zero; cannot normalise.")
            else:
                print(f"Unexpected type for gene_coverage[{name}]: {type(mean_coverage)}")
       
        # Store the normalised coverage in the dictionary
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
    # Initialise variables to store the gene's intervals
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
    parser.add_argument("-n","--norm", required=True, help="Name of the region to use for normalisation describe in the regions file.")
    parser.add_argument("-o", "--outfile", required=True, help="Prefix for output file names")

    #Optional arguments
    ## Plotting
    parser.add_argument("--plot", type=str, default=None, help="Optional: produces a plot of the raw and normalised depth of coverage. Accepts as value the following extensions for the plot file: \
                        png, jpeg, jpg, pdf, svg, eps. If used jointly with --breakpoint, the putative breakpoints will be plotted.")
    parser.add_argument("--plot_threshold", type=float, default= 1.4, help="Treshold value to consider a duplication. Default value is 1.4, i.e. the expected normalised coverage of a heterozygous diploid organism possessing a mono-copy and a duplicated copy of a given gene.")
    parser.add_argument("--plot_interval", type=str, help="A specific genomic interval used for plotting, format chromosome:start-stop")
    parser.add_argument("--plot_proportion", type=float, default= 2, help="Set the plotting interval to X time the size of the total gene span. Default = 2.")
    parser.add_argument("--plot_slw", type=int, default= 1000, help="Sliding window size (pb) used for coverage representation, higher values will result in a smoother depth of coverage representation.")
    parser.add_argument("--plot_min_norm_depth", type=float, help="Minimum normalised depth of coverage to use for plotting.")
    parser.add_argument("--plot_max_norm_depth", type=float, help="Maximum normalised depth of coverage to use for plotting.")
    parser.add_argument("--plot_gene_pos",action='store_true', help="When set, add vertical lines to indicate the duplicated gene position on the plot.")
    parser.add_argument("--plot_force", action='store_true', help="Forces plotting regardless of depth of coverage value.")


    ## Breakpoints detection
    parser.add_argument("--breakpoint", choices=["ruptures", "rollingaverage"], help="Optional: assess putative breakpoints. Two methods are available: 'ruptures' and 'rollingaverage'. \n\
                        'ruptures' uses the ruptures python package (see https://centre-borelli.github.io/ruptures-docs/). Associated options are --window-size, --pen and --model. \n\
                        'rollingaverage' uses successive rolling average to detect shifts in depth of coverage. Associated options are --window-size and --threshold and --passes. \n\
                        Setting this argument will create a '.breakpoints.tsv' output file containing the positions found by the rupture package, or the regions of depth of coverage \
                        shifting found with 'rollingaverage'.")
    parser.add_argument("--bkp_slw", type=int, default=1000, help="Window size for ruptures model and rolling average calculation. Default= 1 kb.")
    parser.add_argument("--bkp_nb", type=int, help="Number of expected breakpoints in this structure, cannot be used alongside --bkp_pen.")
    parser.add_argument("--bkp_pen", type=int, help="Penalty parameter for ruptures algo.predict. A higher value will result in a higher penalty in breakpoint creation.\n\
                        use this option if you don't have a prior on the bkp number in the structure you're looking at.")
    parser.add_argument("--bkp_model", type=str, default="l2", help="Model used by ruptures package. Default ='l2'.")
    parser.add_argument("--bkp_threshold", type=float, default=1.0, help="Threshold for detecting shifts in depth of coverage.")
    parser.add_argument("--bkp_passes", type=int, default=1, help="Number of rolling average passes. Increasing will lessen the overal variation in depth of coverage.")

    ## Genotyping 
    parser.add_argument("--mutation", help="Optional: if set, this option will return the number of reads supporting each nucleotides at the given position(s) in a .mutation.tsv file. Takes as input a tab-delimited file\
                         with the following columns: chromosome, position and an optional third column containing the name of the muations you're screening for (tab-delimited).")

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
    all_coverage = {gene: {'normalised': [], 'mean': [], 'sd': []} for gene in regions_dict.keys()}

    # Process each BAM file
    for bam_file in bam_files:
        try:
            # Compute coverage values for each gene
            bam_name = bam_file[:-len('.bam')]
            print(f"Processing {bam_name}")
            gene_coverage = calculate_mean_and_sd_coverage_for_genes(regions_dict, bam_file)

            # Calculate reference coverage for normalisation
            normalisation_region = regions_dict[args.norm]
            ref_coverage = calculate_mean_and_sd_coverage_for_genes({args.norm: normalisation_region}, bam_file)
            ref_mean = list(ref_coverage.values())[0]['mean']

            # normalise gene coverage
            norm_gene_coverage = region_normalisation(gene_coverage, ref_mean)

            if isinstance(ref_mean, str):
                print(f"Could not compute a normalised depth of coverage in {bam_file}.")

            # Append coverage values to the all_coverage dictionary
            for gene in gene_coverage.keys():
                mean_coverage = gene_coverage[gene]['mean']
                sd_coverage = gene_coverage[gene]['sd']
                normalised_coverage = norm_gene_coverage.get(gene, 'NA')

                all_coverage[gene]['normalised'].append(normalised_coverage)
                all_coverage[gene]['mean'].append(mean_coverage)
                all_coverage[gene]['sd'].append(sd_coverage)

                # Create a DataFrame for easier handling
                coverage_df = pd.DataFrame(all_coverage).T  # Transpose to have genes as rows

                # Apply threshold check
                if not args.plot:
                    continue # Move on to the next BAM
                else:
                    if args.plot_force:
                        # All genes should be included if plot_force is True
                        filtered_genes = coverage_df
                        filtered_genes = filtered_genes.drop(index=[args.norm])

                    else:
                        # Only include genes that meet the threshold
                        coverage_df['meets_threshold'] = coverage_df['normalised'].apply(lambda x: x[0] > args.plot_threshold if isinstance(x, list) and len(x) > 0 else (x > args.plot_threshold if isinstance(x, (int, float)) else False))
                        filtered_genes = coverage_df[coverage_df['meets_threshold']]

                    # Process filtered genes
                    for gene in filtered_genes.index:
                        mean_coverage = filtered_genes.at[gene, 'mean']
                        sd_coverage = filtered_genes.at[gene, 'sd']
                        normalised_coverage = filtered_genes.at[gene, 'normalised']

                        # Depth file creation to produce a plot
                        if args.plot_interval and args.plot:
                            str_incr_span = args.plot_interval
                            samfile = pysam.AlignmentFile(bam_file, "rb")
                            depth_data = [(p.pos, p.n) for p in samfile.pileup(region=str_incr_span)]

                        elif args.plot_proportion and args.plot:
                            totalspan = get_totalspan(args.region, gene)
                            if totalspan is None:
                                print(f"Warning: No intervals found for gene {gene} in BED file.")
                                continue
                            chromosome, plotStart, plotStop = totalspan
                            # Adjust plotStart and plotStop according to args.plot_proportion
                            gene_length = plotStop - plotStart
                            # Calculate extension based on the proportion of total gene length
                            extension = gene_length * args.plot_proportion
                            plotStart -= int(extension)
                            plotStop += int(extension)
                            # Create depth file
                            str_incr_span = f"{chromosome}:{plotStart}-{plotStop}"
                            samfile = pysam.AlignmentFile(bam_file, "rb")
                            depth_data = [(p.pos, p.n) for p in samfile.pileup(region=str_incr_span)]

                        else:
                            continue

                        samfile.close()
                        d = pd.DataFrame(depth_data, columns=["pos", "depth"])  # Convert depth data into a DataFrame

                        # normalise depth data
                        d['depth'] = pd.to_numeric(d['depth'], errors='coerce')
                        d["norm"] = d.depth/ref_mean

                        # Smooth the normalised depth of coverage using numpy window average
                        d["moving_average"] = window_average(d["norm"].to_numpy(), args.plot_slw)

                        # Remove under- or over-covered bases
                        if args.plot_min_norm_depth:
                            d = d[d["norm"] >= args.plot_min_norm_depth]

                        if args.plot_max_norm_depth:
                            d = d[d["norm"] <= args.plot_max_norm_depth]

                        # Breakpoints identification
                        if args.breakpoint:
                            print(f"Breakpoint method: {args.breakpoint}")

                            if args.breakpoint == "ruptures":
                                if args.bkp_nb:
                                    print("Nb of breakpoints specified")
                                elif args.bkp_pen:
                                    print("Penalty specified")
                                else:
                                    print("Error: neither bkp_nb nor bkp_pen is set")

                                algo = rpt.Window(model=args.bkp_model, width=args.bkp_slw).fit(d["depth"].to_numpy().reshape(-1, 1))

                                if args.bkp_pen:
                                    print(f"Penalty set to {args.bkp_pen}")
                                    result = algo.predict(pen=args.bkp_pen)
                                elif args.bkp_nb:
                                    print(f"Expecting {args.bkp_nb} breakpoints")
                                    result = algo.predict(n_bkps=args.bkp_nb)
                                else:
                                    result = []

                                # Write genomic positions of detected breakpoints to the output file
                                if result:
                                    with open(f"{bam_name}.{gene}.breakpoints_ruptures.tsv", "w") as f:
                                        f.write("Chromosome\tPosition\tNumber\n")
                                        for line_num, b in enumerate(result[:-1], start=1):
                                            f.write(f"{args.region.split(':')[0]}\t{int(d.loc[b]['pos'])}\t{line_num}\n")

                            elif args.breakpoint == "rollingaverage":
                                significant_shifts = detect_shifts(d, args.bkp_threshold, args.bkp_slw, args.bkp_passes)

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
                        plt.plot(d.pos, d.norm, label="Normalised Depth of Coverage", color="lightsteelblue", linestyle='none', marker='.')
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
            print(f"{bam_file}: {e}")

    # Write the results to a single TSV file
    coverage_output_file = f"{args.outfile}_coverage.tsv"
    with open(coverage_output_file, 'w') as f:
        # Remove the .bam extension from bam_file names
        bam_names_without_ext = [bam_name.replace('.bam', '') for bam_name in bam_files]

        # Create header with proper grouping, using names without the .bam extension
        header = "#FORMAT MEAN_DEPTH:SD:NORMALISED_DEPTH" + "\n" + "gene_name\t" + "\t".join([f"{bam_name}" for bam_name in bam_names_without_ext]) + "\n"
        f.write(header)

        # Loop through all genes and write their coverage data
        for gene, coverages in all_coverage.items():
            # Initialize a list to hold line parts for the current gene
            line_parts = [gene]

            # Loop through the coverage data for each BAM file
            for mean, sd, norm in zip(coverages['mean'], coverages['sd'], coverages['normalised']):
                # Prepare the combined string for the current BAM file
                mean_str = str(round(mean, 2)) if isinstance(mean, float) else str(mean)
                sd_str = str(round(sd, 2)) if isinstance(sd, float) else str(sd)
                norm_str = str(round(norm,1)) if norm != 'NA' else 'NA'

                # Combine mean, SD, and normalised values into a single string
                combined_str = f"{mean_str}:{sd_str}:{norm_str}"
                line_parts.append(combined_str)

            # Join the line_parts with tabs and write to the file
            line = "\t".join(line_parts) + "\n"
            f.write(line)

    if args.mutation:
        # Initialise a dictionary to hold mutation counts per mutation_id across all BAM files
        mutation_data = {}

        # Read mutation positions and chromosomes from the provided file
        with open(args.mutation, 'r') as mutation_file:
            for line in mutation_file:
                parts = line.strip().split('\t')
                chromosome = parts[0]
                position = parts[1]
                mutation_name = parts[2] if len(parts) > 2 else "noname"  # Default to "noname"

                # Create a unique identifier for each mutation
                mutation_id = f"{chromosome}:{position}:{mutation_name}"

                # Initialize an empty dictionary for each mutation_id to store BAM file counts
                mutation_data[mutation_id] = {}

                for bam_file in bam_files:
                    try:
                        nucleotide_counts, total_depth = calculate_nucleotide_counts(bam_file, chromosome, position)

                        # Store counts in the dictionary
                        mutation_data[mutation_id][bam_file] = {
                            'A': nucleotide_counts['A'],
                            'T': nucleotide_counts['T'],
                            'C': nucleotide_counts['C'],
                            'G': nucleotide_counts['G'],
                            'depth': total_depth
                        }

                    except Exception as e:
                        print(f"Error processing {mutation_id} in {bam_file}: {e}")
                        # If there's an error, store 'NA' to indicate missing data
                        mutation_data[mutation_id][bam_file] = {
                            'A': 'NA',
                            'T': 'NA',
                            'C': 'NA',
                            'G': 'NA',
                            'depth': 'NA'
                        }

        # Write the mutation data to the output file in a matrix format
        mutation_output_file = f"{args.outfile}_mutations.tsv"
        with open(mutation_output_file, 'w') as out_file:
            # Header row: mutation IDs followed by BAM file names
            out_file.write("mutation_id\t" + "\t".join(bam_files) + "\n")

            # For each mutation, write mutation counts and depths for each BAM file
            for mutation_id, bam_counts in mutation_data.items():
                row = [mutation_id]
                for bam_file in bam_files:
                    if bam_file in bam_counts:
                        counts = bam_counts[bam_file]
                        # Format counts as A=5;T=10;C=3;G=2;depth=20
                        count_str = f"A={counts['A']};T={counts['T']};C={counts['C']};G={counts['G']};depth={counts['depth']}"
                    else:
                        count_str = "A=NA;T=NA;C=NA;G=NA;depth=NA"
                    row.append(count_str)
                out_file.write("\t".join(row) + "\n")

        print(f"Mutation data saved to {mutation_output_file}")
    print(f"Coverage data saved to {coverage_output_file}")

    elapsed_time = round(time.time() - start_time,1)
    print(f"Elapsed time: {elapsed_time} seconds")



if __name__ == "__main__":
    main()
