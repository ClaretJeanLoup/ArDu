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

        
# Function to normalize the depth of coverage per region
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
                print(f"Mean coverage for gene '{name}' is zero; cannot normalize.")
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



def main():
    parser = argparse.ArgumentParser(description="Process BAM files for coverage and mutation analysis.")
    parser.add_argument("-b", "--bam", required=True, help="File with list of BAM files.")
    parser.add_argument("-r", "--region", required=True, help="Regions file in BED format.")
    parser.add_argument("-n", "--norm", required=True, help="Normalization region.")
    parser.add_argument("-o", "--outfile", required=True, help="Output file prefix.")
    parser.add_argument("--plot", action="store_true", help="Enable plotting.")
    parser.add_argument("--plot_force", action="store_true", help="Force plotting all genes.")
    parser.add_argument("--plot_threshold", type=float, default=0.0, help="Threshold for plotting.")
    parser.add_argument("--plot_interval", help="Interval for plotting.")
    parser.add_argument("--plot_proportion", type=float, help="Proportion for plotting.")
    parser.add_argument("--plot_min_norm_depth", type=float, help="Minimum normalized depth for plotting.")
    parser.add_argument("--plot_max_norm_depth", type=float, help="Maximum normalized depth for plotting.")
    parser.add_argument("--plot_slw", type=int, default=10, help="Window size for smoothing plot.")
    parser.add_argument("--plot_gene_pos", action="store_true", help="Plot gene positions.")
    parser.add_argument("--breakpoint", help="Breakpoint detection method.")
    parser.add_argument("--bkp_model", help="Breakpoint model.")
    parser.add_argument("--bkp_slw", type=int, help="Breakpoint window size.")
    parser.add_argument("--bkp_pen", type=float, help="Breakpoint penalty.")
    parser.add_argument("--bkp_nb", type=int, help="Number of breakpoints.")
    parser.add_argument("--bkp_threshold", type=float, help="Breakpoint threshold.")
    parser.add_argument("--bkp_passes", type=int, default=1, help="Number of passes for breakpoint detection.")
    parser.add_argument("--mutation", help="Mutation file.")
    args = parser.parse_args()

    start_time = time.time()  # Start timing the execution

    # Check if the provided argument is a valid file path
    with open(args.bam, 'r') as f:
        bam_files = f.read().splitlines()

    # Create dictionary for the regions to screen 
    regions_dict = get_regions(args.region)

    # Check if normalization region is present in the regions file
    if args.norm not in regions_dict.keys():
        raise ValueError(f"Normalization region '{args.norm}' not found in the regions file.")

    # Dictionary to store the coverage for all genes across all BAM files
    all_coverage = {gene: {'normalised': [], 'mean': [], 'sd': []} for gene in regions_dict.keys()}

    # Process each BAM file
    for bam_file in bam_files:
        try:
            # Compute coverage values for each gene
            bam_name = bam_file[:-len('.bam')]
            gene_coverage = calculate_mean_and_sd_coverage_for_genes(regions_dict, bam_file)
            
            # Calculate reference coverage for normalization
            normalisation_region = regions_dict[args.norm]
            ref_coverage = calculate_mean_and_sd_coverage_for_genes({args.norm: normalisation_region}, bam_file)
            ref_mean = list(ref_coverage.values())[0]['mean']
            
            # Normalize gene coverage
            norm_gene_coverage = region_normalisation(gene_coverage, ref_mean)

            if isinstance(ref_mean, str):
                print(f"Could not compute a normalized depth of coverage in {bam_file}.")

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
                continue  # Move to the next BAM file
            else:
                if args.plot_force:
                    # All genes should be included if plot_force is True
                    filtered_genes = coverage_df
                    filtered_genes = filtered_genes.drop(index=[args.norm])
                else:
                    # Only include genes that meet the threshold
                    coverage_df['meets_threshold'] = coverage_df['normalised'].apply(lambda x: x[0] > args.plot_threshold if isinstance(x, list) and x else x > args.plot_threshold)
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

                    # Normalize depth data
                    d['depth'] = pd.to_numeric(d['depth'], errors='coerce')
                    d["norm"] = d.depth / ref_mean

                    # Smooth the normalized depth of coverage using numpy window average
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
                                    f.write(f"{args.region.split(':')[0]}\t{d.iloc[b].pos}\t{line_num}\n")
                    elif args.breakpoint == "manseg":
                        results = detect_shifts(d, passes=args.bkp_passes, threshold=args.bkp_threshold)
                        if results:
                            with open(f"{bam_name}.{gene}.breakpoints_manseg.tsv", "w") as f:
                                f.write("Chromosome\tPosition\tNumber\n")
                                for result in results:
                                    for num, start, end in result:
                                        f.write(f"{args.region.split(':')[0]}\t{d.iloc[start].pos}\t{num}\n")

                    # Generate plot
                    fig, ax = plt.subplots()
                    d.plot(x="pos", y="moving_average", ax=ax)
                    ax.set_xlabel("Position")
                    ax.set_ylabel("Normalized Depth of Coverage")

                    if args.plot_gene_pos:
                        for i in range(len(d)):
                            ax.axvline(d.iloc[i].pos, color="gray", linestyle="--", alpha=0.3)

                    fig.savefig(f"{bam_name}.{gene}.png")
                    plt.close(fig)

        except Exception as e:
            print(f"Error processing {bam_file}: {e}")
            continue

    # Save the coverage data to a TSV file
    coverage_df.to_csv(f"{args.outfile}.coverage.tsv", sep='\t')

    # Mutation analysis (if specified)
    if args.mutation:
        mutation_positions = get_regions(args.mutation)
        mutation_counts = {mut: {bam: 0 for bam in bam_files} for mut in mutation_positions.keys()}

        for bam_file in bam_files:
            bam_name = bam_file[:-len('.bam')]
            try:
                counts = calculate_nucleotide_counts(bam_file, mutation_positions)
                for mut in counts.keys():
                    mutation_counts[mut][bam_name] = counts[mut]
            except Exception as e:
                print(f"Error processing mutations in {bam_file}: {e}")
                continue
        
        # Convert mutation counts to DataFrame and save
        mutation_df = pd.DataFrame(mutation_counts).T
        mutation_df.to_csv(f"{args.outfile}.mutation.tsv", sep='\t')

    elapsed_time = time.time() - start_time
    print(f"Elapsed time: {elapsed_time:.2f} seconds")

if __name__ == "__main__":
    main()
