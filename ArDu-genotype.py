# Loading required packages 
import argparse
import os
import pysam
import subprocess
import time

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


# Main function
def main():
    # Mandatory arguments
    parser = argparse.ArgumentParser(description="Genotype sample by screening diagnostic point mutations")
    parser.add_argument("-b", "--bam", help="Input a list of BAM files (one per line)")
    parser.add_argument("-o", "--outfile", required=True, help="Prefix for output file names")
    parser.add_argument("-m","--mutation", required=True, help="Takes as input a tab-delimited file\
                         with the following columns: chromosome, position and an optional third column containing the name of the muations you're screening for.")


    args = parser.parse_args()

    start_time = time.time()  # Start timing the execution

    # Check if the provided argument is a valid file path
    with open(args.bam, 'r') as f:
        bam_files = f.read().splitlines()

        # Initialize a dictionary to hold mutation counts per mutation_id across all BAM files
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

    elapsed_time = round(time.time() - start_time,1)
    print(f"Elapsed time: {elapsed_time} seconds")

if __name__ == "__main__":
    main()