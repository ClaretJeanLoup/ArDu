import argparse
import pandas as pd
import matplotlib.pyplot as plt

# Define the function to process the coverage data
def process_coverage_data(input_file):
    # Load the data (assuming it's in a tab-separated file)
    df = pd.read_csv(input_file, sep='\t', index_col=0)
    
    # Step 2: Calculate mean depth and standard deviation across samples for each locus
    df['mean_depth'] = df.mean(axis=1)
    df['std_depth'] = df.std(axis=1)
    
    # Step 3: Normalize the depth of coverage for each sample
    normalized_df = df.iloc[:, :-2].div(df['mean_depth'], axis=0)

    # Step 4: Plot distributions
    # Plot distribution of mean depth of coverage
    plt.figure(figsize=(10, 6))
    df['mean_depth'].plot(kind='hist', bins=30, alpha=0.7, color='blue', edgecolor='black')
    plt.title('Distribution of Mean Depth of Coverage')
    plt.xlabel('Mean Depth')
    plt.ylabel('Frequency')
    plt.show()

    # Plot distribution of normalized depth of coverage for the first sample
    plt.figure(figsize=(10, 6))
    normalized_df.iloc[:, 0].plot(kind='hist', bins=30, alpha=0.7, color='green', edgecolor='black')
    plt.title('Distribution of Normalized Depth for Sample 1')
    plt.xlabel('Normalized Depth')
    plt.ylabel('Frequency')
    plt.show()

    # Step 5: Create output file with statistics for each sample
    output = pd.DataFrame()

    # Add mean, standard deviation, and normalized depth for each sample
    for sample in df.columns[:-2]:  # Excluding the 'mean_depth' and 'std_depth' columns
        output[f'{sample}:mean'] = df['mean_depth']
        output[f'{sample}:sd'] = df['std_depth']
        output[f'{sample}:normalized'] = normalized_df[sample]

    # Save the output to a tab-separated file
    output.to_csv('coverage_statistics.txt', sep='\t')
    print("Output file 'coverage_statistics.txt' has been saved.")

# Set up argparse for command-line argument parsing
def main():
    parser = argparse.ArgumentParser(description="Process coverage data for depth analysis.")
    
    # Add an argument for input file
    parser.add_argument('-i', '--input', required=True, help="Path to the input coverage data file (tab-separated).")
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Call the function with the input file
    process_coverage_data(args.input)

# Entry point to run the script
if __name__ == "__main__":
    main()
