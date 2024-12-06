#!/usr/bin/env python3

import pandas as pd
import os
import re
import argparse

# Function to find matching directories with variants_unique_ann.csv files
def find_matching_files(directory):
    matching_sets = {}
    
    # Walk through all subdirectories
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file == 'variants_unique_ann.csv':
                # Extract the parent directory name
                parent_dir = os.path.basename(os.path.dirname(root))
                
                # Match only directories ending with _E3_1, _E3_2, or _E3_3
                match = re.match(r"(.*)(_E3_[123])$", parent_dir)
                if match:
                    sample_name = match.group(1)
                    suffix = match.group(2)
                    
                    # Collect files with similar sample names and specific suffixes
                    if sample_name not in matching_sets:
                        matching_sets[sample_name] = []
                    matching_sets[sample_name].append((os.path.join(root, file), suffix))
    
    return matching_sets

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Filter and merge variant tables based on matching samples.")
parser.add_argument('-d', '--directory', type=str, required=True, help="Root directory to search for files")
args = parser.parse_args()

# Find matching files in the given directory
matching_files = find_matching_files(args.directory)

# If no matching files are found, exit
if not matching_files:
    exit(1)

# Process each matching set of files
for sample_name, file_list in matching_files.items():
    # Sort file list based on suffix (_E3_1, _E3_2, _E3_3) to maintain order
    file_list = sorted(file_list, key=lambda x: x[1])

    # List to store dataframes after filtering for merging
    dataframes = []
    
    for file_path, suffix in file_list:
        # Load the CSV file
        df = pd.read_csv(file_path)
        
        # Filter based on protein_start range for each suffix
        if suffix == '_E3_1':
            df_filtered = df[(df['protein_start'] >= 187) & (df['protein_start'] <= 210)]
        elif suffix == '_E3_2':
            df_filtered = df[(df['protein_start'] >= 211) & (df['protein_start'] <= 241)]

        # Save the filtered dataframe to a new file in the same directory
        filtered_file_path = os.path.join(os.path.dirname(file_path), 'variants_unique_ann_filtered.csv')
        df_filtered.to_csv(filtered_file_path, index=False)

        # Calculate MAF using 'ct' and 'depth' columns
        df_filtered = df_filtered.copy()
        df_filtered['MAF'] = df_filtered['ct'] / df_filtered['depth']

        # Append the filtered dataframe to the list for merging
        dataframes.append(df_filtered)

    # Concatenate dataframes for vertical merging
    merged_df = pd.concat(dataframes, ignore_index=True)

    # Create the output filename dynamically based on the sample name
    output_file = f"merged_table_{sample_name}_E3.csv"
    
    # Save the merged DataFrame to a CSV file
    merged_df.to_csv(output_file, index=False)