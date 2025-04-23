#!/usr/bin/env python3

import os
import glob
import csv
import logging
import sys
import argparse

def remove_na_columns(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    for csv_file in glob.glob(os.path.join(input_dir, '*.csv')):
        try:
            output_file = os.path.join(output_dir, os.path.basename(csv_file))
            
            with open(csv_file, 'r') as f:
                reader = list(csv.reader(f))
                if not reader:
                    continue
                    
                header = reader[0]
                # Find all columns containing NA
                na_columns = [i for i in range(len(header)) if any(row[i] == 'NA' for row in reader[1:])]
                
            # Remove NA columns
            new_rows = []
            for row in reader:
                new_rows.append([cell for i, cell in enumerate(row) if i not in na_columns])
                
            # Write the new file without NA columns
            with open(output_file, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerows(new_rows)
                
            # Count removed and remaining columns
            removed_count = len(na_columns)
            remaining_count = len(header) - removed_count
            print(f"Processed {csv_file}: Removed {removed_count} NA columns. Remaining columns: {remaining_count}.")
        except Exception as e:
            print(f"Error processing file {csv_file}: {e}")

def main():
    parser = argparse.ArgumentParser(description='Remove columns with NA values from CSV files.')
    parser.add_argument('--input', type=str, required=True,
                       help='Input directory containing CSV files.')
    parser.add_argument('--output', type=str, required=True,
                       help='Output directory for processed CSV files.')
    args = parser.parse_args()
    
    remove_na_columns(args.input, args.output)
    print("NA column removal completed.")

if __name__ == "__main__":
    main()
