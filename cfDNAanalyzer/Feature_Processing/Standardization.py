#!/usr/bin/env python3

import pandas as pd
import os
import argparse
import logging
import sys
from sklearn.preprocessing import RobustScaler
from sklearn.preprocessing import QuantileTransformer
from sklearn.preprocessing import MinMaxScaler

def format_number(value):
    if isinstance(value, (int, float)):
        if value.is_integer():
            return int(value)
        formatted = "{0:.6f}".format(value).rstrip('0').rstrip('.') if '.' in "{0:.6f}".format(value) else str(value)
        return formatted
    return value

def z_score_standardization(input_dir, output_dir, method):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Iterate over all CSV files in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith('.csv'):
            file_path = os.path.join(input_dir, filename)
            print(f"Processing file: {filename}")

            df = pd.read_csv(file_path)

            # Retain 'sample' and 'label' columns
            if 'sample' not in df.columns or 'label' not in df.columns:
                logging.warning(f"File '{filename}' is missing 'sample' or 'label' columns. Skipping this file.")
                continue

            sample_label_df = df[['sample', 'label']]
            other_cols = df.columns.difference(['sample', 'label'])

            # Select numeric columns (excluding 'sample' and 'label')
            numeric_cols = df[other_cols].select_dtypes(include=['float64', 'int64']).columns

            # Remove columns that are all zeros, contain NaN, or have the same value
            cols_to_drop = []
            for col in numeric_cols:
                if df[col].isna().any() or df[col].nunique() == 1:
                    cols_to_drop.append(col)
                    # print(f"Column '{col}' in file '{filename}' has been removed (contains NaN or has the same value).")
            numeric_cols = numeric_cols.difference(cols_to_drop)

            # If no numeric columns to standardize, skip the file
            if len(numeric_cols) == 0:
                logging.warning(f"No numeric columns to standardize in file '{filename}'. Skipping this file.")
                continue

            # Standardize the remaining numeric columns
            if method == 'Zscore':
                standardized_values = (df[numeric_cols] - df[numeric_cols].mean()) / df[numeric_cols].std()
                standardized_df = pd.concat([sample_label_df, standardized_values], axis=1)
            elif method == 'MMS':
                scaler = MinMaxScaler() 
                standardized_values = scaler.fit_transform(df[numeric_cols])
                standardized_df = pd.concat([
                    sample_label_df,
                    pd.DataFrame(standardized_values, columns=numeric_cols)
                ], axis=1)
            elif method == 'RS':
                scaler = RobustScaler() 
                standardized_values = scaler.fit_transform(df[numeric_cols])
                standardized_df = pd.concat([
                    sample_label_df,
                    pd.DataFrame(standardized_values, columns=numeric_cols)
                ], axis=1)
            elif method == 'QT':
                scaler = QuantileTransformer() 
                standardized_values = scaler.fit_transform(df[numeric_cols])
                standardized_df = pd.concat([
                    sample_label_df,
                    pd.DataFrame(standardized_values, columns=numeric_cols)
                ], axis=1)

            # Apply number formatting to all numeric columns
            for col in numeric_cols:
                standardized_df[col] = standardized_df[col].apply(format_number)

            # Construct the output filename by adding the 'std_' prefix
            output_filename = f"{filename}"
            output_file_path = os.path.join(output_dir, output_filename)

            # Save the standardized DataFrame to the output directory
            standardized_df.to_csv(output_file_path, index=False)
            print(f"Standardized file saved: {output_file_path}")

    print("Standardization complete.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Standardization of CSV files.')
    parser.add_argument('--input_dir', type=str, required=True, help='Path to the input folder containing CSV files')
    parser.add_argument('--output_dir', type=str, required=True, help='Path to the output folder for standardized CSV files')
    parser.add_argument('--method', type=str, choices=['Zscore', 'MMS', 'RS', 'QT'], default='Zscore', help='Type of standardization method.')
    args = parser.parse_args()
    z_score_standardization(args.input_dir, args.output_dir, args.method)
