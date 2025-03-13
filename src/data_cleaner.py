import pandas as pd
import numpy as np
import glob
import os
import re

def clean_data(raw_file):
    """
    Cleans a single raw data file and saves the filtered data.
    For example, converts samples/Sample1.txt to samples/filtered_sample1.xy.
    """
    df = pd.read_csv(raw_file, sep="\t", header=None)
    df = df[df[0] > 10]  # Filter condition on the first column
    window_size = 5

    # Compute rolling median and MAD for the intensity (column 1)
    df["rolling_median"] = df[1].rolling(window=window_size, center=True).median()
    df["rolling_mad"] = df[1].rolling(window=window_size, center=True).apply(
        lambda x: np.median(np.abs(x - np.median(x))), raw=True
    )
    df["difference"] = np.abs(df[1] - df["rolling_median"])
    threshold_multiplier = 3
    epsilon = 1e-6
    df_filtered = df[df["difference"] <= threshold_multiplier * (df["rolling_mad"] + epsilon)].copy()

    # Keep only the original two columns (2θ and intensity)
    df_filtered = df_filtered[[0, 1]]
    
    # Create a filtered filename using the sample number from the raw file name.
    base = os.path.basename(raw_file)
    match = re.search(r'Sample(\d+)', base, re.IGNORECASE)
    sample_num = match.group(1) if match else "unknown"
    filtered_filename = os.path.join(os.path.dirname(raw_file), f"filtered_sample{sample_num}.xy")
    
    # Save the cleaned data.
    df_filtered.to_csv(filtered_filename, sep='\t', index=False, header=False)
    print(f"Saved filtered data to {filtered_filename}")

def clean_all_raw_data(example_file):
    """
    Finds all raw data files in the same directory as example_file that match
    'Sample*.txt' and cleans each one.
    """
    directory = os.path.dirname(example_file)
    raw_files = glob.glob(os.path.join(directory, "Sample*.txt"))
    for rf in raw_files:
        clean_data(rf)
