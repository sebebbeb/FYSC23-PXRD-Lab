#!/usr/bin/env python3
"""
FYSC23 Powder XRD Analysis Script

This script processes powder X-ray diffraction (XRD) data to determine the lattice 
constant and automatically decide if the sample has a face-centered cubic (fcc) or 
body-centered cubic (bcc) structure based on linear regression quality. It also provides
routines to calculate the crystallite size of the sample using the powerxrd module.

Features:
  - Loads diffraction data (expected as two columns: 2θ in degrees and intensity, tab-separated).
  - Optionally cleans raw data files using a rolling median filter to remove outliers.
  - Automatically selects the filtered data file corresponding to the provided raw file.
  - Performs peak detection and linear regression of sin²θ versus theoretical Q values.
  - Calculates the lattice constant and identifies the best crystal structure (fcc or bcc).
  - Includes optional test routines for additional analyses (e.g., back-subtraction, Scherrer peak analysis).

Usage:
  To analyze a pre-cleaned data file:
      python main.py samples/filtered_sample1.xy

  To clean raw data files (e.g., samples/Sample1.txt, samples/Sample2.txt, etc.) and then 
  analyze the filtered data corresponding to the given raw file:
      python main.py samples/Sample1.txt --clean

  To also run size calculations of the sample using the powerxrd module:
      python main.py samples/Sample1.txt --clean --size

Arguments:
  datafile      Path to the input data file (raw or filtered).
  --wavelength  X-ray wavelength in Å (default: 0.7107 for Mo Kα). 
  ** If you need another wavelength these need to be manually changed in powderxrd_patch.py in src **
  --clean       Flag to clean raw data files and use the filtered versions for analysis.
  --size        Flag to run additional routines for size calculations using the powerxrd module.
"""

import argparse
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.stats import linregress
import os
import src.powderxrd_patch as powderxrd_patch
from src.data_cleaner import clean_all_raw_data

def load_data(filename):
    """
    Load diffraction data from a file.
    Expected file format: two columns (2θ in degrees, intensity) separated by tabs.
    """
    data = np.loadtxt(filename, delimiter='\t')
    two_theta = data[:, 0]
    intensity = data[:, 1]
    return two_theta, intensity

def detect_peaks(two_theta, intensity, height_frac=0.001, prominence=1):
    """
    Detect peaks in the diffraction pattern.
    Returns:
      peak_angles: detected 2θ positions,
      peak_intensities: intensities at those positions.
    """
    height_threshold = np.max(intensity) * height_frac
    peaks, _ = find_peaks(intensity, height=height_threshold, prominence=prominence)
    return two_theta[peaks], intensity[peaks]

def get_allowed_reflections(structure):
    """
    Return allowed reflections for cubic crystals as a list of tuples (Q, Miller label).
    For fcc:
      (111): Q=3, (200): Q=4, (220): Q=8, (311): Q=11, (222): Q=12, ...
    For bcc:
      (110): Q=2, (200): Q=4, (211): Q=6, (220): Q=8, (310): Q=10, ...
    """
    structure = structure.lower()
    if structure == 'fcc':
        reflections = [(3, '(111)'), (4, '(200)'), (8, '(220)'), (11, '(311)'), (12, '(222)')]
    elif structure == 'bcc':
        reflections = [(2, '(110)'), (4, '(200)'), (6, '(211)'), (8, '(220)'), (10, '(310)')]
    else:
        raise ValueError("Structure type must be 'fcc' or 'bcc'.")
    return reflections

def perform_regression(sin2_theta, reflections, wavelength):
    """
    Performs a linear regression of sin²θ versus theoretical Q values.
    Returns:
      slope, intercept, r_squared, and the calculated lattice constant.
    """
    n_reflections = min(len(reflections), len(sin2_theta))
    Q_theoretical = np.array([reflections[i][0] for i in range(n_reflections)])
    sin2_theta_used = sin2_theta[:n_reflections]
    regression = linregress(Q_theoretical, sin2_theta_used)
    slope = regression.slope
    intercept = regression.intercept
    r_squared = regression.rvalue**2
    # Calculate lattice constant using a = λ/(2*sqrt(slope))
    a = wavelength / (2 * np.sqrt(slope))
    return slope, intercept, r_squared, a, Q_theoretical, sin2_theta_used

def main():
    os.makedirs("figs", exist_ok=True)
    parser = argparse.ArgumentParser(
        description="Determine lattice constant and crystal structure from XRD data."
    )
    parser.add_argument(
        'datafile',
        help="Raw data file (e.g., samples/Sample1.txt)"
    )
    parser.add_argument(
        '--wavelength',
        type=float,
        default=0.7107,
        help="X-ray wavelength in Å (default: 0.7107 for Mo Kα)"
    )
    parser.add_argument(
        '--clean',
        action='store_true',
        help="Clean raw data files and use filtered data for analysis."
    )
    parser.add_argument(
        '--size',
        action='store_true',
        help="Run size calculation functions after analysis."
    )
    args = parser.parse_args()

    # If --clean flag is provided, clean all raw files in the directory.
    if args.clean:
        clean_all_raw_data(args.datafile)
        # Update the datafile to use the corresponding filtered file.
        base = os.path.basename(args.datafile)
        match = re.search(r'Sample(\d+)', base, re.IGNORECASE)
        if match:
            sample_num = match.group(1)
            filtered_file = os.path.join(os.path.dirname(args.datafile), f"filtered_sample{sample_num}.xy")
            args.datafile = filtered_file
            print(f"Using filtered data file: {args.datafile}")
        else:
            print("Warning: Could not determine sample number from filename; using raw file.")

    # Proceed with analysis using args.datafile.
    two_theta, intensity = load_data(args.datafile)
    peak_angles, _ = detect_peaks(two_theta, intensity)
    peak_angles = np.sort(peak_angles)
    theta_deg = peak_angles / 2.0
    theta_rad = np.deg2rad(theta_deg)
    sin2_theta = np.sin(theta_rad)**2

    structures = ['fcc', 'bcc']
    results = {}
    for structure in structures:
        reflections = get_allowed_reflections(structure)
        slope, intercept, r_squared, a, Q_theoretical, sin2_used = perform_regression(
            sin2_theta, reflections, args.wavelength
        )
        results[structure] = {
            'slope': slope,
            'intercept': intercept,
            'r_squared': r_squared,
            'lattice_constant': a,
            'Q_theoretical': Q_theoretical,
            'sin2_used': sin2_used,
            'miller_labels': [r[1] for r in reflections[:len(Q_theoretical)]]
        }
        print(f"Structure: {structure.upper()}")
        print(f"  Slope       = {slope:.5e}")
        print(f"  Intercept   = {intercept:.5e}")
        print(f"  R-squared   = {r_squared:.5f}")
        print(f"  Lattice constant a = {a:.5f} Å\n")

    best_structure = max(results, key=lambda s: results[s]['r_squared'])
    best_result = results[best_structure]
    print(f"Automatically determined structure: {best_structure.upper()}")
    print(f"Calculated lattice constant a = {best_result['lattice_constant']:.5f} Å")

    # Plot the linear fit for the best structure.
    plt.figure()
    plt.plot(best_result['Q_theoretical'], best_result['sin2_used'], 'o', label='Data')
    Q_fit = np.linspace(min(best_result['Q_theoretical']), max(best_result['Q_theoretical']), 100)
    sin2_fit = best_result['slope'] * Q_fit + best_result['intercept']
    plt.plot(Q_fit, sin2_fit, '-', label=f'Linear fit (r² = {best_result["r_squared"]:.5f})')
    plt.xlabel('Q=(h²+k²+l²)')
    plt.ylabel('sin²θ')
    plt.title(f'Linear Regression for Sample {sample_num}, {best_structure.upper()} Structure')
    plt.grid()
    plt.legend()
    plt.savefig(f"figs/lin_reg_sample{sample_num}.pdf")
    plt.show()

    # Plot the diffraction pattern with annotated peaks.
    plt.figure()
    plt.plot(two_theta, intensity, label=f'Sample {sample_num} data')
    n_annotate = len(best_result['miller_labels'])
    for angle, label in zip(peak_angles[:n_annotate], best_result['miller_labels']):
        intensity_at_peak = np.interp(angle, two_theta, intensity)
        plt.annotate(label, xy=(angle, intensity_at_peak), xytext=(0, 10),
                     textcoords='offset points', ha='center', fontsize=10, color='red')
    plt.xlabel('2 $\\theta$ (deg)')
    plt.ylabel('Intensity (a.u.)')
    bottom, _ = plt.ylim()
    plt.ylim(bottom, 1.1 * max(intensity))
    plt.title(f'Diffraction Pattern with Indexed Peaks for Sample {sample_num}')
    plt.grid()
    plt.legend()
    plt.savefig(f"figs/diff_pattern_sample{sample_num}.pdf")
    plt.show()

    # Run size calculations if --size flag is provided.
    if args.size:
        powderxrd_patch.backsub_multiplt() # This function does not actually provide any calculations, it's just nice to have/see
        powderxrd_patch.all_peaks(args.datafile)

if __name__ == "__main__":
    main()


