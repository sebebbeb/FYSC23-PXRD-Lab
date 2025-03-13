# FYSC23 Powder XRD Analysis

This project provides a Python-based tool for analyzing powder X-ray diffraction (XRD) data. It determines the lattice constant of a sample and automatically decides whether the sample has a face-centered cubic (fcc) or body-centered cubic (bcc) structure based on linear regression quality. Additionally, the project includes options for cleaning raw data files and running crystallite size calculations using the `powerxrd` module.

## Features

- **Lattice Constant Determination:**  
  Calculates the lattice constant from XRD data using linear regression on sin²θ versus theoretical Q values.
  
- **Crystal Structure Identification:**  
  Automatically selects between fcc and bcc structures based on regression quality.
  
- **Data Cleaning:**  
  Provides functions to clean raw XRD data (e.g., `Sample1.txt`, `Sample2.txt`, etc.) using a rolling median filter to remove outliers, saving the cleaned files as `filtered_sample<X>.xy`.
  
- **Crystallite Size Calculations:**  
  Includes size calculations routines for additional analysis steps such as back-subtraction, all peaks detection, and Scherrer peak analysis.

## Requirements

- Python 3.x (tested on 3.12.6)
- [NumPy](https://numpy.org/)
- [Pandas](https://pandas.pydata.org/)
- [Matplotlib](https://matplotlib.org/)
- [SciPy](https://www.scipy.org/)
- [powerxrd](https://github.com/andrewrgarcia/powerxrd/tree/main)

You can install the required packages using pip:

```bash
pip install -r requirements.txt
```

## Installation

```bash
git clone https://github.com/sebebbeb/FYSC23-PXRD-Lab.git
cd FYSC23-PXRD-Lab
```

## Usage

The main script is `main.py`. Assuming the format of the data doesn't change in later iterations of the course. I'd recommend running the complete anlysis command-line option:
### Complete Analysis
This cleans all raw data files found in the specified directory, saves them as `filtered_sample<X>.xy` and uses corresponding filtered file for the analysis. This analysis includes everything mentioned in the Features section above.
```bash
python main.py samples/Sample1.txt --clean --size
```

### Basic Analysis
Here an already cleaned data is required
```bash
python main.py samples/filtered_sample1.xy
```

### Clean Raw Data + Basic Analysis 
If you have raw data files(e.g. `Sample1.txt`, `Sample2.txt`, etc) and want to clean them before analysis (recommended), use the `--clean` flag. This will:
- Clean all raw data files found in the specified directory.
- Save cleaned files as `filtered_sample<X>.xy.`
- Use the corresponding filtered file for further analysis.
```bash
python main.py samples/Sample1.txt --clean
```

## Project Structure
- `main.py`:\
The main script that loads data, optionally cleans raw files, performs analysis (peak detection, regression, lattice constant calculation), and generates plots.

- `data_cleaner.py`:\
Contains functions for cleaning raw XRD data files. The cleaning functions process raw files (e.g., `Sample1.txt`) and save the output as filtered files (e.g., `filtered_sample1.xy`).

- `powderxrd_patch.py`:\
Contains additional functions for calculating the size of the crystallites and patching behaviors of the `powerxrd` library, including routines like `test_allpeaks`.

- `samples/`:\
A folder containing raw and filtered data files.

## Troubleshooting
- File Paths:\
Ensure that the file paths provided to the script are correct. Relative paths are interpreted from the script's working directory.

- Dependencies:\
Verify that all required packages are installed and that you are using a compatible Python version.

- Data Format:\
The script expects data files to be in a two-column, tab-separated format (2θ and intensity). Ensure your data conforms to this format.

- Peaks Not Being Recognized:\
This is something that might occur. I'd recommend changing the `height_frac` and `prominence` in the `detect_peaks` function in `main.py`.\
If your problem lies in the Gaussian fits and Scherrer calculations take a look in `src/powderxrd_patch.py` and add the peaks manually in the `test_allpeaks` function. one peak around 2θ=35(deg) is already added manually. 

- Dead Pixels/Spots Being Counted:\
If this is the case the data has somewhat changed since the code was written. The window size in the `clean_data` function of `src/data_cleaner.py` might need to be changed then. Furthermore, if this is the case, then you might need to change/delete the manually added peak in `test_allpeaks` as mentioned above. 

## Additional Notes
I've tried to make this program as modifiable as possible, I encourage you to modify all you see fit. 


