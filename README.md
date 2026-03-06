# AZFP MATLAB Processor

A MATLAB toolbox for processing and analyzing data from the ASL Environmental Sciences **Acoustic Zooplankton Fish Profiler (AZFP)**. This repository provides scripts to process raw acoustic data, integrate it with glider telemetry, and generate standardized outputs for analysis.

## Project Overview

This project provides a comprehensive suite of scripts and functions for:
- Initial processing of raw AZFP data (`.01*` files).
- Data averaging (time and range bins) and filtering (noise floor, seafloor echoes).
- Frequency differencing for target identification.
- Integration with external datasets (e.g., glider telemetry, multinet samples).
- Advanced analysis such as dive averaging and seafloor echo processing.

The core processing logic is based on the **AzfpMatlabToolbox_v18** by ASL Environmental Sciences Inc.

## Requirements

- **MATLAB:** Tested with version R2024a (TODO: Verify minimum version required).
- **Toolbox Dependency:** This project requires `AzfpMatlabToolbox_v18` (included in `src/azfp`).
- **Memory:** Large datasets may require increasing MATLAB's Java Heap Memory (Home -> Preferences -> General -> Java Heap Memory).

## Installation and Setup

1.  **Clone the Repository:**
    ```bash
    git clone https://github.com/halyna-klymentieva/azfp_processor.git
    cd azfp_processor/src
    ```

2.  **Data Preparation:** 
    - Place raw `.01*` and `.XML` data files in the `data/` folder.
    - Place glider telemetry data (`.mat` format) in the `gliderData/` folder.
    - Ensure an `output/` directory exists at the root.

## Usage Guide

### Basic Workflow

1.  **Initialize Parameters:** Open `src/AZFP_Initial_Processing.m` and configure the `User-defined config variables` section:
    - `config.dates`: Dates to process (e.g., `['25-07-24'; '25-07-25']`).
    - `config.xmlFileName`: The `.XML` configuration file name (e.g., `'25072317.XML'`).
    - `config.gliderFileName`: The glider data file name (e.g., `'cabot_20250723_213_delayed_0660_8583_ae5f.mat'`).
    - `config.calibrationOffsets`: Frequency-specific calibration adjustments.
    - `config.farFieldCutOffRange`: Far field per frequency cut-off range.
    - `config.dataCutOffDepth`: Data cut-off depth.
    - `config.surfaceNoiseDepthMin`/`Max`: Surface noise depth range.
    - `config.maxDepth`: Maximum depth for processing.
    
2.  **Process Data:** Run `src/AZFP_Initial_Processing.m` in MATLAB. This script:
    - Loads raw data and configuration.
    - Applies calibration and filtering (Near-field, Far-field, Seafloor).
    - Generates processed `.mat` files in the `output/` directory for each day.
    - Produces visualization plots (Echograms, Sv Histograms, etc.).

3.  **Merge Results:** Run `src/Merge_Dives_Data.m` to aggregate multiple days of processed data into a single `0-merged-dives.mat` file and generate combined plots.

## Key Scripts and Components

### Core Processing
- **`src/AZFP_Initial_Processing.m`**: Primary entry point for daily processing and visualization.
- **`src/Merge_Dives_Data.m`**: Merges multiple daily outputs into a single dataset.
- **`src/azfp/`**: Contains the core ASL toolbox functions (`LoadAZFP.m`, `ProcessAZFP.m`, `PlotAZFP.m`, `LoadAZFPxml.m`).
- **`src/functions/procesAZFPRawData1Day.m`**: Orchestrates the loading and processing of one day of data.

### Analysis & Utilities (`src/functions/`)
- **Filtering:** `filterAZFPData.m`, `filterNoiseFloor.m`, `filterSeafloorEchoes.m`, `filterFarField.m`, `filterNearField.m`, `filterSurfaceNoise.m`, `filterSurfacePings.m`.
- **Plotting:** `drawAndSaveFigures.m`, `getFigureDayVNight1.m`, `getFigureDayVNight2.m`, `getFigureMaskedSvAllFrqs.m`, `getFigureAllFreqSvHystogram.m`, `getFigureDBDiffAllFrqs.m`, `getFigureMedianSvAllFrqs.m`, `getFigureSeafloorDecibelStrength.m`.
- **Aggregation:** `mergeDives.m`, `aggregateDivesData.m`, `aggregateVerticalResolution.m`.
- **Integration:** `timeAlignAZFPToGliderData.m`, `loadGliderData.m`, `gliderDataFilter.m`.
- **Other Utilities:** `calculateDiveBottomDepth.m`, `calcNumericalDensity.m`, `saveAZFPData.m`, `saveDiveData.m`, `standardSphereCallibration.m`.

## Project Structure

```text
azfp_processor/
├── data/               # Input: Raw AZFP data (.01*) and XML files
├── gliderData/         # Input: Auxiliary glider telemetry data (.mat)
├── output/             # Output: Processed .mat files and generated plots
├── src/                # Source code
│   ├── azfp/           # Core ASL AZFP Matlab Toolbox (v18)
│   ├── functions/      # Utility, filtering, and analysis functions
│   ├── AZFP_Initial_Processing.m # Main entry point script for daily processing
│   └── Merge_Dives_Data.m        # Merges daily output files
├── LICENSE             # MIT License
└── README.md           # Project documentation
```

## Environment Variables

(None required. Configuration is handled within the MATLAB scripts' `User-defined config variables` sections.)

## Tests

(TODO: No formal test suite exists yet. Verification is currently performed by inspecting output plots and `.mat` files in the `output/` directory.)

## Credits

- **Original Toolbox:** Dave Billenness, ASL Environmental Sciences Inc.
- **Modifications & Analysis Scripts:** Kim Davies, Andrea Mesquita, Scott Loranger, Delphine Mossman, Halyna Klymentieva, and others at UNB.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
