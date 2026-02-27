# AZFP MATLAB Toolbox - UNB

A MATLAB toolbox for processing and analyzing data from the ASL Environmental Sciences **Acoustic Zooplankton Fish Profiler (AZFP)**.

## Project Overview

This project provides a comprehensive suite of scripts and functions for:
- Initial processing of raw AZFP data (`.01*` files).
- Data averaging (time and range bins).
- Frequency differencing for target identification.
- Integration with external datasets (e.g., glider telemetry, multinet samples).
- Advanced analysis such as dive averaging and seafloor echo processing.

The core processing logic is based on the **AzfpMatlabToolbox_v18** by ASL Environmental Sciences Inc.

## Installation and Setup

1. **Add to MATLAB Path:** Ensure all folders within `src/` (and its subdirectories) are added to your MATLAB path. Many scripts use `addpath(genpath(...))` to include dependencies.
2. **Toolbox Dependency:** This project requires the `AzfpMatlabToolbox_v18`.
3. **Seawater Library:** Includes `seawater_ver3_3.1` for calculating acoustic parameters like sound speed and absorption.

## Key Components

### Core Processing
- **`AZFP_Initial_Processing.m`**: The primary entry point for raw data processing. Loads `.01*` and `.XML` files, applies calibration parameters, and outputs processed `.mat` files. (Refactored to Version 2 by Halyna).
- **`ProcessAZFP.m`**: High-level wrapper function to load and process AZFP files based on a `Parameters` structure.
- **`LoadAZFP.m`**: Low-level function to read binary AZFP data.
- **`ParametersAZFP.m`**: Defines default processing and plotting parameters.

### Analysis & Utilities
- **`AZFP_Differencing.m`**: Performs frequency differencing (e.g., 200kHz minus 130kHz) in linear space and averages results into depth bins.
- **`AZFP_Integration.m`**: Aligns AZFP data with external position data (e.g., from a glider) and multinet deployment logs.
- **`AZFP_Dive_Averaging.m`**: Calculates median acoustic returns across different glider dives.
- **`PlotAZFP.m`**: Generates echograms and plots for Sv, TS, Counts, or Temperature/Tilts.
- **`find_bottom.m` / `removeBottom.m`**: Functions for detecting and removing seafloor echoes.

### Supporting Functions
- **`readULS6.m`**: Specialized reader for ULS6 instruments.
- **`cmocean.m`**: Perceptually uniform colormaps for oceanography.
- **`ddm2dd.m`**: Converts degrees and decimal minutes to decimal degrees.

## Usage Guide

### Basic Workflow
1.  **Initialize Parameters:** Open `AZFP_Initial_Processing.m` and configure your paths and processing parameters (e.g., `Bins2Avg`, `Time2Avg`, `Salinity`).
2.  **Process Data:** Run `AZFP_Initial_Processing.m`. It will prompt you to select the data folder and the XML configuration file.
3.  **Analyze & Plot:** Use `PlotAZFP.m` to visualize results or `AZFP_Differencing.m` for multi-frequency analysis.

*Note: Ensure you have run the initial processing before attempting to run differencing or integration scripts, as they rely on the processed `.mat` output.*

## Credits

- **Original Toolbox:** Written by Dave Billenness, ASL Environmental Sciences Inc. (dbillenness@aslenv.com).
- **Modifications & Analysis Scripts:** Scott Loranger, Delphine Mossman, and others at UNB.
- **Colormaps:** `cmocean` by Chad Greene.
