%Version 1 created in 2020 by Kim Davies and Delphine Mossman

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update log:

% Jan 7, 2021 AZFP Noise Floor Estimate
% Scott Loranger

% 31 August 2021 Far Field Noise Floor removal code
% Modified Nov 19, 2021
% Delphine Mossman

% Version 2 created by Kim Davies
% Oct 4, 2024
% Code clean up and made some pieces more efficient for processing longer
% deployments.

% Halyna's version
% Feb 19, 2026
%Code clean up; include the code to save output from each day

%% Prepare your workspace and file directories
addpath(genpath(pwd));
clc
clear variables
close all

%% User-defined config variables
config.dates = ['25-07-24'; '25-07-25'; '25-07-26'; '25-07-27'; '25-07-28'; '25-07-29'];
% config.dates = '25-07-24';
config.xmlFileName = '25070417.XML';
config.gliderVariableName = 'cabot_20250723_213_delayed';
config.gliderFileName = 'cabot_20250723_213_delayed_0660_8583_ae5f.mat';
config.calibrationOffsets = [-0.72, -4.60, -4.63, -0.20];
% Far field per frequency cut-off range, see filterFarField function
config.farFieldCutOffRange = [30, 10, 10, 5];
% Surface-noise band removal (filterSurfaceNoise). Switched off: no noise
% band is visible in the 2025 data, and the old step only ran from dive 95.
% To use it again, set e.g. [10 14] (depth range in m, applied to all pings).
config.surfaceNoiseBand = [];
config.maxDepth = 105;
% Binned CTD file (CTD_processor) used for sound speed and absorption.
% '' = fixed values from the ASL toolbox (5 degC, salinity 33, 100 dbar)
config.ctdFile = fullfile(pwd, '..', '..', 'CTD_processor', 'CTD_processed', 'ctd_binned_1dbar_Jul24-29.csv');

%% Getting output of AZFP raw data processing
% before starting, you may want to increase the amount of memory that
% MATLAB can use.  Select Home - Preferences - General - Java Heap Memory
% and use the scale bar to increase memory.
%
% If output file alteady exist AZFP processing will not happen - data
% will be loaded from existing output
dates = string(config.dates);
config.saveWithDateName = 1;
for i = 1:length(dates) % for each dive
    [~, Dives, bottomDepth] = procesAZFPRawData1Day(config, dates{i});
    % draw figures
    config.dateOfData = dates{i};
    drawAndSaveFigures(Dives, bottomDepth, config)
end
clear i