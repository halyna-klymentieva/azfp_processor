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
% Always work from the src folder: data/, output/, gliderData/ and the CTD
% file are found relative to it. (Run the whole script with Run / F5.)
cd(fileparts(mfilename('fullpath')));
addpath(genpath(pwd));
rehash path % pick up new or changed files (MATLAB can miss changes in OneDrive folders)
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
% Noise removal:
%  'range'    = range-dependent noise (De Robertis & Higginbottom 2007), per
%               sample, before the far-field cut (filterRangeNoise.m). Noise
%               level re-estimated every noiseBlockPings pings, smoothed over
%               noiseSmoothBlocks blocks. Samples < snrThreshold dB above the
%               noise -> 0. Needs the data-*.mat caches to be rebuilt.
%  'constant' = original: one noise floor per frequency and day, subtracted
%               from the 1 m bins (filterNoiseFloor.m).
config.noiseMethod = 'range';
config.noiseBlockPings = 20;
config.noiseSmoothBlocks = 11;
% Samples less than this many dB above the noise are treated as no detectable
% backscatter. 10 dB = De Robertis & Higginbottom 2007, McLennan et al. 2025.
% (With 'constant': [] = only samples below the noise floor.)
config.snrThreshold = 10;
% What samples below the noise floor (or SNR threshold) become: 0 or NaN
config.belowNoiseValue = 0;
% 1 = samples outside the copepod dB-difference windows count as 0
% (density per m^3 of water); 0 = left out as NaN (original behaviour)
config.maskedOutAsZero = 0;
% Krill check for the 455 and 769 kHz masks (unmaskedMaskedComparison.m): reject
% samples whose 200-130 kHz difference is <= krillCheck dB (Calanus 5.1-6.6 dB,
% krill ~0 dB; McQuinn et al. 2013). [] = off.
config.krillCheck = 3;
% Per-depth dive profiles (Sv, masked Sv, density): 'median' (original) or
% 'mean'. With the zeros above, medians are often 0; 'mean' = average per m^3.
config.profileStatistic = 'median';
% Near-surface exclusion: 1 m depth bins 1..N (m) are set to NaN in every dive
% (bubbles / glider at the surface; spikes fade by ~6 m). 5 = keep data deeper
% than 5.5 m. [] = keep everything.
config.surfaceExclusion = 5;
% Seafloor (filterSeafloorEchoes.m): onset = shallowest sample deeper than
% seafloorMinDepth with Sv >= seafloorThreshold (dB); seafloor depth per ping =
% shallowest onset of the 4 frequencies, moving median over seafloorSmoothPings
% pings; each frequency is cut at its own onset (if within seafloorMaxJump m of
% that depth, else at that depth) minus seafloorMargin m.
config.seafloorThreshold = -35;
config.seafloorMinDepth = 45;
config.seafloorMargin = 0.5;
config.seafloorSmoothPings = 101;
config.seafloorMaxJump = 3;
% Spike removal on the 1 m bins (filterSpikes.m): bins > spikeThreshold dB above
% the median of spikeWindowPings x spikeWindowDepth neighbours (median floored at
% spikeFloorDb), or with Sv > spikeCeilingDb, become NaN. spikeThreshold = [] = off.
config.spikeThreshold = 15;
config.spikeWindowPings = 11;
config.spikeWindowDepth = 3;
config.spikeFloorDb = -75;
config.spikeCeilingDb = -40;
config.maxDepth = 105;
% Colour scale (dB) for the Sv figures (median Sv and masked Sv)
config.svLimits = [-100, -50];
% Glider instrument schedule (CTD_processor): used to mark camera dives on figures
config.scheduleFile = fullfile(pwd, '..', '..', 'CTD_processor', 'CTD_processed', 'instrument_schedule_by_cast.csv');
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