function Output = filterAZFPData(Output,config)
%FILTERAZFPDATA Summary of this function goes here
%   Detailed explanation goes here

%% Standard Sphere Calibration application - Halyna used 2025 callibration data: 130 kHz (-0.72), 200 kHZ (-4.60), 455 kHZ (-4.63),769 kHZ (-0.20),

% Echosounder 59016 (Davies) calibration offset from standard sphere
% calibration conducted June 2024.  Davies and Mesquita have
% calibration files.
Output = standardSphereCallibration(config.calibrationOffsets, Output);

%% Trim Transmit Pulse and Near Field from Sv data (Step 5 in tutorial)
Output = filterNearField(Output);

%% Load in glider data
[gliderDepth, gliderTime, ~] = loadGliderData(config.gliderFileFullName, config.gliderVariableName);

%% glider data filtering for extra pings during inactivity
% Output = gliderDataFilter(gliderData, Output);

%% Time - align the glider and AZFP pressure data (Step 2 in tutorial)
[Output, fixedDepth] = timeAlignAZFPToGliderData(Output, gliderTime, gliderDepth);
fixedDepth(end) = [];
%% Far field noise cut off  (Steps 3 and 4 in the tutorial)
Output = filterFarField(Output, fixedDepth, config.farFieldCutOffRange);
%% Remove pings at the surface of the ocean because these often have bubbles in them
[Output, Dives] = filterSurfacePings(Output);

%% Histogram of Sv data for all frequencies
% getFigureAllFreqSvHystogram(Output)

%% Make a histogram of the seafloor data decibel strengths
% getFigureSeafloorDecibelStrength(Output)

%% Remove seafloor echoes (Step 1 in the tutorial)
Output = filterSeafloorEchoes(Output, Dives);
%% Remove surface noise
Output = filterSurfaceNoise(Output, Dives, config.surfaceNoiseDepthMin, config.surfaceNoiseDepthMax);

end
