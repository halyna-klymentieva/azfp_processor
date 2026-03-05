function [output, fixedDepth] = timeAlignAZFPToGliderData(azfpData, gliderTime, gliderDepth)
%TIMEALINEAZFPTOGLIDERDATA Time - align the glider and AZFP pressure data (Step 2 in tutorial)
%   Detailed explanation goes here
arguments (Input)
    azfpData
    gliderTime
    gliderDepth
end

arguments (Output)
    output
    fixedDepth
end
fprintf('Aligning the glider and AZFP pressure data...\n');
tic
% for each date in the echosounder Output file, create a timeindex entry equal to
% the index # of where the minimum difference between each recorded echosounder time
% stamp and every non-NaN glider time stamp is; this is time-aligning the glider and
% echosounder data and assumes no clock drift between glider and AZFP
timeindex = zeros(1, length(azfpData(1).Date));
for ii = 1:length(azfpData(1).Date)
    [~, timeindex(ii)] = min(abs(azfpData(1).Date(ii)-gliderTime));
end
fixedDepth = gliderDepth(timeindex);

% Create a new variable in the echosounder Output structure called Depth,
% and make it equal to the echosounder Range (i.e. transducer ping depth) plus
% the glider depth at that time
% Then make three more Depths of equal value, so you have one per frequency
azfpData(1).Depth = azfpData(1).Range(1, :) + fixedDepth;
azfpData(2).Depth = azfpData(1).Depth(:, 1:size(azfpData(2).Sv, 2));
azfpData(3).Depth = azfpData(1).Depth(:, 1:size(azfpData(3).Sv, 2));
azfpData(4).Depth = azfpData(1).Depth(:, 1:size(azfpData(4).Sv, 2));
toc

output = azfpData;
end
