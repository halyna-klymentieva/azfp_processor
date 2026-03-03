function output = filterFarField(azfpData,fixedDepth, farFieldCutOffRange)
%FILTERFARFIELD Far field noise cut off  (Steps 3 and 4 in the tutorial)
%   Detailed explanation goes here
arguments (Input)
    azfpData
    fixedDepth
    farFieldCutOffRange
end

arguments (Output)
    output
end

% To start getting an idea of where the noise floor cutoff range is
% figure
% for ii=2:50:2000  % can change to plot more or less data.
%     scatter(Output(1).Range(2,:),Output(1).Sv(ii,:),'k') % frequency change line
%     % plots range vs frequency-dependent Sv
%     hold on
%     xlabel('range')
%     ylabel('Sv')
% end

% For Bay of Fundy, these values are: 130 kHz = 75 m, 200 kHz = 50 m, 455 kHz = 35 m, 769 kHz = 20 m
% these are our "eyeballed" values from the figures
% Remove these depths from the analysis

% For Baffin Bay mission; 15 m on the 130 kHz is all we are getting.

for i = 1:length(azfpData) % for each frequency
    J = find(azfpData(i).Range(1, :) >= farFieldCutOffRange(i)); % find the frequency-specific
    % far field data
    azfpData(i).Sv(:, J) = []; % eliminate the Sv in the far field
    azfpData(i).Range(:, J) = []; % eliminate the far field ranges
    % Redo the depth calculation to account for the far field being cut off
    azfpData(i).Depth = azfpData(i).Range(1, :) + fixedDepth;
    % azfpData(i).PingDepth = azfpData(i).Range(1,:) + azfpData(1).Depth;
end

output = azfpData;
end