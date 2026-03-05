function output = filterNearField(azfpData)
%FILTERNEARFIELD Trim Transmit Pulse and Near Field from Sv data (Step 5 in tutorial)
% Using 1 m as the calculated Rb for the highest frequency (769 kHz) is ~2 m;
% therefore this should eliminate the near-field data from all four
% frequencies
arguments (Input)
    azfpData
end

arguments (Output)
    output
end

fprintf('Filtering near field 2m...\n');
tic
I = find(azfpData(1).Range(1, :) <= 2);
for i = 1:length(azfpData)
    azfpData(i).Sv(:, I) = [];
    azfpData(i).Range(:, I) = [];
end
toc


output = azfpData;
end