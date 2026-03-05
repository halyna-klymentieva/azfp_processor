function [output,dive] = filterSurfacePings(azfpData)
%FILTERSURFACEPINGS Remove pings at the surface of the ocean 
% because these often have bubbles in them
arguments (Input)
    azfpData    
end

arguments (Output)
    output
    dive
end

fprintf('Removing pings at the surface...\n');
tic
cnt = length(azfpData(1).Date);

azfpData(1).Sv(cnt, :) = []; % May 2025:
azfpData(2).Sv(cnt, :) = []; % Andréa activated lines 272-277 to match # of elements across matrices
azfpData(3).Sv(cnt, :) = [];
azfpData(4).Sv(cnt, :) = [];

azfpData(1).Date(cnt, :) = [];

% Find the indices of each dive in the AZFP data (part of Step 2 in the tutorial)
StartDive = find([1; diff(azfpData(1).Depth(:, 1)) < -10]);

cc = 0; % number of distinct dives to go into the structure
for DD = 1:length(StartDive) - 1
    if (StartDive(DD+1) - 1 - StartDive(DD)) > 50
        cc = cc + 1;
        Dive(cc).Index = [StartDive(DD); StartDive(DD+1) - 1];
    end
end
Dive(cc+1).Index = [StartDive(DD+1); length(azfpData(1).Depth(:, 1))]; % manual entry of the last dive

for DD = 1:cc - 1
    Dive(DD).Index(2) = Dive(DD+1).Index(1);
end
Dive(end).Index(2) = length(azfpData(1).Depth);
toc

output = azfpData;
dive = Dive;
end