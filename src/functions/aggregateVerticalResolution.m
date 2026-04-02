function [output, dives] = aggregateVerticalResolution(azfpData, maxDepth)
%AGGREGATEVERTICALRESOLUTION Average 10 cm vertical resolution into 1 m depth bins
%  to make the matrices smaller for better storage space
[m, ~] = size(azfpData(1).Depth);
P(1).avg_sv = NaN(m, maxDepth);
P(2).avg_sv = NaN(m, maxDepth);
P(3).avg_sv = NaN(m, maxDepth);
P(4).avg_sv = NaN(m, maxDepth);
fprintf('Aggregating data for 1m vertical resolution...\n');
tic
for ii = 1:length(azfpData) % for each frequency
    fprintf('Processing frequency #%d: *',ii);
    for pp = 1:size(azfpData(ii).Depth, 1) % for each ping
        Xw = 10.^(azfpData(ii).Sv(pp, :) ./ 10); % data
        id = round(azfpData(ii).Depth(pp, :)); % index
        idx = unique(id);
        mn = accumarray(id', Xw', [], @mean);
        mn(mn == 0) = [];
        P(ii).avg_sv(pp, idx) = mn;
        clear id idx mm Xw
        if (mod(pp,10000) == 0) 
            fprintf('*');
        end
    end
    P(ii).avg_sv = P(ii).avg_sv(:, 1:maxDepth);
    fprintf('\n');
end

% Find the indices of each dive in the AZFP data (part of Step 2 in the tutorial)
StartDive = find([1; diff(azfpData(1).Depth(:, 1)) < -10]);

cc = 0;
% First, define the dives based on frequency 1

%% Step 6: Organized Dive Structure (Stable Version)- Halyna's version
fprintf('Organize Dive Structure...\n');
Dive = struct('Index', {}, 'P', {});
for DD = 1:length(StartDive)
    % Determine the end index for this dive
    if DD == length(StartDive)
        end_idx = size(P(1).avg_sv, 1);
    else
        end_idx = StartDive(DD+1) - 1;
    end

    % Only process dives longer than 50 pings
    if (end_idx - StartDive(DD)) > 50
        cc = cc + 1;
        % Store the indices so we know exactly where this dive is
        Dive(cc).Index = [StartDive(DD), end_idx];

        % Now pull data for ALL frequencies (jj) into this dive
        for jj = 1:length(azfpData)
            Dive(cc).P(jj).sv = P(jj).avg_sv(StartDive(DD):end_idx, :);
            % Create the median profile (msv)
            Dive(cc).P(jj).msv = median(Dive(cc).P(jj).sv, 1, 'omitnan');
        end
    end
end
toc
fprintf('Successfully created %d dives with all frequencies aligned.\n', cc);

output = P;
dives = Dive;
end
