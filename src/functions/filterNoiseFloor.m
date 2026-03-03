function dives = filterNoiseFloor(maxDepth, azfpData, azfpAggregateData, Dives)
%FILTERNOISEFLOOR Moving average to determine noise floor for each bin (Steps 3 and 4) and noise removal

% According to Scott, each frequency should have its own noise floor due to
% frequency dependencies, differences in the conditions, etc
% We assume that the minimum Sv in each frequency is equivalent to the
% noise floor for that frequency

nFreq = numel(azfpData);
dbins = 1:maxDepth;
d_int = 10; % depth interval to average over

for i = 1:nFreq % for each frequency
    % preallocate enough space
    M(i).AvgSv = nan * ones(size(Dives, 2), length(dbins));
    for f = 1:size(Dives, 2) % for each dive
        % take the mean of d_int Sv values at a time and put them in the M
        % structure
        % any means that include a NaN are set to NaN
        % (need to include the NaNs here for depth window calculations
        % later)
        temp = movmean(Dives(f).P(i).msv, d_int, 'includenan', 'Endpoints', 'discard');
        M(i).AvgSv(f, 1:length(temp)) = temp;
    end
end

NoiseFloor = nan(nFreq,1);       
divenum   = nan(nFreq,1);       
DepthWindow = nan(nFreq, d_int + 1);
% Then find the minimum noise interval for each frequency
for i = 1:nFreq % for each frequency
    % find the minimum Sv value in the moving average and its index, not
    % counting any NaN values
    [N, index] = min(M(i).AvgSv, [], 'all', 'linear', 'omitnan');

    % raw minimum value
    NoiseFloor(i) = N;
    % dive number for each frequency where the minimum is located
    [D, J] = ind2sub(size(M(i).AvgSv), index);
    divenum(i) = D;

    while J >= 187
        J = J - 1;
    end

    % depth interval for each frequency where the minimum Sv is located
    DepthWindow(i, :) = dbins(J:J+d_int);
end

% Subtract the noise floor from the avg_sv structures
for i = 1:nFreq % for each frequency
    % subtract the frequency-dependent noise floor from avg_sv
    azfpAggregateData(i).avg_sv = azfpAggregateData(i).avg_sv - NoiseFloor(i);
end

for j = 1:size(Dives, 2) - 1 % for each dive
    for k = 1:nFreq % for each frequency
        % subtract the frequency-dependent noise floor
        Dives(j).P(k).sv = Dives(j).P(k).sv - NoiseFloor(k);
        % recalculate the median
        Dives(j).P(k).msv = median(Dives(j).P(k).sv, 1, 'omitnan');
    end
end

%% % add dive start and end time to each dive
% Find the indices of each dive in the AZFP data (part of Step 2 in the tutorial)
StartDive = find([1; diff(azfpData(1).Depth(:, 1)) < -10]);
for DD = 1:length(Dives)
    Dives(DD).starttime = azfpData(1).Date(StartDive(DD));
    Dives(DD).endtime = azfpData(1).Date(Dives(DD).Index(2));
end

%% db differencing
[~, n] = size(Dives);
for ii = 1:n 
    Dives(ii).P(1).Diff = real(10*log10(Dives(ii).P(2).sv)) - real(10*log10(Dives(ii).P(1).sv));
    Dives(ii).P(2).Diff = real(10*log10(Dives(ii).P(3).sv)) - real(10*log10(Dives(ii).P(2).sv));
    Dives(ii).P(3).Diff = real(10*log10(Dives(ii).P(4).sv)) - real(10*log10(Dives(ii).P(3).sv));
end

%% averaging the db differences per profile
for ii = 1:n 
    Dives(ii).P(1).mDiff = real(10*log10(median(10.^((Dives(ii).P(1).Diff) ./ 10), 'omitnan')));
    Dives(ii).P(2).mDiff = real(10*log10(median(10.^((Dives(ii).P(2).Diff) ./ 10), 'omitnan')));
    Dives(ii).P(3).mDiff = real(10*log10(median(10.^((Dives(ii).P(3).Diff) ./ 10), 'omitnan')));
end

dives = Dives;
end