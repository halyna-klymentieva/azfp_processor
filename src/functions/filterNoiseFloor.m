function dives = filterNoiseFloor(maxDepth, azfpData, azfpAggregateData, Dives, config)
% config.snrThreshold (dB, optional): after the noise floor is subtracted,
%   samples whose signal is less than this many dB above the noise floor are
%   set to 0 (no detectable backscatter). [] = only values below the noise
%   floor (negative) are set to 0. NaN (no data) stays NaN.
if nargin < 5
    config = struct();
end
snrThreshold = [];
if isfield(config, 'snrThreshold')
    snrThreshold = config.snrThreshold;
end
% config.belowNoiseValue: what weak samples become: 0 (no detectable
% backscatter, default) or NaN (treated as no data, left out of statistics)
belowNoiseValue = 0;
if isfield(config, 'belowNoiseValue') && ~isempty(config.belowNoiseValue)
    belowNoiseValue = config.belowNoiseValue;
end
belowNoiseLabel = 'set to 0';
if isnan(belowNoiseValue)
    belowNoiseLabel = 'set to NaN';
end
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

% config.noiseMethod = 'range': noise was already subtracted per sample with a
% range-dependent noise level and the SNR threshold (filterRangeNoise.m).
% Then no constant noise floor is subtracted here; 1 m bins whose samples were
% all below the SNR threshold (bin mean 0) become belowNoiseValue.
rangeMode = isfield(config, 'noiseMethod') && strcmpi(config.noiseMethod, 'range');
if rangeMode
    fprintf('Noise already removed per sample (range-dependent); constant noise floor not subtracted.\n');
    NoiseFloor(:) = 0;
end

% Subtract the noise floor from the avg_sv structures
for i = 1:nFreq % for each frequency
    % subtract the frequency-dependent noise floor from avg_sv
    azfpAggregateData(i).avg_sv = azfpAggregateData(i).avg_sv - NoiseFloor(i);
end

nZero = zeros(nFreq, 1); nValid = zeros(nFreq, 1);
for j = 1:size(Dives, 2) % for each dive (previously the last dive was skipped)
    for k = 1:nFreq % for each frequency
        % subtract the frequency-dependent noise floor (linear units)
        sv = Dives(j).P(k).sv - NoiseFloor(k);
        % below the noise floor, or below the SNR threshold = no detectable
        % backscatter -> 0 (previously negative values were kept and later
        % turned into fake dB values by real(10*log10(x)))
        if rangeMode
            weak = sv <= 0;
        elseif isempty(snrThreshold)
            weak = sv < 0;
        else
            weak = sv < NoiseFloor(k) * 10^(snrThreshold / 10);
        end
        sv(weak) = belowNoiseValue;
        Dives(j).P(k).sv = sv;
        nZero(k) = nZero(k) + sum(weak(:));
        nValid(k) = nValid(k) + sum(~isnan(sv(:))) + sum(weak(:) & isnan(belowNoiseValue));
        % recalculate the median
        Dives(j).P(k).msv = profileStatistic(sv, config);
    end
end
if rangeMode
    fprintf('1 m bins with no sample above the SNR threshold, %s: %s %% of bins\n', ...
        belowNoiseLabel, mat2str(round(100 * nZero' ./ max(nValid', 1))));
elseif isempty(snrThreshold)
    fprintf('Noise floor (dB): %s; %s (below noise): %s %% of samples\n', ...
        mat2str(round(10*log10(NoiseFloor'), 1)), belowNoiseLabel, mat2str(round(100 * nZero' ./ max(nValid', 1))));
else
    fprintf('Noise floor (dB): %s; %s (SNR < %g dB): %s %% of samples\n', ...
        mat2str(round(10*log10(NoiseFloor'), 1)), belowNoiseLabel, snrThreshold, mat2str(round(100 * nZero' ./ max(nValid', 1))));
end

%% % add dive start and end time to each dive
% Start and end times come from each dive's own first and last ping
% (Dives(DD).Index). Previously the start time was taken from a separate
% list of depth jumps that also contained the short segments (<= 50 pings)
% skipped when the dives were built, so one skipped segment would shift the
% start time of every later dive. (No segment was skipped in the 2025 data.)
for DD = 1:length(Dives)
    Dives(DD).starttime = azfpData(1).Date(Dives(DD).Index(1));
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