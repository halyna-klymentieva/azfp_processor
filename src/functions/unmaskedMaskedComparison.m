function output = unmaskedMaskedComparison(dives, config)
% config.maskedOutAsZero (optional, default 0):
%   1 = samples outside the dB-difference window count as 0 (no copepods),
%       so medians/densities are per volume of water;
%   0 = they are NaN and left out (original behaviour), so medians are
%       among copepod-like samples only.
%   Samples where the dB difference cannot be calculated (no data at one
%   of the two frequencies) are NaN in both cases.
if nargin < 2
    config = struct();
end
maskedOutAsZero = isfield(config, 'maskedOutAsZero') && config.maskedOutAsZero;
%UNMASKEDMASKEDCOMPARISON Summary of this function goes here

%  June 2025: Added by Andréa to execute all processing steps in a single run

% First get the dB difference window
% Values below are for copepods between 1.27 and 2.99 mm in length, from
% Joe's spreadsheet

% windows are likely too small; play with these values until the matching
% matrix looks like the patches in the echogram
% 130-200 make 0-7 dB and see if that helps
% ignore 769 kHz for now

% Controlled parameter tuning based on MultiNet data
% Do correlations with windows in 200-455 kHz of 1 dB, 5 dB, 10 dB
% dB_Diff_Lower = [7.4, 13.7, 7.8];
dB_Diff_Lower = [5.1, 3.1, -0.8];
% dB_Diff_Upper = [7.5, 14.2, 8.8];
dB_Diff_Upper = [6.6, 8.9, 1.1];

% is the 455 kHz data "real" or just noise? Calibration issues? Offset or
% dynamic range

% pick a transect, look at the bottom value, see what the values are as a
% pseudo calibration
% if the bottom depth values are off, we will need to do a calibration
% correction; bottom is flat, broad, frequency-independent
% histogram of 1 m above bottom to 2 m below bottom (and right at the bottom)
% for each frequency, see how similar the values are (or how different)
% gives us insight into the sensitivity
% if the dynamic range of the different frequencies is off, this becomes
% trickier
%% Create the binary filter matrix

% frequency 1 > frequency 2 as a masking matrix, to start (and vice versa)

for i = 1:length(dives) % for each dive
    for j = 1:3 % for each frequency difference and dB window
        %if isfield(dives(i).P(j), "Diff") && ~isempty(dives(i).P(j).Diff)
        dB_Diff = dives(i).P(j).Diff; %dives.P.Diff is already in dB space
        dives(i).P(j).mask = (dB_Diff > dB_Diff_Lower(j)) & (dB_Diff < dB_Diff_Upper(j));
        %end
    end
end

% Then multiply masking matrix by Sv to get masked observed Sv

for i = 1:length(dives)
    dives(i).P(1).masked = dives(i).P(1).sv; % 130 kHz is not masked
end

for i = 1:length(dives)
    for j = 1:3 % for each frequency difference mask
        target = dives(i).P(j+1).sv;
        other = dives(i).P(j).sv;
        masked = target .* dives(i).P(j).mask;
        % dB difference unknown (one frequency has no data) -> NaN;
        % both frequencies 0 (nothing detected) -> known, 0
        known = ~isnan(dives(i).P(j).Diff) | (target == 0 & other == 0);
        masked(~known) = NaN;
        if ~maskedOutAsZero
            masked(masked == 0) = NaN; % original behaviour
        end
        dives(i).P(j+1).masked = masked;
    end
end
%% averaging the masked sv per profile (using median instead of mean)

[~, n] = size(dives);
for ii = 1:n
    dives(ii).P(1).mMasked = median(dives(ii).P(1).masked, 'omitnan');
    dives(ii).P(2).mMasked = median(dives(ii).P(2).masked, 'omitnan');
    dives(ii).P(3).mMasked = median(dives(ii).P(3).masked, 'omitnan');
    dives(ii).P(4).mMasked = median(dives(ii).P(4).masked, 'omitnan');
end

output = dives;
end
