%% Day vs night median density profile (one frequency), non-camera dives
% One graph: median numerical density against depth for day and night,
% with shading showing the variation between days.
%   line    = median of the daily median profiles
%   shading = range (min-max) of the daily median profiles
% For each day (or night) the dives are combined into one median profile
% first, so every day counts the same however many dives it has.
%
% Uses output/0-merged-dives.mat (run AZFP_Initial_Processing and
% Merge_Dives_Data first) and the camera schedule from CTD_processor;
% camera dives are left out.
% Day = 07:00-20:59, night = 23:00-04:59 local time (Dives(i).tod, set in
% plotDayNight.m); twilight dives are not used. A night is labelled with the
% date it starts (night of 24 Jul = 24 Jul 23:00 to 25 Jul 04:59).
%
% Output: output/figure-density-<freq>-day-night-median.png

cd(fileparts(mfilename('fullpath')));
addpath(genpath(pwd));
rehash path
clc
clear variables
close all

%% Settings
config.outputFolder = fullfile(pwd, '..', 'output');
config.scheduleFile = fullfile(pwd, '..', '..', 'CTD_processor', 'CTD_processed', 'instrument_schedule_by_cast.csv');
config.maxDepth     = 105;

P_index     = 3;          % Dives(i).P(j): 2 = 200 kHz, 3 = 455 kHz, 4 = 769 kHz
freqLabel   = '455 kHz';
depthRange  = [5, 80];    % m, depths shown (top 5 m left out: surface spikes)
minDives    = 3;          % dives needed for a daily profile, and at each depth of it
minDays     = 2;          % daily profiles needed at a depth to draw it
xMax        = [];         % ind./m^3, x-axis limit ([] = automatic; e.g. 12000)

%% Load merged dives, keep non-camera dives
S = load(fullfile(config.outputFolder, '0-merged-dives.mat'), 'Dives');
Dives = S.Dives;
cam = getCameraDives(Dives, config);
if isempty(cam)
    error('Camera schedule not found: %s', config.scheduleFile);
end
keep = find(~cam(:)');

tod       = string({Dives(keep).tod});
tLocal    = [Dives(keep).starttime] - 3/24;       % datenum, local time (UTC-3)
groupDate = floor(tLocal - 6/24);                 % night after midnight -> date it started
depths    = (1:config.maxDepth)';
inRange   = depths >= depthRange(1) & depths <= depthRange(2);

% density matrix: depth x dive
N = NaN(config.maxDepth, numel(keep));
for k = 1:numel(keep)
    p = Dives(keep(k)).P(P_index).mN(:);
    n = min(numel(p), config.maxDepth);
    N(1:n, k) = p(1:n);
end

%% Daily median profiles and the summary for day and night
periods = ["Day", "Night"];
cols    = {[0.93, 0.55, 0.00], [0.10, 0.20, 0.55]};
names   = {'day', 'night'};
dayWord = {'days', 'nights'};

fig = figure(1);
try, theme(fig, 'light'); end                     % white background (R2025a+ dark theme)
set(fig, 'Color', 'w', 'InvertHardcopy', 'off', 'Position', [100, 100, 600, 800])
hold on
h = gobjects(1, 2);
fprintf('%s, non-camera dives (median of daily medians):\n', freqLabel);
for p = 1:2
    sel   = tod == periods(p);
    dates = unique(groupDate(sel));
    daily = NaN(config.maxDepth, 0);
    nDives = 0;
    for d = 1:numel(dates)
        s = sel & groupDate == dates(d);
        if sum(s) < minDives, continue; end
        dm = median(N(:, s), 2, 'omitnan');
        dm(sum(~isnan(N(:, s)), 2) < minDives) = NaN;   % per depth: need minDives dives with data
        daily(:, end+1) = dm; %#ok<SAGROW>
        nDives = nDives + sum(s);
    end
    md = median(daily, 2, 'omitnan');
    lo = min(daily, [], 2);
    hi = max(daily, [], 2);
    ok = inRange & sum(~isnan(daily), 2) >= minDays;

    fill([lo(ok); flipud(hi(ok))], [depths(ok); flipud(depths(ok))], cols{p}, ...
        'FaceAlpha', 0.22, 'EdgeColor', 'none', 'HandleVisibility', 'off')
    h(p) = plot(md(ok), depths(ok), '-', 'Color', cols{p}, 'LineWidth', 2.2, ...
        'DisplayName', sprintf('%s (%d %s, %d dives)', names{p}, size(daily, 2), dayWord{p}, nDives));

    [pk, ipk] = max(md .* ok);
    fprintf('  %-5s: %d %s, %d dives; peak %.0f ind./m^3 at %d m; mean over %d-%d m %.0f ind./m^3\n', ...
        names{p}, size(daily, 2), dayWord{p}, nDives, pk, depths(ipk), depthRange, mean(md(ok), 'omitnan'));
end
hold off

set(gca, 'YDir', 'reverse', 'Box', 'on')
ylim(depthRange)
if ~isempty(xMax), xlim([0, xMax]); else, xlim([0, inf]); end
grid on
xlabel('Numerical density (ind. m^{-3})')
ylabel('Depth (m)')
title(sprintf('%s, non-camera dives, 24-29 July 2025', freqLabel))
legend(h, 'Location', 'east')
text(0.98, 0.02, {'line = median of daily median profiles', 'shading = range (min-max) between days'}, ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
    'FontSize', 8, 'Color', [0.3, 0.3, 0.3])

drawnow
outName = sprintf('figure-density-%s-day-night-median.png', strrep(freqLabel, ' kHz', ''));
print(fig, '-dpng', fullfile(config.outputFolder, outName), '-r300')
fprintf('Saved %s\n', outName);
