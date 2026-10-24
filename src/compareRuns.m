%% Compare two processing runs side by side
% Draws, for one day, the dive-median Sv and the masked Sv of two runs next
% to each other on the SAME fixed colour scale, plus the difference B - A
% (dB) per depth bin. Figures go to output/compare/.
% Needs only the dives-<date>.mat file of each run (no reprocessing).
%
% Run with Run / F5.

cd(fileparts(mfilename('fullpath')));
addpath(genpath(pwd));
rehash path % pick up new or changed files (MATLAB can miss changes in OneDrive folders)
clear variables
close all

%% Settings: which runs to compare
dateOfData = '25-07-24';

runA.folder = fullfile('..', 'output_1');            % original code
runA.label  = 'A: original (5 degC, S 33)';

runB.folder = fullfile('..', 'output');              % current run
runB.label  = 'B: CTD + NaN + median';

% other runs you have, e.g.:
%   fullfile('..', 'output_noNa_median')  (zeros, SNR 10 dB, median)
%   fullfile('..', 'output_noNA_mean')    (zeros, SNR 10 dB, mean)

maxDepth = 105;
svLimits = [-100, -50];    % dB, same as figures 1 and 3 (config.svLimits)
diffLimits = [-6, 6];      % dB, colour range of the difference panels
outFolder = fullfile('..', 'output', 'compare');

%% Load both runs
A = load(fullfile(runA.folder, "dives-" + dateOfData + ".mat"), 'Dives');
B = load(fullfile(runB.folder, "dives-" + dateOfData + ".mat"), 'Dives');
cA = getChartsData(A.Dives);
cB = getChartsData(B.Dives);
fprintf('%s: %d dives;  %s: %d dives\n', runA.label, numel(A.Dives), runB.label, numel(B.Dives));

if ~isfolder(outFolder)
    mkdir(outFolder);
end

freqs = {'130', '200', '455', '769'};
kinds = {'Sv', 'masked'};
kindTitles = {'Dive-median Sv', 'Masked Sv'};
diffMap = blueWhiteRed(256);

for k = 1:numel(kinds)
    fig = figure('Position', [50, 50, 1500, 1050], 'Color', 'w');
    tl = tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('%s, %s', kindTitles{k}, dateOfData));
    fprintf('\n%s (B - A, dB, over bins with data in both runs):\n', kindTitles{k});

    for i = 1:4
        fA = cA.([kinds{k} '_' freqs{i}]);
        fB = cB.([kinds{k} '_' freqs{i}]);

        ax = nexttile;
        drawPanel(ax, fA, svLimits, jet(256), maxDepth);
        title(ax, sprintf('%s kHz - %s', freqs{i}, runA.label));
        cb = colorbar(ax); ylabel(cb, 'Sv (dB re 1 m^{-1})');

        ax = nexttile;
        drawPanel(ax, fB, svLimits, jet(256), maxDepth);
        title(ax, sprintf('%s kHz - %s', freqs{i}, runB.label));
        cb = colorbar(ax); ylabel(cb, 'Sv (dB re 1 m^{-1})');

        ax = nexttile;
        if isequal(size(fA), size(fB))
            d = fB - fA;
            d(~isfinite(d)) = NaN;
            drawPanel(ax, d, diffLimits, diffMap, maxDepth);
            title(ax, sprintf('%s kHz - difference B - A', freqs{i}));
            cb = colorbar(ax); ylabel(cb, 'dB');

            both = ~isnan(fA) & ~isnan(fB) & isfinite(fA) & isfinite(fB);
            onlyA = sum(~isnan(fA(:)) & isnan(fB(:)));
            onlyB = sum(isnan(fA(:)) & ~isnan(fB(:)));
            fprintf(' %s kHz: median %+.2f, 5-95%% range %+.2f to %+.2f; bins only in A: %d, only in B: %d\n', ...
                freqs{i}, median(d(both)), prctile(d(both), 5), prctile(d(both), 95), onlyA, onlyB);
        else
            axis(ax, 'off');
            text(ax, 0.5, 0.5, sprintf('Different number of dives (%d vs %d)', size(fA, 1), size(fB, 1)), ...
                'HorizontalAlignment', 'center');
        end
    end

    outFile = fullfile(outFolder, sprintf('compare-%s-%s.png', kinds{k}, dateOfData));
    drawnow;
    print(fig, '-dpng', outFile, '-r0');
    fprintf(' saved %s\n', outFile);
end

%% local functions
function drawPanel(ax, M, limits, cmap, maxDepth)
imagesc(ax, M', 'AlphaData', ~isnan(M'));
set(ax, 'Color', [0.85, 0.85, 0.85]); % grey = no data
colormap(ax, cmap);
clim(ax, limits);
ylim(ax, [0, maxDepth]);
xlabel(ax, 'Dive Number');
ylabel(ax, 'Depth (m)');
end

function cmap = blueWhiteRed(n)
h = floor(n / 2);
up = linspace(0, 1, h)';
down = linspace(1, 0, n - h)';
cmap = [up, up, ones(h, 1); ones(n - h, 1), down, down];
end
