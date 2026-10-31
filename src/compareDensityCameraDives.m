%% Camera effect on numerical density (ind./m^3): camera vs non-camera dives
% Compares density profiles from dives with the shadowgraph camera on and
% off, at 455 and 769 kHz (no camera noise at these frequencies, so any
% difference is the night-time deficit near the glider), for day and night,
% for all days together and for each day separately.
%
% Run AZFP_Initial_Processing and Merge_Dives_Data first (uses
% output/0-merged-dives.mat). Camera on/off per dive comes from the glider
% instrument schedule in CTD_processor.
%
% Density per dive and depth = Dives(i).P(j).mN (masked sv / sigma_bs,
% per-depth median or mean over pings, see calcNumericalDensity.m).
% Dives are grouped by local date (UTC-3); a night (23:00-04:59) is labelled
% with the date on which it starts (night of 24 Jul = 24 Jul 23:00 to
% 25 Jul 04:59). Twilight dives ("Other") are not used.
%
% Two extra checks separate a real acoustic deficit from the copepod mask:
%   unmasked density = median sv over pings / sigma_bs (no dB-difference
%   mask), and the fraction of samples that pass the mask. The masks are
%   455-200 kHz (for 455 kHz) and 769-455 kHz (for 769 kHz), so camera noise
%   at 200 kHz and a deficit that differs between frequencies both change
%   which samples pass.
%
% Outputs (in output/):
%   figure-density-camera-effect-profiles.png  all days: profiles + difference
%   figure-density-camera-effect-by-day.png    per day: mean density, camera vs not
%   figure-density-camera-effect-daily-profiles.png  per day profiles
%   figure-density-camera-effect-ratios.png    camera / no camera per day:
%                                              masked vs unmasked density
%   figure-density-camera-effect-mask.png      all days: unmasked profiles and
%                                              fraction of samples passing the mask
%   density_camera_effect_by_day.csv           summary table (also printed)

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
summaryDepth        = [5, 100];   % m, depth range for the mean density per dive group
                                  % (top 5 m left out: surface bubbles/noise)
minDives            = 3;          % fewer dives in a group -> no profile/summary

P_indices  = [3, 4];              % Dives(i).P(j): 3 = 455 kHz, 4 = 769 kHz
freqLabels = {'455 kHz', '769 kHz'};
sigma_bs   = [7.6e-11, 7.9e-11];  % as in calcNumericalDensity.m (455, 769 kHz)
colCam     = [0.85, 0.2, 0.2];    % camera on
colNoCam   = [0.1, 0.35, 0.8];    % camera off

%% Load merged dives and camera state
S = load(fullfile(config.outputFolder, '0-merged-dives.mat'), 'Dives');
Dives = S.Dives;
nDives = numel(Dives);
cam = getCameraDives(Dives, config);
if isempty(cam)
    error('Camera schedule not found: %s', config.scheduleFile);
end
cam = cam(:)';

tod = string({Dives.tod});                          % "Day" / "Night" / "Other"
tLocal = [Dives.starttime] - 3/24;                  % datenum, local time (UTC-3)
groupDate = floor(tLocal - 6/24);                   % night after midnight -> previous date
depths = (1:config.maxDepth)';

%% Matrices (depth x dive): masked density, unmasked density, mask pass fraction
N = cell(1, 2); U = cell(1, 2); F = cell(1, 2);
for f = 1:2
    j = P_indices(f);
    [M, Mu, Mf] = deal(NaN(config.maxDepth, nDives));
    for k = 1:nDives
        P = Dives(k).P(j);
        p = P.mN(:);                              % masked density profile
        n = min(numel(p), config.maxDepth);
        M(1:n, k) = p(1:n);
        u = P.msv(:) / sigma_bs(f);               % unmasked: median sv / sigma_bs
        n = min(numel(u), config.maxDepth);
        Mu(1:n, k) = u(1:n);
        valid = sum(~isnan(P.sv), 1);             % samples with a value (pings x depth)
        pass  = sum(~isnan(P.masked) & P.masked > 0, 1);
        fr = pass(:) ./ valid(:);
        fr(valid(:) == 0) = NaN;
        n = min(numel(fr), config.maxDepth);
        Mf(1:n, k) = fr(1:n);
    end
    N{f} = M; U{f} = Mu; F{f} = Mf;
end

periods = ["Day", "Night"];
inSummary = depths >= summaryDepth(1) & depths <= summaryDepth(2);

% mean over depth of the median profile, using only depths where both
% camera and non-camera profiles exist (so the two are compared like for like)
pairedMean = @(a, b) deal(mean(a(~isnan(a) & ~isnan(b) & inSummary)), ...
                          mean(b(~isnan(a) & ~isnan(b) & inSummary)));

%% Figure 1: all days, day and night profiles + difference
fig1 = figure(1);
try, theme(fig1, 'light'); end
set(fig1, 'Color', 'w', 'InvertHardcopy', 'off', 'Position', [80, 80, 1300, 850])
tl1 = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl1, 'Camera effect on numerical density, 24-29 July 2025 (median and 25-75 % range across dives)')

fprintf('\nAll days together (mean density %g-%g m, ind./m^3):\n', summaryDepth);
fprintf(' %-8s %-6s  n cam  n off   camera   no camera   difference   camera / no camera\n', 'freq', 'period');
allRows = {};
for f = 1:2
    dLine = cell(1, 2);
    for p = 1:2
        nexttile((f - 1) * 3 + p)
        hold on
        selC = tod == periods(p) & cam;
        selN = tod == periods(p) & ~cam;
        [mC, qC] = drawProfile(N{f}(:, selC), depths, colCam,   sprintf('camera (n = %d)', sum(selC)));
        [mN, qN] = drawProfile(N{f}(:, selN), depths, colNoCam, sprintf('no camera (n = %d)', sum(selN)));
        hold off
        setProfileXLim([qC, qN], depths, summaryDepth(1))
        set(gca, 'YDir', 'reverse'); ylim([0, config.maxDepth]); box on
        xlabel('Ind. m^{-3}'); ylabel('Depth (m)')
        title(sprintf('%s, %s', freqLabels{f}, lower(periods(p))))
        legend('Location', 'southeast')
        dLine{p} = mC - mN;
        [a, b] = pairedMean(mC, mN);
        fprintf(' %-8s %-6s  %5d  %5d   %7.0f   %9.0f   %+10.0f   %8.2f\n', freqLabels{f}, ...
            periods(p), sum(selC), sum(selN), a, b, a - b, a / b);
        [au, bu] = pairedMean(median(U{f}(:, selC), 2, 'omitnan'), median(U{f}(:, selN), 2, 'omitnan'));
        [af, bf] = pairedMean(median(F{f}(:, selC), 2, 'omitnan'), median(F{f}(:, selN), 2, 'omitnan'));
        allRows(end+1, :) = {freqLabels{f}, char(periods(p)), a / b, au / bu, af, bf}; %#ok<SAGROW>
    end
    nexttile((f - 1) * 3 + 3)
    hold on
    xline(0, 'k-', 'HandleVisibility', 'off')
    plot(dLine{1}, depths, '-', 'Color', [1, 0.6, 0], 'LineWidth', 1.5, 'DisplayName', 'day')
    plot(dLine{2}, depths, '-', 'Color', [0, 0, 0.5], 'LineWidth', 1.5, 'DisplayName', 'night')
    hold off
    dv = [dLine{1}, dLine{2}]; dv = dv(depths >= summaryDepth(1), :);
    dm = max(abs(dv(:)), [], 'omitnan');
    if isfinite(dm) && dm > 0, xlim(1.05 * dm * [-1, 1]); end
    set(gca, 'YDir', 'reverse'); ylim([0, config.maxDepth]); box on
    xlabel('Camera - no camera (ind. m^{-3})'); ylabel('Depth (m)')
    title(sprintf('%s, difference of medians', freqLabels{f}))
    legend('Location', 'southeast')
end
drawnow
print(fig1, '-dpng', fullfile(config.outputFolder, 'figure-density-camera-effect-profiles.png'), '-r300')

fprintf('\nMasked vs unmasked, all days (camera / no camera; mask pass fraction %g-%g m):\n', summaryDepth);
fprintf(' %-8s %-6s  ratio masked  ratio unmasked   pass camera  pass no camera\n', 'freq', 'period');
for r = 1:size(allRows, 1)
    fprintf(' %-8s %-6s  %12.2f  %14.2f   %11.2f  %14.2f\n', allRows{r, :});
end

%% Per day: summary table
dates = unique(groupDate(tod ~= "Other"));
rows = {};
fprintf('\nPer day (mean density %g-%g m, ind./m^3; night = night starting on that date):\n', summaryDepth);
fprintf(' %-7s %-8s %-6s  n cam  n off   camera   no camera   difference   camera / no camera\n', 'date', 'freq', 'period');
for d = 1:numel(dates)
    for p = 1:2
        for f = 1:2
            selC = groupDate == dates(d) & tod == periods(p) & cam;
            selN = groupDate == dates(d) & tod == periods(p) & ~cam;
            [a, b, au, bu, af, bf] = deal(NaN);
            if sum(selC) >= minDives && sum(selN) >= minDives
                mC = median(N{f}(:, selC), 2, 'omitnan');
                mN = median(N{f}(:, selN), 2, 'omitnan');
                [a, b] = pairedMean(mC, mN);
                [au, bu] = pairedMean(median(U{f}(:, selC), 2, 'omitnan'), median(U{f}(:, selN), 2, 'omitnan'));
                [af, bf] = pairedMean(median(F{f}(:, selC), 2, 'omitnan'), median(F{f}(:, selN), 2, 'omitnan'));
            end
            fprintf(' %-7s %-8s %-6s  %5d  %5d   %7.0f   %9.0f   %+10.0f   %8.2f\n', ...
                datestr(dates(d), 'dd mmm'), freqLabels{f}, periods(p), sum(selC), sum(selN), a, b, a - b, a / b);
            rows(end+1, :) = {datetime(dates(d), 'ConvertFrom', 'datenum'), freqLabels{f}, ...
                char(periods(p)), sum(selC), sum(selN), a, b, a - b, a / b, ...
                au, bu, au / bu, af, bf}; %#ok<SAGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'date_local', 'frequency', 'period', 'n_camera', ...
    'n_no_camera', 'N_camera', 'N_no_camera', 'N_difference', 'ratio_camera_no_camera', ...
    'Nunmasked_camera', 'Nunmasked_no_camera', 'ratio_unmasked', 'pass_frac_camera', 'pass_frac_no_camera'});
T = T(~isnan(T.N_camera), :);     % drop dates/periods without enough dives
dates = unique(datenum(T.date_local))';
writetable(T, fullfile(config.outputFolder, 'density_camera_effect_by_day.csv'));

%% Figure 2: per day, mean density camera vs no camera
fig2 = figure(2);
try, theme(fig2, 'light'); end
set(fig2, 'Color', 'w', 'InvertHardcopy', 'off', 'Position', [80, 80, 1200, 750])
tl2 = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl2, sprintf('Mean density %g-%g m per day: camera vs non-camera dives', summaryDepth))
x = 1:numel(dates);
for f = 1:2
    for p = 1:2
        nexttile((f - 1) * 2 + p)
        sel = strcmp(T.frequency, freqLabels{f}) & strcmp(T.period, char(periods(p)));
        xs = x(1:sum(sel));
        bar(xs, [T.N_no_camera(sel), T.N_camera(sel)], 'grouped')
        colororder(gca, [colNoCam; colCam])
        ylabel('Ind. m^{-3}')
        dd = datenum(T.date_local(sel));
        xticks(xs); xticklabels(cellstr(datestr(dd, 'dd mmm')))
        title(sprintf('%s, %s', freqLabels{f}, lower(periods(p))))
        r = T.ratio_camera_no_camera(sel);
        top = max([T.N_no_camera(sel), T.N_camera(sel)], [], 2);
        ylim([0, 1.15 * max(top)])
        for i = 1:numel(xs)
            if ~isnan(r(i))
                text(xs(i), top(i), sprintf('%.2f', r(i)), 'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom', 'FontSize', 8)
            end
        end
        if f == 1 && p == 1
            legend({'no camera', 'camera'}, 'Location', 'northoutside', 'Orientation', 'horizontal')
        end
        box on
    end
end
xlabel(tl2, 'Date (night = night starting on that date); numbers = camera / no camera')
drawnow
print(fig2, '-dpng', fullfile(config.outputFolder, 'figure-density-camera-effect-by-day.png'), '-r300')

%% Figure 3: per day profiles (rows: 455 day, 455 night, 769 day, 769 night)
fig3 = figure(3);
try, theme(fig3, 'light'); end
set(fig3, 'Color', 'w', 'InvertHardcopy', 'off', 'Position', [40, 40, 1700, 1100])
tl3 = tiledlayout(4, numel(dates), 'TileSpacing', 'tight', 'Padding', 'compact');
title(tl3, 'Daily median density profiles: camera (red) vs non-camera (blue) dives')
for f = 1:2
    for p = 1:2
        for d = 1:numel(dates)
            nexttile
            selC = groupDate == dates(d) & tod == periods(p) & cam;
            selN = groupDate == dates(d) & tod == periods(p) & ~cam;
            hold on
            pv = NaN(config.maxDepth, 2);
            if sum(selC) >= minDives, pv(:, 1) = median(N{f}(:, selC), 2, 'omitnan'); plot(pv(:, 1), depths, '-', 'Color', colCam, 'LineWidth', 1.2); end
            if sum(selN) >= minDives, pv(:, 2) = median(N{f}(:, selN), 2, 'omitnan'); plot(pv(:, 2), depths, '-', 'Color', colNoCam, 'LineWidth', 1.2); end
            hold off
            setProfileXLim(pv, depths, summaryDepth(1))
            set(gca, 'YDir', 'reverse', 'FontSize', 7); ylim([0, config.maxDepth]); box on
            title(sprintf('%s %s, %s (%d/%d)', freqLabels{f}, lower(periods(p)), ...
                datestr(dates(d), 'dd mmm'), sum(selC), sum(selN)), 'FontSize', 7)
            if d == 1, ylabel('Depth (m)'); end
        end
    end
end
xlabel(tl3, 'Ind. m^{-3}   (numbers in titles = camera / non-camera dives)')
drawnow
print(fig3, '-dpng', fullfile(config.outputFolder, 'figure-density-camera-effect-daily-profiles.png'), '-r300')

%% Figure 4: camera / no camera ratio per day, masked vs unmasked density
fig4 = figure(4);
try, theme(fig4, 'light'); end
set(fig4, 'Color', 'w', 'InvertHardcopy', 'off', 'Position', [80, 80, 1200, 480])
tl4 = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl4, sprintf('Camera / no-camera mean density %g-%g m: with and without the copepod mask', summaryDepth))
pcol = struct('Day', [1, 0.6, 0], 'Night', [0, 0, 0.5]);
for f = 1:2
    nexttile
    hold on
    yline(1, 'k-', 'HandleVisibility', 'off')
    for p = 1:2
        sel = strcmp(T.frequency, freqLabels{f}) & strcmp(T.period, char(periods(p)));
        dd = datenum(T.date_local(sel));
        c = pcol.(char(periods(p)));
        plot(dd, T.ratio_camera_no_camera(sel), '-o', 'Color', c, 'MarkerFaceColor', c, ...
            'LineWidth', 1.5, 'DisplayName', [lower(char(periods(p))), ', masked (copepods)'])
        plot(dd, T.ratio_unmasked(sel), '--s', 'Color', c, 'LineWidth', 1.2, ...
            'DisplayName', [lower(char(periods(p))), ', unmasked (all scatterers)'])
    end
    hold off
    set(gca, 'YScale', 'log'); ylim([0.05, 60]); yticks([0.05 0.1 0.2 0.5 1 2 5 10 20 50]); box on; grid on
    xticks(dates); xticklabels(cellstr(datestr(dates, 'dd mmm')))
    ylabel('Camera / no camera')
    title(freqLabels{f})
    if f == 2, legend('Location', 'southoutside', 'NumColumns', 2); end
end
xlabel(tl4, 'Date (night = night starting on that date); 1 = no camera effect')
drawnow
print(fig4, '-dpng', fullfile(config.outputFolder, 'figure-density-camera-effect-ratios.png'), '-r300')

%% Figure 5: all days, unmasked density profiles and mask pass fraction
fig5 = figure(5);
try, theme(fig5, 'light'); end
set(fig5, 'Color', 'w', 'InvertHardcopy', 'off', 'Position', [80, 80, 1300, 850])
tl5 = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl5, 'All days: unmasked density and fraction of samples passing the copepod mask (median, 25-75 %)')
for f = 1:2
    for p = 1:2
        selC = tod == periods(p) & cam;
        selN = tod == periods(p) & ~cam;
        nexttile((f - 1) * 4 + p)
        hold on
        [~, qC] = drawProfile(U{f}(:, selC), depths, colCam,   'camera');
        [~, qN] = drawProfile(U{f}(:, selN), depths, colNoCam, 'no camera');
        hold off
        setProfileXLim([qC, qN], depths, summaryDepth(1))
        set(gca, 'YDir', 'reverse'); ylim([0, config.maxDepth]); box on
        xlabel('Ind. m^{-3} (unmasked)'); ylabel('Depth (m)')
        title(sprintf('%s %s, unmasked', freqLabels{f}, lower(periods(p))))
        legend('Location', 'southeast')
        nexttile((f - 1) * 4 + 2 + p)
        hold on
        drawProfile(F{f}(:, selC), depths, colCam,   'camera');
        drawProfile(F{f}(:, selN), depths, colNoCam, 'no camera');
        hold off
        set(gca, 'YDir', 'reverse'); ylim([0, config.maxDepth]); xlim([0, 1]); box on
        xlabel('Fraction passing mask'); ylabel('Depth (m)')
        title(sprintf('%s %s, mask', freqLabels{f}, lower(periods(p))))
        legend('Location', 'southeast')
    end
end
drawnow
print(fig5, '-dpng', fullfile(config.outputFolder, 'figure-density-camera-effect-mask.png'), '-r300')

fprintf('\nFigures and density_camera_effect_by_day.csv saved in %s\n', config.outputFolder);

%% Local function
function [md, q3] = drawProfile(M, depths, c, name)
% median profile (line) and 25-75 % range (shading) across dives
md = median(M, 2, 'omitnan');
q3 = md;
if size(M, 2) < 2, plot(md, depths, '-', 'Color', c, 'DisplayName', name); return; end
q1 = prctile(M, 25, 2);
q3 = prctile(M, 75, 2);
ok = ~isnan(md) & ~isnan(q1) & ~isnan(q3);
fill([q1(ok); flipud(q3(ok))], [depths(ok); flipud(depths(ok))], c, ...
    'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off')
plot(md, depths, '-', 'Color', c, 'LineWidth', 1.5, 'DisplayName', name)
end

function setProfileXLim(vals, depths, minDepth)
% x-axis from 0 to the largest value below minDepth (surface spikes ignored)
v = vals(depths >= minDepth, :);
xm = max(v(:), [], 'omitnan');
if ~isempty(xm) && isfinite(xm) && xm > 0, xlim([0, 1.05 * xm]); end
end
