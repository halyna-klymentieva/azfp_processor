%% Numerical density (ind./m^3) for dives WITHOUT the shadowgraph camera
% Uses the merged dives file output/0-merged-dives.mat (run
% AZFP_Initial_Processing and Merge_Dives_Data first) and the glider
% instrument schedule from CTD_processor (camera on/off per descent).
%
% Why: on camera dives the AZFP data are affected (night: fewer scatterers
% near the glider at all frequencies; day: camera noise at 130/200 kHz),
% so densities are taken from non-camera dives only.
%
% Density per dive and depth = Dives(i).P(j).mN (masked sv / sigma_bs,
% per-depth median or mean over pings, see calcNumericalDensity.m).
% Frequencies: 200, 455, 769 kHz (130 kHz has no density).
%
% Outputs (in output/):
%   figure-density-noncamera-sections.png   density vs dive and depth
%   figure-density-noncamera-profiles.png   day vs night median profiles
%   density_noncamera_dives.csv             all values, one row per dive and depth

cd(fileparts(mfilename('fullpath')));
addpath(genpath(pwd));
rehash path
clc
clear variables
close all

%% Settings
config.outputFolder  = fullfile(pwd, '..', 'output');
config.scheduleFile  = fullfile(pwd, '..', '..', 'CTD_processor', 'CTD_processed', 'instrument_schedule_by_cast.csv');
config.maxDepth      = 105;
config.densityLimits = [0, 10000];   % colour scale for the sections (ind./m^3)

P_indices  = [2, 3, 4];
freqLabels = {'200 kHz', '455 kHz', '769 kHz'};

%% Load merged dives
S = load(fullfile(config.outputFolder, '0-merged-dives.mat'), 'Dives');
Dives = S.Dives;
nDives = numel(Dives);

%% Camera on/off per dive -> keep non-camera dives
cam = getCameraDives(Dives, config);
if isempty(cam)
    error('Camera schedule not found: %s', config.scheduleFile);
end
keep = find(~cam);
fprintf('%d dives: %d with camera on, %d without camera (used here).\n', ...
    nDives, sum(cam), numel(keep));

tStart = [Dives(keep).starttime];                 % datenum, UTC
tod    = string({Dives(keep).tod});               % "Day" / "Night" / "Other"
depths = (1:config.maxDepth)';
nk     = numel(keep);

fprintf('Non-camera dives: %d day, %d night, %d twilight.\n', ...
    sum(tod == "Day"), sum(tod == "Night"), sum(tod == "Other"));

%% Density matrices (depth x dive), one per frequency
N = cell(1, 3);
for f = 1:3
    j = P_indices(f);
    M = NaN(config.maxDepth, nk);
    for k = 1:nk
        p = Dives(keep(k)).P(j).mN(:);
        n = min(numel(p), config.maxDepth);
        M(1:n, k) = p(1:n);
    end
    N{f} = M;
end

%% Figure 1: density sections (non-camera dives only)
fig1 = figure(1);
try, theme(fig1, 'light'); end   % white background (MATLAB R2025a+ dark theme)
set(fig1, 'Color', 'w', 'InvertHardcopy', 'off')
set(fig1, 'Position', [100, 100, 1200, 900])

% x-axis ticks at the first dive of each day
dayStart = [1, find(diff(floor(tStart)) > 0) + 1];
isNight  = find(tod == "Night");
isDay    = find(tod == "Day");

for f = 1:3
    ax = subplot(3, 1, f);
    imagesc(1:nk, depths, N{f})
    set(ax, 'YDir', 'reverse')
    caxis(config.densityLimits)
    colormap(ax, hot)
    cb = colorbar;
    cb.Label.String = 'Ind. m^{-3}';
    hold on
    % day/night marks below the data: orange = day, dark blue = night
    plot(isDay,   (config.maxDepth + 3) * ones(size(isDay)),   's', 'MarkerSize', 3, ...
        'MarkerFaceColor', [1, 0.6, 0], 'MarkerEdgeColor', 'none')
    plot(isNight, (config.maxDepth + 3) * ones(size(isNight)), 's', 'MarkerSize', 3, ...
        'MarkerFaceColor', [0, 0, 0.5], 'MarkerEdgeColor', 'none')
    hold off
    ylim([0, config.maxDepth + 6])
    yticks(10:10:config.maxDepth)
    ylabel('Depth (m)')
    xticks(dayStart)
    xticklabels(cellstr(datestr(tStart(dayStart), 'dd mmm')))
    title([freqLabels{f}, ' - non-camera dives only'])
end
xlabel('Non-camera dives (orange = day, dark blue = night)')

drawnow
print(fig1, '-dpng', fullfile(config.outputFolder, 'figure-density-noncamera-sections.png'), '-r300')

%% Figure 2: day vs night profiles (median and 25-75 % range across dives)
fig2 = figure(2);
try, theme(fig2, 'light'); end   % white background (MATLAB R2025a+ dark theme)
set(fig2, 'Color', 'w', 'InvertHardcopy', 'off')
set(fig2, 'Position', [100, 100, 1200, 500])
cols = struct('Day', [1, 0.6, 0], 'Night', [0, 0, 0.5]);
periods = ["Day", "Night"];

for f = 1:3
    subplot(1, 3, f)
    hold on
    h = gobjects(1, 2);
    for p = 1:2
        sel = tod == periods(p);
        if ~any(sel), continue; end
        M  = N{f}(:, sel);
        md = median(M, 2, 'omitnan');
        q1 = prctile(M, 25, 2);
        q3 = prctile(M, 75, 2);
        ok = ~isnan(md) & ~isnan(q1) & ~isnan(q3);
        c  = cols.(char(periods(p)));
        fill([q1(ok); flipud(q3(ok))], [depths(ok); flipud(depths(ok))], c, ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off')
        h(p) = plot(md, depths, '-', 'Color', c, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('%s (n = %d)', periods(p), sum(sel)));
    end
    hold off
    set(gca, 'YDir', 'reverse')
    ylim([0, config.maxDepth])
    xlabel('Ind. m^{-3}')
    ylabel('Depth (m)')
    title(freqLabels{f})
    legend('Location', 'southeast')
    box on
end
sgtitle('Non-camera dives: median (line) and 25-75 % range (shading)')

drawnow
print(fig2, '-dpng', fullfile(config.outputFolder, 'figure-density-noncamera-profiles.png'), '-r300')

%% Export: one row per non-camera dive and 1 m depth bin
T = table( ...
    repelem(keep(:), config.maxDepth), ...
    repelem(datetime(tStart(:), 'ConvertFrom', 'datenum', 'TimeZone', 'UTC'), config.maxDepth), ...
    repelem(tod(:), config.maxDepth), ...
    repmat(depths, nk, 1), ...
    N{1}(:), N{2}(:), N{3}(:), ...
    'VariableNames', {'dive', 'start_time_utc', 'tod', 'depth_m', ...
                      'N_200kHz', 'N_455kHz', 'N_769kHz'});
writetable(T, fullfile(config.outputFolder, 'density_noncamera_dives.csv'));
fprintf('Saved %d rows to density_noncamera_dives.csv\n', height(T));
