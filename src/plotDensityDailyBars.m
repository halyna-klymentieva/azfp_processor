%% Daily day/night bar chart of numerical density (one frequency), non-camera dives
% For each dive: mean density over depthRange (from Dives(i).P(j).mN).
% For each date and period (day / night), barStatistic sets the bars:
%   'mean'   : bar = mean over dives, error bar = standard deviation
%   'median' : bar = median over dives, error bar = 25th-75th percentile
% Circles = outlier dives (below Q1 - 1.5 IQR or above Q3 + 1.5 IQR within
% that day or night, Tukey 1977). Outliers are only shown, not removed.
%
% Uses output/0-merged-dives.mat (run AZFP_Initial_Processing and
% Merge_Dives_Data first) and the camera schedule from CTD_processor;
% camera dives are left out.
% Day = 07:00-20:59, night = 23:00-04:59 local time (Dives(i).tod, set in
% plotDayNight.m); twilight dives are not used. A night is labelled with the
% date it starts (night of 24 Jul = 24 Jul 23:00 to 25 Jul 04:59).
%
% Outputs: output/figure-density-<freq>-daily-bars.png and
%          output/density_<freq>_daily_bars.csv (also printed)

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

P_index    = 3;           % Dives(i).P(j): 2 = 200 kHz, 3 = 455 kHz, 4 = 769 kHz
freqLabel  = '455 kHz';
depthRange = [5, 80];     % m, depths averaged for each dive (top 5 m left out)
minDives   = 3;           % dives needed for a bar
plotType   = 'box';       % 'box' (box plots) or 'bar' (bars + error bars)
barStatistic = 'median';  % for 'bar': 'median' (error bars 25-75 %) or 'mean' (error bars SD)
showDives  = true;        % for 'box': also draw every dive as a small grey dot
diveValue  = 'integrated';  % one value per dive: 'mean' (ind./m^3 over depthRange) or
                           % 'integrated' (sum over the 1 m bins = ind./m^2; bins without
                           % data or with no copepod-like samples count as 0)
useLog10   = true;         % true = plot log10(value + 1), the usual way to show abundance;
                           % box, quartiles, outliers and means are then on the log scale
showMean   = true;        % for 'box': mark the mean of each day / night with a white diamond

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

% one value per dive over depthRange: mean (ind./m^3) or integrated (ind./m^2)
vDive = NaN(1, numel(keep));
for k = 1:numel(keep)
    p = Dives(keep(k)).P(P_index).mN(:);
    z = (1:numel(p))';
    pz = p(z >= depthRange(1) & z <= depthRange(2));
    if strcmpi(diveValue, 'integrated')
        if any(~isnan(pz)), vDive(k) = sum(pz, 'omitnan'); end   % 1 m bins -> ind./m^2
    else
        vDive(k) = mean(pz, 'omitnan');
    end
end

if useLog10
    vDive = log10(vDive + 1);
end

%% Statistics per date and period
periods = ["Day", "Night"];
cols    = {[0.93, 0.55, 0.00], [0.10, 0.20, 0.55]};
offsets = [-0.2, 0.2];
use     = tod == "Day" | tod == "Night";
dates   = unique(groupDate(use));
nd      = numel(dates);
[mu, sd, nn, eLo, eHi] = deal(NaN(nd, 2));
outl = cell(nd, 2);
rows = {};
for p = 1:2
    for d = 1:nd
        v = vDive(tod == periods(p) & groupDate == dates(d));
        v = v(~isnan(v));
        if numel(v) < minDives, continue; end
        q   = prctile(v, [25, 75]);
        iqr_ = q(2) - q(1);
        o   = v(v < q(1) - 1.5 * iqr_ | v > q(2) + 1.5 * iqr_);
        if strcmpi(barStatistic, 'median')
            mu(d, p) = median(v); eLo(d, p) = mu(d, p) - q(1); eHi(d, p) = q(2) - mu(d, p);
        else
            mu(d, p) = mean(v); eLo(d, p) = std(v); eHi(d, p) = std(v);
        end
        sd(d, p) = std(v); nn(d, p) = numel(v); outl{d, p} = o;
        rows(end+1, :) = {datetime(dates(d), 'ConvertFrom', 'datenum'), char(periods(p)), ...
            numel(v), mean(v), std(v), median(v), q(1), q(2), numel(o), mat2str(round(o))}; %#ok<SAGROW>
    end
end

%% Figure
fig = figure(1);
try, theme(fig, 'light'); end                     % white background (R2025a+ dark theme)
set(fig, 'Color', 'w', 'InvertHardcopy', 'off', 'Position', [100, 100, 900, 550])
hold on
hb = gobjects(1, 2);
ho = gobjects(0);
if strcmpi(plotType, 'box')
    % box = 25-75 %, line = median, whiskers = furthest dive within
    % 1.5 IQR, circles = outliers (beyond 1.5 IQR)
    for p = 1:2
        sel = find(tod == periods(p) & ~isnan(vDive));
        xi  = arrayfun(@(g) find(dates == g), groupDate(sel));
        use = ~isnan(nn(xi, p))';                 % groups with >= minDives dives
        xv  = xi(use) + offsets(p);
        yv  = vDive(sel(use));
        hb(p) = boxchart(xv, yv, 'BoxWidth', 0.35, 'BoxFaceColor', cols{p}, ...
            'BoxFaceAlpha', 0.45, 'WhiskerLineColor', [0.2, 0.2, 0.2], ...
            'MarkerStyle', 'o', 'MarkerColor', 'k', 'DisplayName', lower(char(periods(p))));
        if showDives
            jit = (rand(size(xv)) - 0.5) * 0.18;
            scatter(xv + jit, yv, 9, [0.45, 0.45, 0.45], 'filled', ...
                'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off')
        end
        if showMean
            ok = ~isnan(mu(:, p))';
            xm = (1:nd) + offsets(p);
            if strcmpi(barStatistic, 'median')   % mu holds medians -> recompute means
                mm = arrayfun(@(d) mean(vDive(tod == periods(p) & groupDate == dates(d)), 'omitnan'), 1:nd);
            else
                mm = mu(:, p)';
            end
            hm = plot(xm(ok), mm(ok), 'd', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', ...
                'MarkerSize', 7, 'LineStyle', 'none', 'DisplayName', 'mean');
            if p == 1, hm.HandleVisibility = 'off'; end
            ho = hm;
        end
    end
else
for p = 1:2
    x = (1:nd) + offsets(p);
    ok = ~isnan(mu(:, p))';
    hb(p) = bar(x(ok), mu(ok, p)', 0.38, 'FaceColor', cols{p}, 'EdgeColor', 'none', ...
        'DisplayName', lower(char(periods(p))));
    errorbar(x(ok), mu(ok, p)', eLo(ok, p)', eHi(ok, p)', 'k', 'LineStyle', 'none', 'CapSize', 8, ...
        'LineWidth', 1, 'HandleVisibility', 'off');
    for d = find(ok)
        o = outl{d, p};
        if ~isempty(o)
            ho = plot(x(d) * ones(size(o)), o, 'o', 'MarkerFaceColor', 'w', ...
                'MarkerEdgeColor', 'k', 'MarkerSize', 6, 'DisplayName', 'outlier (1.5 \times IQR)');
        end
    end
end
end
hold off
xticks(1:nd)
xticklabels(cellstr(datestr(dates, 'dd mmm')))
xlim([0.4, nd + 0.6])
% y-axis: from 0 to 8 % above the highest dive or error bar, so outliers
% are not cut off at the top edge
if strcmpi(plotType, 'box')
    yTop = max(vDive(ismember(tod, periods)), [], 'omitnan');
else
    yTop = max([mu(:) + eHi(:); cell2mat(outl(:)')'], [], 'omitnan');
end
if useLog10
    yBot = min(vDive(ismember(tod, periods)), [], 'omitnan');
    if ~isempty(yTop) && isfinite(yTop) && isfinite(yBot)
        pad = 0.08 * (yTop - yBot);
        ylim([yBot - pad, yTop + pad]);
    end
elseif ~isempty(yTop) && isfinite(yTop) && yTop > 0
    ylim([0, 1.08 * yTop]);
end
grid on; set(gca, 'XGrid', 'off', 'Box', 'on')
if strcmpi(diveValue, 'integrated')
    yl = sprintf('integrated density %d-%d m (ind. m^{-2}', depthRange);
else
    yl = sprintf('dive mean density %d-%d m (ind. m^{-3}', depthRange);
end
if useLog10
    ylabel(['log_{10}(', yl, ' + 1))'])
else
    ylabel([upper(yl(1)), yl(2:end), ')'])
end
xlabel('Date (night = night starting on that date)')
if strcmpi(plotType, 'box')
    title(sprintf('%s, non-camera dives: one box per day / night (n = dives)', freqLabel))
    subtitle('box = 25-75 %, line = median, diamond = mean, whiskers = 1.5 IQR, circles = outliers, dots = dives', 'FontSize', 8)
elseif strcmpi(barStatistic, 'median')
    title(sprintf('%s, non-camera dives: daily median, error bars = 25-75 %%', freqLabel))
else
    title(sprintf('%s, non-camera dives: daily mean \\pm SD', freqLabel))
end
if isempty(ho), legend(hb, 'Location', 'northeast'); else, legend([hb, ho], 'Location', 'northeast'); end

drawnow
tag = strrep(freqLabel, ' kHz', '');
if strcmpi(plotType, 'box'), figTag = 'box'; else, figTag = ['bars-', lower(barStatistic)]; end
if strcmpi(diveValue, 'integrated'), figTag = [figTag, '-integrated']; end
if useLog10, figTag = [figTag, '-log10']; end
print(fig, '-dpng', fullfile(config.outputFolder, sprintf('figure-density-%s-daily-%s.png', tag, figTag)), '-r300')

%% Table
T = cell2table(rows, 'VariableNames', {'date_local', 'period', 'n_dives', 'mean', 'sd', ...
    'median', 'q25', 'q75', 'n_outliers', 'outlier_values'});
T = sortrows(T, {'period', 'date_local'});
disp(T)
csvName = sprintf('density_%s_daily_bars_%s.csv', tag, lower(diveValue));
if useLog10, csvName = strrep(csvName, '.csv', '_log10.csv'); end
writetable(T, fullfile(config.outputFolder, csvName));
fprintf('Saved figure-density-%s-daily-%s.png and %s\n', tag, figTag, csvName);
