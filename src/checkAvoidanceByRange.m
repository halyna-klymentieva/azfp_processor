%% Avoidance check: does Sv depend on distance from the glider? (non-camera dives)
% The processed densities use only 2-10 m below the glider (200/455 kHz,
% 2-5 m at 769 kHz; see config.farFieldCutOffRange). If animals move away
% from the glider in daylight, daytime Sv close to the glider should be
% lower than further away IN THE SAME WATER DEPTH, while at night it should
% not change with range.
%
% For each ping of a NON-camera descent, Sv is averaged (linear) in range
% bins from the transducer, separately for fixed water-depth bands, so every
% range bin looks at the same water. Medians over pings are compared for day
% (07:00-20:59) and night (23:00-04:59 local, UTC-3; as plotDayNight.m).
% Each band has its own bottom margin: samples closer than that to the
% seafloor (glider m_water_depth) are left out. A band bottom of Inf means
% "down to the seafloor minus the margin". Samples shallower than surfaceMin
% are always left out.
%
% Noise subtraction (De Robertis & Higginbottom 2007 style): the raw Sv is
% signal + noise, and noise grows with range like the range correction
% (20 log10 R + 2 alpha R). For each frequency and day the noise level N0 is
% estimated from the non-camera pings: Sv minus the range correction is
% averaged in blocks of noiseBlockPings pings x each range bin (water column
% only), the minimum over range is taken for each block, and N0 = median of
% these minima. The noise Sv(R) = N0 + 20 log10 R + 2 alpha R is subtracted
% in linear units from every sample (negative results -> 0). Corrected
% medians that are 0 (noise-dominated) are shown as NaN. The share of samples
% at least snrThreshold dB above the noise is also reported.
%
% Raw files are re-read with ProcessAZFP (no filters, no calibration offset,
% fixed sound speed - these cancel in day/night and range comparisons).
% Takes a few minutes per day. Run with Run / F5.
% Outputs (output/compare/): figure-avoidance-by-range.png and
% avoidance_by_range.csv (also printed).

cd(fileparts(mfilename('fullpath')));
addpath(genpath(pwd));
rehash path
clear variables
close all

%% Settings
dates = {'25-07-24', '25-07-25', '25-07-26', '25-07-27', '25-07-28', '25-07-29'};
xmlFileName = '25070417.XML';
gliderVariableName = 'cabot_20250723_213_delayed';
gliderFileName = 'cabot_20250723_213_delayed_0660_8583_ae5f.mat';
scheduleFile = fullfile(pwd, '..', '..', 'CTD_processor', 'CTD_processed', 'instrument_schedule_by_cast.csv');
sourceFolder = fullfile(pwd, '..', 'data');
outFolder = fullfile(pwd, '..', 'output', 'compare');
rangeEdges   = [2, 5, 10, 20, 30, 50];   % m from the transducer
% water-depth bands (m): [top, bottom]; Inf = down to the seafloor
depthBands   = [10, 30; 30, 55; 55, Inf; 55, Inf];
bandMargin   = [3, 3, 1, 3];             % m above the seafloor left out, per band
surfaceMin   = 6;                        % m, samples shallower than this are left out
noiseBlockPings = 20;                    % pings per block for the noise estimate
snrThreshold = 10;                       % dB above noise counted as "clear signal"
alphaDefault = [0.034, 0.047, 0.116, 0.267];   % dB/m, used if Output(k).seaAbs is missing
freqLabels = {'130', '200', '455', '769'};
if ~isfolder(outFolder), mkdir(outFolder); end

nR = numel(rangeEdges) - 1;
nB = size(depthBands, 1);
rMid = (rangeEdges(1:end-1) + rangeEdges(2:end)) / 2;
bandName = cell(1, nB);
for b = 1:nB
    if isinf(depthBands(b, 2))
        bandName{b} = sprintf('%g m to %g m above seafloor', depthBands(b, 1), bandMargin(b));
    else
        bandName{b} = sprintf('%g-%g m', depthBands(b, :));
    end
end

%% Glider depth and water depth (time in datenum)
S = load(fullfile(pwd, '..', 'gliderData', gliderFileName), gliderVariableName);
g = S.(gliderVariableName); clear S
tg = g.time(:) / 86400 + datenum(1970, 1, 1);
okD = ~isnan(g.depth(:));
[tD, iu] = unique(tg(okD)); dG = g.depth(okD); dG = dG(iu);
wd = g.m_water_depth(:);
okW = isfinite(wd) & wd > 10 & wd < 150;
[tW, iw] = unique(tg(okW)); wW = wd(okW); wW = wW(iw);
clear g

%% Camera schedule (descents): columns by position, 'end' is a keyword
opts = detectImportOptions(scheduleFile);
vn = opts.VariableNames;
opts = setvartype(opts, vn([2 3 9]), 'char');
sch = readtable(scheduleFile, opts);
isDown = strcmp(sch.(vn{9}), 'down');
fmt = 'yyyy-MM-dd''T''HH:mm:ss''Z''';
cStart = datenum(datetime(sch.(vn{2})(isDown), 'InputFormat', fmt));
cEnd   = datenum(datetime(sch.(vn{3})(isDown), 'InputFormat', fmt));
sg = sch.(vn{14})(isDown);
if iscell(sg), camOn = strcmpi(sg, 'TRUE'); else, camOn = logical(sg); end

%% Read each day and collect per-ping means
nC  = nB * nR;
PM  = cell(1, 4);          % raw: pings x (nB*nR) mean linear sv
PMc = cell(1, 4);          % noise-subtracted
PS  = cell(1, 4);          % share of samples >= snrThreshold above noise
NP  = cell(1, 4);          % pings x nR, any depth (1st-percentile noise reference)
PER = cell(1, 4);          % 1 = day, 2 = night
for k = 1:4
    PM{k} = zeros(0, nC); PMc{k} = zeros(0, nC); PS{k} = zeros(0, nC);
    NP{k} = zeros(0, nR); PER{k} = zeros(0, 1);
end
N0dB = NaN(4, numel(dates));

for dd = 1:numel(dates)
    Parameters.ProcDir = 0;
    Parameters.folderName = sourceFolder;
    Parameters.datafilename = getFilenamesByDate(dates{dd}, sourceFolder);
    Parameters.xmlfilename = fullfile(sourceFolder, xmlFileName);
    Parameters.Salinity = 33;
    Parameters.Bins2Avg = 1;
    Parameters.Time2Avg = 1;
    Parameters.Pressure = 100;
    Parameters.Plot = 0;
    Parameters.Channel = 1;
    Parameters.Value2Plot = 2;
    Parameters.NoiseFloor = [12564; 9660; 13604; 10780];
    Parameters.Orientation = 1;
    Parameters.UseTiltCorr = 0;
    fprintf('Reading %d raw files for %s...\n', size(Parameters.datafilename, 1), dates{dd});
    [Output, ~] = ProcessAZFP(Parameters);

    [t, order] = sort(Output(1).Date(:));
    nPing = numel(t);
    zG  = interp1(tD, dG, t);
    wdp = interp1(tW, wW, t);
    pingCam = NaN(nPing, 1);
    for i = 1:numel(cStart)
        in = t >= cStart(i) & t <= cEnd(i) + 60/86400;
        pingCam(in) = camOn(i);
    end
    localHour = hour(datetime(t, 'ConvertFrom', 'datenum') - hours(3));
    per = zeros(nPing, 1);
    per(localHour >= 7 & localHour <= 20) = 1;
    per(localHour >= 23 | localHour <= 4) = 2;
    use = pingCam == 0 & per > 0 & ~isnan(zG) & ~isnan(wdp);
    nU = sum(use);
    fprintf('  %d non-camera day/night pings used\n', nU);

    for k = 1:4
        R = Output(k).Range(1, :);
        if isfield(Output(k), 'seaAbs') && ~isempty(Output(k).seaAbs)
            alpha = Output(k).seaAbs(1);
        else
            alpha = alphaDefault(k);
        end
        tvgLin = 10.^((20 * log10(R) + 2 * alpha * R) / 10);   % 1 x samples
        lin = 10.^(Output(k).Sv(order(use), :) / 10);
        Z = zG(use) + R;                                         % pings x samples
        water = Z >= surfaceMin & Z <= wdp(use) - min(bandMargin);

        % --- noise level N0 (Sv minus range correction), per day ---
        P = lin ./ tvgLin;
        P(~water) = NaN;
        nBlk = floor(nU / noiseBlockPings);
        blkMin = NaN(nBlk, 1);
        for q = 1:nBlk
            rows = (q - 1) * noiseBlockPings + (1:noiseBlockPings);
            bm = NaN(1, nR);
            for r = 1:nR
                cols = R >= rangeEdges(r) & R < rangeEdges(r + 1);
                if any(cols)
                    blk = P(rows, cols);
                    bm(r) = mean(blk(:), 'omitnan');
                end
            end
            blkMin(q) = min(bm, [], 'omitnan');
        end
        N0 = median(blkMin, 'omitnan');
        N0dB(k, dd) = 10 * log10(N0);
        clear P
        noiseLin = N0 * tvgLin;                                  % 1 x samples
        corr = max(lin - noiseLin, 0);
        clear_ = lin >= noiseLin * 10^(snrThreshold / 10);

        % --- per-ping means in each band x range bin ---
        [pm, pmc, ps] = deal(NaN(nU, nC));
        c = 0;
        for b = 1:nB
            okZ = Z >= surfaceMin & Z <= wdp(use) - bandMargin(b);
            inB = okZ & Z >= depthBands(b, 1) & Z < depthBands(b, 2);
            for r = 1:nR
                c = c + 1;
                cols = R >= rangeEdges(r) & R < rangeEdges(r + 1);
                if ~any(cols), continue; end
                msk = inB(:, cols);
                m = lin(:, cols);   m(~msk) = NaN;   pm(:, c)  = mean(m, 2, 'omitnan');
                m = corr(:, cols);  m(~msk) = NaN;   pmc(:, c) = mean(m, 2, 'omitnan');
                m = double(clear_(:, cols)); m(~msk) = NaN; ps(:, c) = mean(m, 2, 'omitnan');
            end
        end
        np = NaN(nU, nR);
        for r = 1:nR
            cols = R >= rangeEdges(r) & R < rangeEdges(r + 1);
            if ~any(cols), continue; end
            m = lin(:, cols); m(~water(:, cols)) = NaN;
            np(:, r) = mean(m, 2, 'omitnan');
        end
        PM{k} = [PM{k}; pm]; PMc{k} = [PMc{k}; pmc]; PS{k} = [PS{k}; ps];
        NP{k} = [NP{k}; np]; PER{k} = [PER{k}; per(use)];
        clear lin corr clear_ Z water m
    end
    clear Output
end

fprintf('\nNoise level N0 (Sv minus 20logR + 2aR, dB) per day:\n');
for k = 1:4
    fprintf('  %s kHz: %s\n', freqLabels{k}, sprintf('%7.1f', N0dB(k, :)));
end

%% 1st-percentile noise reference (as before) and the subtracted noise curve
noiseRef = NaN(4, nR);
for k = 1:4
    for r = 1:nR
        v = NP{k}(:, r); v = v(~isnan(v) & v > 0);
        if numel(v) >= 50, noiseRef(k, r) = 10 * log10(prctile(v, 1)); end
    end
end

%% Medians and table
[med, medc] = deal(NaN(4, nB, nR, 2));     % freq x band x range x period
[cnt, shr]  = deal(NaN(4, nB, nR, 2));
rows = {};
for k = 1:4
    c = 0;
    for b = 1:nB
        for r = 1:nR
            c = c + 1;
            for p = 1:2
                sel = PER{k} == p & ~isnan(PM{k}(:, c));
                cnt(k, b, r, p) = sum(sel);
                if sum(sel) < 50, continue; end
                v = PM{k}(sel, c);  med(k, b, r, p) = 10 * log10(median(v(v > 0)));
                v = PMc{k}(sel, c); mv = median(v);
                if mv > 0, medc(k, b, r, p) = 10 * log10(mv); end
                shr(k, b, r, p) = 100 * mean(PS{k}(sel, c), 'omitnan');
            end
            rows(end+1, :) = {freqLabels{k}, bandName{b}, ...
                sprintf('%g-%g', rangeEdges(r), rangeEdges(r + 1)), ...
                med(k, b, r, 1), med(k, b, r, 2), medc(k, b, r, 1), medc(k, b, r, 2), ...
                shr(k, b, r, 1), shr(k, b, r, 2), cnt(k, b, r, 1), cnt(k, b, r, 2)}; %#ok<SAGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'freq_kHz', 'depth_band', 'range_m', ...
    'Sv_day_dB', 'Sv_night_dB', 'Sv_day_noisecorr_dB', 'Sv_night_noisecorr_dB', ...
    'pct_day_above_snr', 'pct_night_above_snr', 'n_pings_day', 'n_pings_night'});
writetable(T, fullfile(outFolder, 'avoidance_by_range.csv'));

fprintf('\nMedian Sv (dB) by range from the glider, same water depth, non-camera pings\n');
fprintf('(corr = noise subtracted, NaN = noise-dominated; %%>SNR = share of samples >= %g dB above noise)\n', snrThreshold);
for k = 1:4
    for b = 1:nB
        fprintf('\n %s kHz, water depth %s:\n', freqLabels{k}, bandName{b});
        fprintf('   range (m):    '); fprintf('%9s', strcat(string(rangeEdges(1:end-1)), "-", string(rangeEdges(2:end)))); fprintf('\n');
        fprintf('   day:          '); fprintf('%9.1f', squeeze(med(k, b, :, 1)));  fprintf('\n');
        fprintf('   night:        '); fprintf('%9.1f', squeeze(med(k, b, :, 2)));  fprintf('\n');
        fprintf('   day corr:     '); fprintf('%9.1f', squeeze(medc(k, b, :, 1))); fprintf('\n');
        fprintf('   night corr:   '); fprintf('%9.1f', squeeze(medc(k, b, :, 2))); fprintf('\n');
        fprintf('   night-day corr'); fprintf('%9.1f', squeeze(medc(k, b, :, 2) - medc(k, b, :, 1))); fprintf('\n');
        fprintf('   day %%>SNR:    '); fprintf('%9.0f', squeeze(shr(k, b, :, 1))); fprintf('\n');
        fprintf('   night %%>SNR:  '); fprintf('%9.0f', squeeze(shr(k, b, :, 2))); fprintf('\n');
    end
end

%% Figure: Sv vs range, raw (solid) and noise-subtracted (dashed), day and night
cDay = [0.93, 0.55, 0]; cNight = [0.1, 0.2, 0.55];
fig = figure(1);
try, theme(fig, 'light'); end
set(fig, 'Color', 'w', 'InvertHardcopy', 'off', 'Position', [40, 40, 1600, 1100])
tl = tiledlayout(4, nB, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Median Sv by distance from the glider, same water depth (non-camera dives, 24-29 July 2025)')
subtitle(tl, 'solid = raw, dashed = noise subtracted, grey = noise reference (1st percentile)', 'FontSize', 9)
for k = 1:4
    for b = 1:nB
        nexttile
        hold on
        plot(rMid, squeeze(med(k, b, :, 1)), '-o', 'Color', cDay, 'MarkerFaceColor', cDay, ...
            'LineWidth', 1.5, 'DisplayName', 'day')
        plot(rMid, squeeze(med(k, b, :, 2)), '-o', 'Color', cNight, 'MarkerFaceColor', cNight, ...
            'LineWidth', 1.5, 'DisplayName', 'night')
        plot(rMid, squeeze(medc(k, b, :, 1)), '--s', 'Color', cDay, 'LineWidth', 1.2, ...
            'DisplayName', 'day, noise subtracted')
        plot(rMid, squeeze(medc(k, b, :, 2)), '--s', 'Color', cNight, 'LineWidth', 1.2, ...
            'DisplayName', 'night, noise subtracted')
        plot(rMid, noiseRef(k, :), ':', 'Color', [0.55, 0.55, 0.55], 'LineWidth', 1.2, ...
            'DisplayName', 'noise reference')
        xline(10, ':', 'processing limit', 'Color', [0.4, 0.4, 0.4], 'HandleVisibility', 'off', ...
            'LabelVerticalAlignment', 'bottom', 'FontSize', 7)
        hold off
        grid on; box on
        xlim([0, rangeEdges(end)])
        title(sprintf('%s kHz, %s', freqLabels{k}, bandName{b}), 'FontSize', 9)
        if b == 1, ylabel('Median Sv (dB)'); end
        if k == 4, xlabel('Distance from the transducer (m)'); end
        if k == 1 && b == 1, legend('Location', 'best', 'FontSize', 7); end
    end
end
drawnow
print(fig, '-dpng', fullfile(outFolder, 'figure-avoidance-by-range.png'), '-r300')
fprintf('\nSaved figure-avoidance-by-range.png and avoidance_by_range.csv in %s\n', outFolder);
