%% Camera effect check: AZFP voltage and seafloor echo, camera vs no camera
% For one day, re-reads the raw AZFP files (no filtering, so the seafloor
% echo is still there), tags every ping with the shadowgraph camera state of
% its glider descent, and compares:
%   1. AZFP main battery and transmit battery voltage (per ping);
%   2. seafloor echo strength (peak Sv near the bottom) at each frequency,
%      compared at the same range from the transducer.
% If the camera lowers the AZFP's transmit/receive level, the voltage and the
% seafloor echo should both be lower on camera descents.
% Figures go to output/compare/. Takes a few minutes (raw files are re-read).
% Run with Run / F5.

cd(fileparts(mfilename('fullpath')));
addpath(genpath(pwd));
rehash path % pick up new or changed files (MATLAB can miss changes in OneDrive folders)
clear variables
close all

%% Settings
dateOfData = '25-07-29';
xmlFileName = '25070417.XML';
gliderVariableName = 'cabot_20250723_213_delayed';
gliderFileName = 'cabot_20250723_213_delayed_0660_8583_ae5f.mat';
scheduleFile = fullfile(pwd, '..', '..', 'CTD_processor', 'CTD_processed', 'instrument_schedule_by_cast.csv');
sourceFolder = fullfile(pwd, '..', 'data');
outFolder = fullfile(pwd, '..', 'output', 'compare');
bottomMinDepth = 45;      % m, seafloor searched below this depth (as find_bottom_AM)
bottomMinSv = -50;        % dB, peak must be at least this strong to count as seafloor
rangeClassWidth = 1;      % m, seafloor echoes compared within range classes
freqLabels = {'130', '200', '455', '769'};
if ~isfolder(outFolder), mkdir(outFolder); end

%% 1. Raw AZFP data for the day (same toolbox settings as the main processing)
Parameters.ProcDir = 0;
Parameters.folderName = sourceFolder;
Parameters.datafilename = getFilenamesByDate(dateOfData, sourceFolder);
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
fprintf('Reading %d raw files for %s...\n', size(Parameters.datafilename, 1), dateOfData);
[Output, ~] = ProcessAZFP(Parameters);

[t, order] = sort(Output(1).Date(:));
batMain = Output(1).BatteryMain(order);
batTx = Output(1).BatteryTx(order);
nPing = numel(t);

%% 2. Camera state per ping (from the glider descent the ping belongs to)
% columns: profile_id, start, end, ..., direction (9), ..., sg_on (14).
% 'end' is a MATLAB keyword, so MATLAB renames that column; use positions.
opts = detectImportOptions(scheduleFile);
vn = opts.VariableNames;
opts = setvartype(opts, vn([2 3 9]), 'char');
sch = readtable(scheduleFile, opts);
isDown = strcmp(sch.(vn{9}), 'down');
fmt = 'yyyy-MM-dd''T''HH:mm:ss''Z''';
tStart = datenum(datetime(sch.(vn{2})(isDown), 'InputFormat', fmt));
tEnd = datenum(datetime(sch.(vn{3})(isDown), 'InputFormat', fmt));
sg = sch.(vn{14})(isDown);
if iscell(sg), camOn = strcmpi(sg, 'TRUE'); else, camOn = logical(sg); end

pingCam = NaN(nPing, 1);                 % 1 = camera on, 0 = off, NaN = not in a descent
for i = 1:numel(tStart)
    in = t >= tStart(i) & t <= tEnd(i) + 60/86400;
    pingCam(in) = camOn(i);
end
localHour = hour(datetime(t, 'ConvertFrom', 'datenum') - hours(3));
isNight = localHour >= 23 | localHour <= 4;
isDay = localHour >= 7 & localHour <= 20;
fprintf('%d pings: %d on camera descents, %d on non-camera descents, %d outside descents\n', ...
    nPing, sum(pingCam == 1), sum(pingCam == 0), sum(isnan(pingCam)));

%% 3. Voltage
fprintf('\nAZFP voltage (median, V):\n');
fprintf('                 camera   no camera   difference\n');
groups = {true(nPing, 1), isNight, isDay};
groupNames = {'all day', 'night', 'day'};
for g = 1:3
    c = groups{g} & pingCam == 1; n = groups{g} & pingCam == 0;
    fprintf(' main, %-8s  %7.3f   %7.3f     %+7.3f\n', groupNames{g}, median(batMain(c), 'omitnan'), ...
        median(batMain(n), 'omitnan'), median(batMain(c), 'omitnan') - median(batMain(n), 'omitnan'));
    fprintf(' Tx,   %-8s  %7.3f   %7.3f     %+7.3f\n', groupNames{g}, median(batTx(c), 'omitnan'), ...
        median(batTx(n), 'omitnan'), median(batTx(c), 'omitnan') - median(batTx(n), 'omitnan'));
end

fig1 = figure('Position', [50, 50, 1300, 600], 'Color', 'w');
tl = tiledlayout(2, 1, 'TileSpacing', 'compact');
title(tl, sprintf('AZFP voltage, %s (red = camera descent, blue = no camera, grey = ascent)', dateOfData));
tt = datetime(t, 'ConvertFrom', 'datenum');
for p = 1:2
    nexttile; hold on
    v = batMain; lab = 'Main battery (V)';
    if p == 2, v = batTx; lab = 'Transmit battery (V)'; end
    plot(tt(isnan(pingCam)), v(isnan(pingCam)), '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 3);
    plot(tt(pingCam == 0), v(pingCam == 0), '.', 'Color', [0.1 0.3 0.9], 'MarkerSize', 3);
    plot(tt(pingCam == 1), v(pingCam == 1), '.', 'Color', [0.9 0.1 0.1], 'MarkerSize', 3);
    ylabel(lab); grid on; hold off
end
print(fig1, '-dpng', fullfile(outFolder, ['camera-check-voltage-' dateOfData '.png']), '-r0');

%% 4. Seafloor echo strength
[gliderDepth, gliderTime] = loadGliderData(fullfile(pwd, '..', 'gliderData', gliderFileName), gliderVariableName);
[gT, iu] = unique(gliderTime(:));
gD = gliderDepth(:);
zGlider = interp1(gT, gD(iu), t);

fprintf('\nSeafloor echo (peak Sv deeper than %d m and >= %d dB):\n', bottomMinDepth, bottomMinSv);
fprintf(' freq   n camera  n no-cam   median camera  median no-cam   range-matched difference (camera - no camera)\n');
fig2 = figure('Position', [50, 50, 1300, 800], 'Color', 'w');
tl2 = tiledlayout(2, 2, 'TileSpacing', 'compact');
title(tl2, sprintf('Seafloor echo vs range, %s (red = camera, blue = no camera)', dateOfData));
for k = 1:4
    sv = Output(k).Sv(order, :);
    R = Output(k).Range(1, :);
    Z = zGlider + R;
    sv(~(Z > bottomMinDepth & R > 2)) = -Inf;
    [pk, ib] = max(sv, [], 2);
    rBot = R(ib)';
    isBot = pk >= bottomMinSv & ~isnan(pingCam) & ~isnan(zGlider);
    c = isBot & pingCam == 1; n = isBot & pingCam == 0;

    % compare within range classes so the different ranges do not bias it
    edges = 0:rangeClassWidth:ceil(max(R));
    cls = discretize(rBot, edges);
    d = []; w = [];
    for j = 1:numel(edges) - 1
        cj = c & cls == j; nj = n & cls == j;
        if sum(cj) >= 20 && sum(nj) >= 20
            d(end+1) = median(pk(cj)) - median(pk(nj)); %#ok<SAGROW>
            w(end+1) = min(sum(cj), sum(nj));            %#ok<SAGROW>
        end
    end
    if isempty(d)
        dMean = NaN;
    else
        dMean = sum(d .* w) / sum(w);
    end
    fprintf(' %3s   %8d  %8d   %11.1f   %12.1f        %+6.2f dB (%d range classes)\n', freqLabels{k}, ...
        sum(c), sum(n), median(pk(c)), median(pk(n)), dMean, numel(d));

    nexttile; hold on
    plot(rBot(n), pk(n), '.', 'Color', [0.1 0.3 0.9], 'MarkerSize', 4);
    plot(rBot(c), pk(c), '.', 'Color', [0.9 0.1 0.1], 'MarkerSize', 4);
    xlabel('Range from transducer to seafloor peak (m)'); ylabel('Peak Sv (dB)');
    title(sprintf('%s kHz: camera - no camera = %+.1f dB', freqLabels{k}, dMean)); grid on; hold off
end
print(fig2, '-dpng', fullfile(outFolder, ['camera-check-seafloor-' dateOfData '.png']), '-r0');
fprintf('\nFigures saved in %s\n', outFolder);

%% 5. Range test: camera - no camera by distance from the transducer
% Water column only (5-45 m depth). How to read it:
%   ~0 near, growing positive with range  -> added noise on camera dives
%   same value at every range             -> sensitivity (gain) change
%   most negative near, ~0 far away       -> fewer/weaker scatterers near the glider
fprintf('\nRange test (camera - no camera, dB, water column 5-45 m):\n');
edges = [2 4 8 15 30 60];
for per = 1:2
    if per == 1, sel = isNight; name = 'night'; else, sel = isDay; name = 'day'; end
    for k = 1:4
        sv = Output(k).Sv(order, :); R = Output(k).Range(1, :); Z = zGlider + R;
        lin = 10.^(sv / 10); lin(Z > 45 | Z < 5) = NaN;
        fprintf(' %-5s %3s kHz:', name, freqLabels{k});
        for j = 1:numel(edges) - 1
            cols = R >= edges(j) & R < edges(j+1);
            if ~any(cols), continue; end
            pm = 10*log10(mean(lin(:, cols), 2, 'omitnan'));
            d = median(pm(sel & pingCam == 1), 'omitnan') - median(pm(sel & pingCam == 0), 'omitnan');
            fprintf('  %g-%g m: %+.1f', edges(j), edges(j+1), d);
        end
        fprintf('\n');
    end
end

%% 6. Near-glider test: camera - no camera at 2-8 m range, by glider depth
% If animals avoid the camera, the deficit should appear wherever the glider is
% (at night), and possibly also deep down by day where it is dark.
fprintf('\nNear-glider test (camera - no camera, dB, 2-8 m from transducer, by glider depth):\n');
bands = [0 10 20 30 40 50];
for per = 1:2
    if per == 1, sel = isNight; name = 'night'; else, sel = isDay; name = 'day'; end
    for k = 1:4
        sv = Output(k).Sv(order, :); R = Output(k).Range(1, :);
        cols = R >= 2 & R < 8;
        lin = 10.^(sv(:, cols) / 10);
        lin((zGlider + R(cols)) > 60) = NaN;          % keep clear of the seafloor
        pm = 10*log10(mean(lin, 2, 'omitnan'));
        fprintf(' %-5s %3s kHz:', name, freqLabels{k});
        for b = 1:numel(bands) - 1
            g = sel & zGlider >= bands(b) & zGlider < bands(b+1);
            d = median(pm(g & pingCam == 1), 'omitnan') - median(pm(g & pingCam == 0), 'omitnan');
            fprintf('  %d-%d m: %+.1f', bands(b), bands(b+1), d);
        end
        fprintf('\n');
    end
end
