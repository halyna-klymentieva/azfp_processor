function azfpData = applyCTDSoundCorrection(azfpData, gliderDepth, ctdFile)
%APPLYCTDSOUNDCORRECTION Recalculate range, Sv and depth with CTD sound speed and absorption
%   The ASL toolbox computes Sv with ONE sound speed and ONE absorption
%   coefficient per channel, from a fixed temperature (5 degC when the XML
%   has no temperature coefficients), salinity 33 and pressure 100 dbar.
%   This function replaces those with values from the glider CTD:
%     - each ping is matched (by time) to the CTD cast it was recorded in;
%     - for every sample, sound speed (harmonic mean) and absorption
%       (arithmetic mean) are averaged along the acoustic path, from the
%       glider depth down to the sample depth;
%     - range is rescaled by c_path / c0, and Sv is corrected:
%         Sv_new = Sv + 20log10(R_new/R0) + 2*(alpha_path*R_new - alpha0*R0)
%                     - 10log10(c_local/c0)
%     - Depth is recalculated as glider depth + corrected range.
%   Above the shallowest CTD bin of a cast (3 dbar) the shallowest value is
%   used; below the deepest bin the deepest value is used.
%   TS is not corrected (it is not used later in the processing).
%
%   gliderDepth : glider depth per ping (m), same number of rows as Sv
%   ctdFile     : ctd_binned_1dbar_*.csv from CTD_processor

fprintf('Applying CTD sound speed and absorption...\n');
tic

%% Nominal values the ASL toolbox used
c0 = median(azfpData(1).SoundSpeed, 'omitnan');
if (max(azfpData(1).SoundSpeed) - min(azfpData(1).SoundSpeed)) > 0.5
    warning('Toolbox sound speed varies between files (%.1f-%.1f m/s); using the median %.1f m/s.', ...
        min(azfpData(1).SoundSpeed), max(azfpData(1).SoundSpeed), c0);
end
nCh = numel(azfpData);
alpha0 = zeros(1, nCh);
freq = zeros(1, nCh);
for k = 1:nCh
    alpha0(k) = median(azfpData(k).seaAbs, 'omitnan');
    freq(k) = azfpData(k).Freq;
end
fprintf(' Toolbox values: c0 = %.1f m/s, alpha0 = %s dB/m\n', c0, mat2str(round(alpha0, 4)));

%% Load CTD casts
opts = detectImportOptions(ctdFile);
opts = setvartype(opts, {'time', 'direction'}, 'char');
ctd = readtable(ctdFile, opts);
ctdTime = datenum(datetime(ctd.time, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z'''));

pingTime = azfpData(1).Date(:);
gliderDepth = gliderDepth(:);
nPing = numel(pingTime);

% keep only casts near this day's pings
keep = ctdTime >= min(pingTime) - 1/24 & ctdTime <= max(pingTime) + 1/24;
ctd = ctd(keep, :);
ctdTime = ctdTime(keep);
[profIds, ~, g] = unique(ctd.profile_id);
nProf = numel(profIds);
if nProf == 0
    warning('No CTD casts overlap these pings - CTD correction skipped.');
    return
end
tStart = accumarray(g, ctdTime, [], @min);
tEnd = accumarray(g, ctdTime, [], @max);

%% Match each ping to a cast (the cast whose time span contains it, else the nearest)
pingProf = zeros(nPing, 1);
gap = inf(nPing, 1);
for p = 1:nProf
    d = max(max(tStart(p) - pingTime, pingTime - tEnd(p)), 0); % 0 inside the cast
    better = d < gap;
    gap(better) = d(better);
    pingProf(better) = p;
end
maxGap = 30 / 1440; % 30 min
pingProf(gap > maxGap | isnan(gliderDepth)) = 0;
fprintf(' %d of %d pings matched to %d CTD casts (median gap %.1f s).\n', ...
    sum(pingProf > 0), nPing, nProf, 86400 * median(gap(pingProf > 0)));

%% Depth grid for path integrals
zMax = ceil(max(ctd.depth_ctd)) + 40;
zg = (0:0.5:zMax)';

sumDsv = zeros(1, nCh); minDsv = inf(1, nCh); maxDsv = -inf(1, nCh); nDsv = zeros(1, nCh);
sumDz = zeros(1, nCh);

for p = 1:nProf
    rows = find(pingProf == p);
    if isempty(rows)
        continue
    end
    sel = g == p;
    [zc, iu] = unique(ctd.depth_ctd(sel));
    Tc = ctd.T_corr(sel); Tc = Tc(iu);
    Sc = ctd.SP_tm(sel); Sc = Sc(iu);
    cc = ctd.sound_speed(sel); cc = cc(iu);
    ok = ~isnan(zc) & ~isnan(Tc) & ~isnan(Sc) & ~isnan(cc);
    zc = zc(ok); Tc = Tc(ok); Sc = Sc(ok); cc = cc(ok);
    if numel(zc) < 2
        continue
    end
    T = gridProfile(zc, Tc, zg);
    S = gridProfile(zc, Sc, zg);
    c = gridProfile(zc, cc, zg);
    invC = cumtrapz(zg, 1 ./ c);              % integral of 1/c dz

    zGlider = max(gliderDepth(rows), 0);
    for k = 1:nCh
        a = computeAbsCTD(T, zg, S, freq(k));  % dB/m along depth (P ~ depth in dbar)
        intA = cumtrapz(zg, a);

        R0 = azfpData(k).Range(1, :);         % nominal range (m), 1 x nBins
        Z0 = min(zGlider + R0, zMax);         % nPing x nBins
        dz = Z0 - zGlider;
        cPath = dz ./ (interp1(zg, invC, Z0) - interp1(zg, invC, zGlider));
        aPath = (interp1(zg, intA, Z0) - interp1(zg, intA, zGlider)) ./ dz;

        Rn = R0 .* cPath / c0;
        Zn = zGlider + Rn;
        cLoc = interp1(zg, c, min(Zn, zMax));

        dSv = 20 * log10(Rn ./ R0) + 2 * (aPath .* Rn - alpha0(k) .* R0) - 10 * log10(cLoc / c0);
        good = ~isnan(dSv);
        sv = azfpData(k).Sv(rows, :);
        sv(good) = sv(good) + dSv(good);
        azfpData(k).Sv(rows, :) = sv;
        dep = azfpData(k).Depth(rows, :);
        dep(good) = Zn(good);
        azfpData(k).Depth(rows, :) = dep;

        sumDsv(k) = sumDsv(k) + sum(dSv(good));
        nDsv(k) = nDsv(k) + sum(good(:));
        minDsv(k) = min(minDsv(k), min(dSv(good)));
        maxDsv(k) = max(maxDsv(k), max(dSv(good)));
        sumDz(k) = sumDz(k) + sum(abs(Zn(good) - Z0(good)));
    end
end

for k = 1:nCh
    fprintf(' %3d kHz: Sv change mean %+.2f dB (range %+.2f to %+.2f), mean depth shift %.2f m\n', ...
        freq(k), sumDsv(k) / max(nDsv(k), 1), minDsv(k), maxDsv(k), sumDz(k) / max(nDsv(k), 1));
end
azfpData(1).CTDProfileId = NaN(nPing, 1);
azfpData(1).CTDProfileId(pingProf > 0) = profIds(pingProf(pingProf > 0));
azfpData(1).CTDSoundCorrection = true;
toc
end

function v = gridProfile(zc, vc, zg)
% linear interpolation onto zg; constant (nearest) above and below the cast
v = interp1(zc, vc, zg, 'linear');
v(zg < zc(1)) = vc(1);
v(zg > zc(end)) = vc(end);
end

function seaAbs = computeAbsCTD(T, P, S, Freq)
% Same absorption formula as computeAbs in the ASL toolbox (LoadAZFP.m),
% vectorised over T, P (dbar) and S. Freq in kHz, result in dB/m.
T_K = T + 273.0;
f1 = 1320.0 .* T_K .* exp(-1700 ./ T_K);
f2 = (1.55e7) .* T_K .* exp(-3052 ./ T_K);
k = 1 + P / 10.0;
a = (8.95e-8) .* (1 + T .* ((2.29e-2) - (5.08e-4) .* T));
b = (S / 35.0) .* (4.88e-7) .* (1 + 0.0134 .* T) .* (1 - 0.00103 .* k + (3.7e-7) .* (k .* k));
c = (4.86e-13) .* (1 + T .* ((-0.042) + T .* ((8.53e-4) - T .* 6.23e-6))) .* (1 + k .* (-(3.84e-4) + k .* 7.57e-8));
f = Freq * 1000;
seaAbs = (a .* f1 .* f.^2) ./ (f1.^2 + f.^2) + (b .* f2 .* f.^2) ./ (f2.^2 + f.^2) + c .* f.^2;
end
