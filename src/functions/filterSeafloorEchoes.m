function output = filterSeafloorEchoes(azfpData, dives, config)
%FILTERSEAFLOORECHOES Remove seafloor echoes (rewritten 2026-10-31)
%   1. Seafloor onset per ping and frequency: the shallowest sample deeper
%      than config.seafloorMinDepth (45 m) with Sv >= config.seafloorThreshold
%      (-35 dB). Seafloor echoes here are about -10 to -20 dB; the densest
%      copepod layers reach about -45 dB (single samples), so -35 dB separates
%      them. (The previous detector, find_bottom_AM, took the shallowest peak
%      >= -55 dB, which in dense near-bottom layers was the layer, not the
%      seafloor: the seafloor estimate was then ~10 m too shallow and the real
%      seafloor echo was not removed - up to 1e9 "ind/m3" in the deepest bin.)
%   2. Seafloor depth per ping = shallowest onset over the four frequencies,
%      moving median over config.seafloorSmoothPings pings (101), kept only
%      where a detection exists within the window -> azfpData(1).SeafloorDepth.
%   3. Cut depth per ping and frequency = the frequency's own onset if it is
%      within config.seafloorMaxJump (3 m) of the smoothed seafloor depth
%      (otherwise, e.g. a fish, the smoothed depth); everything deeper than
%      (cut depth - config.seafloorMargin) is removed.

if nargin < 3, config = struct(); end
thr      = getOpt(config, 'seafloorThreshold', -35);
minDep   = getOpt(config, 'seafloorMinDepth', 45);
margin   = getOpt(config, 'seafloorMargin', 0.5);
smoothP  = getOpt(config, 'seafloorSmoothPings', 101);
maxJump  = getOpt(config, 'seafloorMaxJump', 3);

fprintf('Filtering seafloor echoes (onset >= %g dB below %g m, margin %.1f m)...\n', thr, minDep, margin);
tic
nF = length(azfpData);
nPing = size(azfpData(1).Sv, 1);

%% 1. onset per ping and frequency
onset = NaN(nF, nPing);
for k = 1:nF
    dep = azfpData(k).Depth;                       % pings x samples
    hit = azfpData(k).Sv >= thr & dep > minDep;    % NaN/-Inf -> false
    [anyHit, first] = max(hit, [], 2);
    rows = find(anyHit);
    onset(k, rows) = dep(sub2ind(size(dep), rows, first(rows)));
end

%% 2. smoothed seafloor depth
zRaw = min(onset, [], 1, 'omitnan');
zSea = movmedian(zRaw, smoothP, 'omitnan');
hasDet = movsum(~isnan(zRaw), smoothP) > 0;
zSea(~hasDet) = NaN;
azfpData(1).SeafloorDepth = zSea(:);

%% 3. remove
nRemoved = zeros(nF, 1); nOwn = zeros(nF, 1); nOdd = zeros(nF, 1);
for k = 1:nF
    cutDep = zSea;
    own = onset(k, :);
    useOwn = ~isnan(own) & (isnan(zSea) | abs(own - zSea) <= maxJump);
    cutDep(useOwn) = own(useOwn);
    nOwn(k) = sum(useOwn);
    nOdd(k) = sum(~isnan(own) & ~useOwn);
    cut = azfpData(k).Depth >= (cutDep(:) - margin);   % NaN -> false
    nRemoved(k) = sum(cut(:) & ~isnan(azfpData(k).Sv(:)));
    azfpData(k).Sv(cut) = NaN;
end
toc
fprintf(' Seafloor depth known for %.0f %% of pings (median %.1f m).\n', 100 * mean(~isnan(zSea)), median(zSea, 'omitnan'));
fprintf(' Pings with own seafloor onset used: %s; strong echoes > %g m from the seafloor ignored (fish?): %s\n', ...
    mat2str(nOwn'), maxJump, mat2str(nOdd'));
fprintf(' Samples removed: %s per frequency\n', mat2str(nRemoved'));

output = azfpData;
end

function v = getOpt(config, name, default)
v = default;
if isfield(config, name) && ~isempty(config.(name))
    v = config.(name);
end
end
