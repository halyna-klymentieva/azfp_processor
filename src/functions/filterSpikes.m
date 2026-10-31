function Dives = filterSpikes(Dives, config)
%FILTERSPIKES Remove isolated strong 1 m bins (fish, bubbles, noise spikes)
%   Works on the 1 m bin data of each dive (Dives(j).P(k).sv, pings x depth,
%   linear sv), per frequency, before the noise-floor step and the copepod
%   masks. A bin is removed (set to NaN = no data) when
%     - it is more than config.spikeThreshold dB (15) above the median of its
%       neighbourhood (config.spikeWindowPings pings x config.spikeWindowDepth
%       1 m bins, default 11 x 3; the median is floored at
%       config.spikeFloorDb = -75 dB so isolated weak echoes in empty water
%       are not all flagged), or
%     - its Sv is above config.spikeCeilingDb (-40 dB): far stronger than any
%       copepod layer here (-40 dB at 455 kHz = 1.3 million ind/m3), i.e.
%       fish, bubbles or seafloor.
%   Tested on the merged 2025 data (2026-10-31): removes 0.1-0.4 % of the
%   non-zero bins per frequency, including the largest mid-water spikes.
%   Similar in purpose to the fish / outlier removal in McLennan et al. 2025.

thrDb   = getOpt(config, 'spikeThreshold', 15);
nP      = getOpt(config, 'spikeWindowPings', 11);
nZ      = getOpt(config, 'spikeWindowDepth', 3);
floorDb = getOpt(config, 'spikeFloorDb', -75);
ceilDb  = getOpt(config, 'spikeCeilingDb', -40);
if isempty(thrDb)
    return
end
hp = floor(nP / 2); hz = floor(nZ / 2);
fac = 10^(thrDb / 10); fl = 10^(floorDb / 10); ce = 10^(ceilDb / 10);

fprintf('Removing spikes (> %g dB above the local median, or Sv > %g dB)...\n', thrDb, ceilDb);
nF = numel(Dives(1).P);
nRem = zeros(nF, 1); nNz = zeros(nF, 1);
for j = 1:numel(Dives)
    for k = 1:nF
        x = Dives(j).P(k).sv;                 % pings x depth
        if isempty(x), continue; end
        [m, n] = size(x);
        xp = NaN(m + 2*hp, n + 2*hz);
        xp(hp+1:hp+m, hz+1:hz+n) = x;
        stack = NaN(m, n, nP * nZ);
        c = 0;
        for dp = -hp:hp
            for dz = -hz:hz
                c = c + 1;
                stack(:, :, c) = xp(hp+1+dp : hp+m+dp, hz+1+dz : hz+n+dz);
            end
        end
        ref = max(median(stack, 3, 'omitnan'), fl);
        bad = (x > ref * fac) | (x > ce);     % NaN -> false
        nRem(k) = nRem(k) + sum(bad(:));
        nNz(k) = nNz(k) + sum(x(:) > 0);
        x(bad) = NaN;
        Dives(j).P(k).sv = x;
    end
end
fprintf(' Spike bins removed: %s %% of non-zero bins (%s bins)\n', ...
    mat2str(round(100 * nRem' ./ max(nNz', 1), 2)), mat2str(nRem'));
end

function v = getOpt(config, name, default)
v = default;
if isfield(config, name)
    v = config.(name);
end
end
