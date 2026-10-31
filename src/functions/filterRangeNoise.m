function azfpData = filterRangeNoise(azfpData, config)
%FILTERRANGENOISE Range-dependent noise subtraction and SNR threshold
%   (De Robertis & Higginbottom 2007), applied to the single samples before
%   they are averaged into 1 m depth bins.
%
%   The recorded Sv = signal + noise. Noise is constant in received power,
%   so in Sv units it grows with range like the range correction
%   (TVG = 20 log10 R + 2 alpha R). For each frequency:
%     1. P = sv / tvg (linear), i.e. Sv with the range correction removed;
%     2. blocks of config.noiseBlockPings pings: mean P in 1 m range bins,
%        noise of the block = minimum over the range bins;
%     3. the block values are smoothed with a moving median over
%        config.noiseSmoothBlocks blocks (follows camera on/off dives,
%        which are noisier), and expanded to every ping;
%     4. noise Sv(R) = N0 * tvg is subtracted in linear units;
%     5. samples less than config.snrThreshold dB above the noise are set
%        to 0 (no detectable backscatter; Sv = -Inf).
%   Must run BEFORE the far-field cut (the noise is estimated where the
%   long ranges contain noise only) and before the CTD sound correction
%   (the TVG here is the one the ASL toolbox used).
%
%   Adds per frequency: azfpData(k).NoiseN0dB (1 value per ping, dB re
%   Sv - TVG) and azfpData(k).NoiseSNR (the threshold used).

blockPings = 20;
if isfield(config, 'noiseBlockPings') && ~isempty(config.noiseBlockPings)
    blockPings = config.noiseBlockPings;
end
smoothBlocks = 11;
if isfield(config, 'noiseSmoothBlocks') && ~isempty(config.noiseSmoothBlocks)
    smoothBlocks = config.noiseSmoothBlocks;
end
snr = 10;
if isfield(config, 'snrThreshold') && ~isempty(config.snrThreshold)
    snr = config.snrThreshold;
end

fprintf('Range-dependent noise subtraction (blocks of %d pings, SNR >= %g dB)...\n', blockPings, snr);
tic
for k = 1:numel(azfpData)
    R = azfpData(k).Range(1, :);                      % 1 x samples (m)
    alpha = median(azfpData(k).seaAbs, 'omitnan');    % dB/m used by the toolbox
    tvgLin = 10.^((20 * log10(R) + 2 * alpha * R) / 10);
    nPing = size(azfpData(k).Sv, 1);

    % 1 m range bins
    rBin = floor(R);
    [ub, ~, ib] = unique(rBin);
    nRB = numel(ub);

    % --- noise per block of pings ---
    nBlk = ceil(nPing / blockPings);
    blkN0 = NaN(nBlk, 1);
    for q = 1:nBlk
        rows = (q - 1) * blockPings + 1 : min(q * blockPings, nPing);
        P = 10.^(azfpData(k).Sv(rows, :) / 10) ./ tvgLin;     % pings x samples
        colMean = mean(P, 1, 'omitnan');                       % 1 x samples
        ok = ~isnan(colMean);
        if ~any(ok), continue; end
        s = accumarray(ib(ok), colMean(ok)', [nRB, 1]);
        c = accumarray(ib(ok), 1, [nRB, 1]);
        binMean = s ./ c;
        binMean(c == 0) = NaN;
        blkN0(q) = min(binMean, [], 'omitnan');
    end
    blkN0s = movmedian(blkN0, smoothBlocks, 'omitnan');
    dayN0 = median(blkN0, 'omitnan');
    blkN0s(isnan(blkN0s)) = dayN0;
    N0 = repelem(blkN0s, blockPings);
    N0 = N0(1:nPing);                                           % nPing x 1

    % --- subtract and apply the SNR threshold, in chunks of pings ---
    nWeak = 0; nValid = 0;
    chunk = 2000;
    for r0 = 1:chunk:nPing
        rows = r0 : min(r0 + chunk - 1, nPing);
        lin = 10.^(azfpData(k).Sv(rows, :) / 10);
        noise = N0(rows) * tvgLin;                              % rows x samples
        weak = lin < noise * 10^(snr / 10);
        corr = lin - noise;
        corr(weak) = 0;
        valid = ~isnan(lin);
        nWeak = nWeak + sum(weak(:) & valid(:));
        nValid = nValid + sum(valid(:));
        azfpData(k).Sv(rows, :) = 10 * log10(corr);             % 0 -> -Inf, NaN stays NaN
    end
    azfpData(k).NoiseN0dB = 10 * log10(N0);
    azfpData(k).NoiseSNR = snr;

    q = prctile(10 * log10(blkN0), [5 50 95]);
    fprintf(' %g kHz: noise N0 (Sv - TVG) median %.1f dB (5-95%%: %.1f to %.1f); noise Sv at 5 m / 10 m: %.1f / %.1f dB; below SNR: %.0f %% of samples\n', ...
        azfpData(k).Freq, q(2), q(1), q(3), ...
        q(2) + 20*log10(5) + 2*alpha*5, q(2) + 20*log10(10) + 2*alpha*10, 100 * nWeak / max(nValid, 1));
end
toc
end
