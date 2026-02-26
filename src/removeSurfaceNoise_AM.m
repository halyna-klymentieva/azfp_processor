%% Remove surface noise causing streaks in the data

function Sv_clean = removeSurfaceNoise_AM(Sv, binSize)
% Remove strong surface echoes from the top of the echogram

    if nargin < 2
        binSize = 0.3;
    end

    maxDepth = 5;
    surfBufBins = 10; % remove ± this many bins around surface peak
    maxSurfBin = round(maxDepth / binSize); % max bin index to inspect

    Sv_clean = Sv;

    for ping = 1:size(Sv,2)
        sv_col = Sv(1:maxSurfBin, ping);

        if all(isnan(sv_col)), continue; end

        % Find the peak in top 5m (most reflective bin)
        [~, imax] = max(sv_col);

        % Define bounds to mask
        i1 = max(1, imax - surfBufBins);
        i2 = min(maxSurfBin, imax + surfBufBins);

        % Set values to NaN
        Sv_clean(i1:i2, ping) = NaN;
    end
end