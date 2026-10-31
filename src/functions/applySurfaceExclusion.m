function Dives = applySurfaceExclusion(Dives, config)
%APPLYSURFACEEXCLUSION Set the near-surface 1 m depth bins to NaN in every dive
%   Bins 1..config.surfaceExclusion (bin z = depth rounded to the nearest
%   metre, so 5 removes everything shallower than 5.5 m) are set to NaN for
%   all frequencies, before the noise floor, dB differences, masks and
%   densities are calculated. Reason: strong broadband echoes near the
%   surface (bubbles / glider at the surface; few pings per bin) that fade by
%   ~6 m in the 2025 data (see processing log, 2026-10-26).
%   config.surfaceExclusion = [] or 0 -> nothing removed.
if ~isfield(config, 'surfaceExclusion') || isempty(config.surfaceExclusion) || config.surfaceExclusion <= 0
    return
end
cut = round(config.surfaceExclusion);
for j = 1:numel(Dives)
    for k = 1:numel(Dives(j).P)
        sv = Dives(j).P(k).sv;
        n = min(cut, size(sv, 2));
        sv(:, 1:n) = NaN;                       % pings x depth bins
        Dives(j).P(k).sv = sv;
    end
end
fprintf('Surface exclusion: depth bins 1-%d m set to NaN (data shallower than %.1f m removed).\n', cut, cut + 0.5);
end
