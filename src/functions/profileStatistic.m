function m = profileStatistic(x, config)
%PROFILESTATISTIC Per-depth profile of a dive (over pings, i.e. down columns)
%   config.profileStatistic = 'median' (default, original) or 'mean'.
%   NaN (no data) is ignored. Use linear values (sv, density), not dB.
%   With many zeros (no detectable backscatter / no copepods), the median is
%   often 0; the mean gives the average per volume of water.
stat = 'median';
if nargin > 1 && isfield(config, 'profileStatistic') && ~isempty(config.profileStatistic)
    stat = lower(config.profileStatistic);
end
switch stat
    case 'mean'
        m = mean(x, 1, 'omitnan');
    case 'median'
        m = median(x, 1, 'omitnan');
    otherwise
        error('config.profileStatistic must be ''median'' or ''mean'' (got ''%s'').', stat);
end
end
