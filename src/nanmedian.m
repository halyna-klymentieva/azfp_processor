function m = nanmedian(x, dim)
    if nargin < 2
        m = median(x, 'omitnan');
    else
        m = median(x, dim, 'omitnan');
    end
end