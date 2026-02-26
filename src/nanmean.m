function m = nanmean(x, varargin)
    % This replaces the missing Statistics Toolbox function
    m = mean(x, varargin{:}, 'omitnan');
end

