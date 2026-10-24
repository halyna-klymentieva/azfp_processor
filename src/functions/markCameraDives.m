function markCameraDives(ax, config)
%MARKCAMERADIVES Magenta marks under each dive recorded with the camera on
%   Uses config.cameraDives (from getCameraDives). Does nothing if empty.
if ~isfield(config, 'cameraDives') || isempty(config.cameraDives) || ~any(config.cameraDives)
    return
end
d = find(config.cameraDives);
hold(ax, 'on');
plot(ax, d, (config.maxDepth + 3) * ones(size(d)), 'v', 'MarkerSize', 4, ...
    'MarkerFaceColor', [1, 0, 1], 'MarkerEdgeColor', 'none');
hold(ax, 'off');
ylim(ax, [0, config.maxDepth + 6]);
xlabel(ax, 'Dive Number (magenta = camera on)');
end
