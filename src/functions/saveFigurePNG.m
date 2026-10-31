function saveFigurePNG(fig, filePath, resolution)
%SAVEFIGUREPNG Save a figure as PNG without stopping the processing
%   If the figure window was closed (or is no longer valid) while the
%   processing was running, a warning is printed and the run continues.
%   resolution: e.g. '-r0' (screen) or '-r300'.
if nargin < 3
    resolution = '-r0';
end
if isempty(fig) || ~isgraphics(fig)
    warning('Figure not saved (window was closed or is invalid): %s', filePath);
    return
end
try
    drawnow
    print(fig, '-dpng', filePath, resolution)
catch err
    warning('Figure not saved (%s): %s', err.message, filePath);
end
end
