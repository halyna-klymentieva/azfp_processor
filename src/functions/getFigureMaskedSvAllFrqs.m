function getFigureMaskedSvAllFrqs(masked_130, masked_200, masked_455, ...
    masked_769, config)
%GETFIGUREMASKEDSVALLFRQS Plot masked Sv for all frequencies; save a plot
fig = figure(3);
clf

freqLabels = {'130 kHz', '200 kHz', '455 kHz', '769 kHz'};
Sv_all = {masked_130, masked_200, masked_455, masked_769};

for i = 1:4
    subplot(2, 2, i)
    imagesc(Sv_all{i}', 'AlphaData', ~isnan(Sv_all{i}'))
    colormap('jet')
    if isfield(config, 'svLimits') && numel(config.svLimits) == 2, clim(config.svLimits), else, clim([-100, -50]), end % same fixed scale for all panels (one shared colorbar), as in figure 1
    ylim([0, config.maxDepth])
    xlabel('Dive Number')
    markCameraDives(gca, config)
    ylabel('Depth (m)')
    title(freqLabels{i})
end

% Shared colorbar
h = colorbar;
ylabel(h, 'Volume Backscatter (dB re 1 m^{-1})')
h.Position(4) = 0.65;
h.Position(1) = .94 - h.Position(3);
h.Position(2) = 0.5 - h.Position(4) / 2;

% save the figure
if (config.saveWithDateName == 0)
    figure3FileName = "figure3-m-masked-merged.png";
else
    figure3FileName = "figure3-m-masked-" + config.dateOfData + '.png';
end
figure3Path = fullfile(pwd, '..', 'output', figure3FileName);
drawnow
saveFigurePNG(fig, figure3Path, '-r0')
end