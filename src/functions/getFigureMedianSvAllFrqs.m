function getFigureMedianSvAllFrqs(Sv_130, Sv_200, Sv_455, Sv_769, ...
    bottomDepth, config)
%GETFIGUREMEDIANSVALLFRQS Plot median Sv for all frequencies; Save a plot
figure(1)
clf

freqLabels = {'130 kHz', '200 kHz', '455 kHz', '769 kHz'};
Sv_all = {Sv_130, Sv_200, Sv_455, Sv_769};
caxisVals = [-110, -60];

for i = 1:4
    subplot(2, 2, i)
    imagesc(Sv_all{i}', 'AlphaData', ~isnan(Sv_all{i}'))
    colormap('jet')
    clim(caxisVals)
    ylim([0, config.maxDepth])
    xlabel('Ping Number')
    ylabel('Depth (m)')
    title(freqLabels{i})

    hold on
    plot(1:length(bottomDepth), bottomDepth, 'k-', 'LineWidth', 1, 'Color', [0, 0, 0, 0.3])
    hold off
end

% Shared colorbar
h = colorbar;
ylabel(h, 'Volume Backscatter (dB re 1 m^{-1})')
h.Position(4) = 0.65;
h.Position(1) = .94 - h.Position(3);
h.Position(2) = 0.5 - h.Position(4) / 2;

% save the figure
if (config.saveWithDateName == 0) 
    figure1MedianSVFileName = "figure1-median-sv-merged.png";
else 
    figure1MedianSVFileName = "figure1-median-sv-" + config.dateOfData + '.png';
end
figure1MedianSVPath = fullfile(pwd, '..', 'output', figure1MedianSVFileName);
print(gcf, '-dpng', figure1MedianSVPath, '-r0')
end