function getFigureDBDiffAllFrqs(Sv_200_130, Sv_455_200, Sv_769_400, config)
%GETFIGUREDBDIFFALLFRQS Plot dB differences for all frequencies; Save a plot
figure(2);
clf

diffLabels = {'200-130 kHz', '455-200 kHz', '769-455 kHz'};
dbDiff_all = {Sv_200_130, Sv_455_200, Sv_769_400};
caxisVals = [3, 9];

for i = 1:3
    subplot(2, 2, i)
    imagesc(dbDiff_all{i}', 'AlphaData', ~isnan(dbDiff_all{i}'))
    colormap('jet')
    clim(caxisVals)
    xlabel('Ping Number')
    ylabel('Depth (m)')
    title(diffLabels{i})
end

% Shared colorbar
h = colorbar;
ylabel(h, 'dB Difference')
h.Position(4) = 0.65;
h.Position(1) = .94 - h.Position(3);
h.Position(2) = 0.5 - h.Position(4) / 2;

% save the figure
if (config.saveWithDateName == 0)
    figure2FileName = "figure2-db-diff-merged.png";
else
    figure2FileName = "figure2-db-diff-" + config.dateOfData + '.png';
end
figure2Path = fullfile(pwd, '..', 'output', figure2FileName);
print(gcf, '-dpng', figure2Path, '-r0')
end
