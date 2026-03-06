function  getFigureDayVNight2(Dives, config, P_indices, freqLabels)
%GETFIGUREDAYVNIGHT2 Summary of this function goes here
f = figure(5);
set(f, 'Position', [100, 100, 1200, 900], 'WindowStyle', 'docked')

for f = 1:3
    j = P_indices(f);
    allProfiles = [];
    todLabels = [];
    diveIDs = [];

    for i = 1:length(Dives)
        profile = Dives(i).P(j).mN(:);
        allProfiles = [allProfiles, profile];
        todLabels = [todLabels, string(Dives(i).tod)];
        diveIDs = [diveIDs, i];
    end

    clim = [0, 10000];
    normVals = (allProfiles - clim(1)) / (clim(2) - clim(1));
    normVals = min(max(normVals, 0), 1);

    hotMap = hot(256);
    boneMap = bone(256);

    idxVals = round(normVals*255) + 1;
    [rows, cols] = size(allProfiles);
    RGB = ones(rows, cols, 3);

    for col = 1:cols
        if strcmpi(todLabels(col), 'Day')
            cmap = hotMap;
        else
            cmap = boneMap;
        end
        for ch = 1:3
            RGB(:, col, ch) = cmap(idxVals(:, col), ch);
        end
    end

    ax = subplot(3, 1, f);
    image(RGB)
    set(gca, 'YDir', 'reverse')
    yticks(10:10:config.maxDepth)
    yticklabels(string(10:10:config.maxDepth))
    ylabel('Depth (m)')
    xlabel('Dive Profile')
    title(freqLabels{f})

    tickStep = max(1, floor(cols/10));
    xtickIdx = 1:tickStep:cols;
    xtickLabels = diveIDs(xtickIdx);
    xticks(xtickIdx)
    xticklabels(xtickLabels)

    axPos = get(ax, 'Position');
    cbHot = colorbar('Position', [axPos(1) - 0.05, axPos(2), 0.01, axPos(4)]);
    colormap(cbHot, hotMap);
    cbHot.Ticks = linspace(0, 1, 5);
    cbHot.TickLabels = arrayfun(@(x) sprintf('%.0f', clim(1)+x*(clim(2) - clim(1))), cbHot.Ticks, 'UniformOutput', false);
    cbHot.Label.String = 'Daytime';

    cbBone = colorbar('Position', [axPos(1) + axPos(3) + 0.01, axPos(2), 0.01, axPos(4)]);
    colormap(cbBone, boneMap);
    cbBone.Ticks = linspace(0, 1, 5);
    cbBone.TickLabels = arrayfun(@(x) sprintf('%.0f', clim(1)+x*(clim(2) - clim(1))), cbBone.Ticks, 'UniformOutput', false);
    cbBone.Label.String = 'Nighttime';
end

% save the figure
if (config.saveWithDateName == 0)
    figure2FileName = "figure2-day-v-night-2-merged.png";
else
    figure2FileName = "figure2-day-v-night-2-" + config.dateOfData + '.png';
end
figure2Path = fullfile(pwd, '..', 'output', figure2FileName);
print(gcf, '-dpng', figure2Path, '-r0')
end