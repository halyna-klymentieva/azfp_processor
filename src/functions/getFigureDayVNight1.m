function getFigureDayVNight1(Dives, config, P_indices, freqLabels)
%GETFUGIREDAYVNIGHT1 Summary of this function goes here
f = figure(4);
set(f, 'Position', [100, 100, 1200, 900], 'WindowStyle', 'docked')
for f = 1:3
    j = P_indices(f);

    dayProfiles = [];
    nightProfiles = [];

    for i = 1:length(Dives)
        profile = Dives(i).P(j).mN(:);
        if strcmpi(Dives(i).tod, 'Day')
            dayProfiles = [dayProfiles, profile];
        elseif strcmpi(Dives(i).tod, 'Night')
            nightProfiles = [nightProfiles, profile];
        end
    end

    allVals = [dayProfiles(:); nightProfiles(:)];
    clim = [min(allVals, [], 'omitnan'), min(max(allVals, [], 'omitnan'), 1e4)];


    % Day
    subplot(3, 2, (f - 1)*2+1)
    if ~isempty(dayProfiles)
        imagesc(dayProfiles)
        set(gca, 'YDir', 'reverse')
        yticks(10:10:config.maxDepth)
        yticklabels(string(10:10:config.maxDepth))
        xlabel('Dive Profile')
        ylabel('Depth (m)')
        title(['Daytime - ', freqLabels{f}])
        caxis(clim)
        colormap(hot)
        colorbar
    else
        text(0.5, 0.5, 'No day data', 'HorizontalAlignment', 'center')
        axis off
    end

    % Night
    subplot(3, 2, (f - 1)*2+2)
    if ~isempty(nightProfiles)
        imagesc(nightProfiles)
        set(gca, 'YDir', 'reverse')
        yticks(10:10:config.maxDepth)
        yticklabels(string(10:10:config.maxDepth))
        xlabel('Dive Profile')
        ylabel('Depth (m)')
        title(['Nighttime - ', freqLabels{f}])
        caxis(clim)
        colormap(hot)
        colorbar
    else
        text(0.5, 0.5, 'No night data', 'HorizontalAlignment', 'center')
        axis off
    end
end

% save the figure
if (config.saveWithDateName == 0)
    figure2FileName = "figure2-day-v-night-1-merged.png";
else
    figure2FileName = "figure2-day-v-night-1-" + config.dateOfData + '.png';
end
figure2Path = fullfile(pwd, '..', 'output', figure2FileName);
print(gcf, '-dpng', figure2Path, '-r0')
end