function drawAndSaveFigures(Dives, bottomDepth, config)
%% Draw and save figures
fprintf('Building figuires...\n')
tic
% make matrices of each Sv variable
chartData = getChartsData(Dives);
%  Plot median Sv for all frequencies; Save a plot
getFigureMedianSvAllFrqs(chartData.Sv_130, chartData.Sv_200, ...
    chartData.Sv_455, chartData.Sv_769, bottomDepth, config)

% Plot dB differences for all frequencies; Save a plot
getFigureDBDiffAllFrqs(chartData.Sv_200_130, chartData.Sv_455_200, ...
    chartData.Sv_769_400, config)

% Plot masked Sv for all frequencies; save a plot
getFigureMaskedSvAllFrqs(chartData.masked_130, chartData.masked_200, ...
    chartData.masked_455, chartData.masked_769, config)

%% save variables

% filename = strcat("d:/AZFP_GSL_2025/AZFP_processed_data/", "July25_Dive.mat");
% save(filename, 'Dive', 'chartData', 'bottomDepth', '-v7.3');

%% Plot numeric abundances for all frequencies
%Approximate civil twilight in July 2024: 7AM-8PM (day), 11PM-4AM (night)

%Numeric abundances were calculated from linearized masked Sv values divided
%by corresponding sigma_bs and averaged by dive profile

P_indices = [2, 3, 4];
freqLabels = {'200 kHz', '455 kHz', '769 kHz'};

% figure('Position', [100, 100, 1200, 900])
% depths = 1:config.maxDepth;
%
% for f = 1:3
%     j = P_indices(f);
%
%     dayProfiles = [];
%     nightProfiles = [];
%
%     for i = 1:length(Dives)
%         profile = Dives(i).P(j).mN(:);
%         if strcmpi(Dives(i).tod, 'Day')
%             dayProfiles = [dayProfiles, profile];
%         elseif strcmpi(Dives(i).tod, 'Night')
%             nightProfiles = [nightProfiles, profile];
%         end
%     end
%
%     allVals = [dayProfiles(:); nightProfiles(:)];
%     clim = [min(allVals,[],'omitnan'), max(allVals,[],'omitnan')];
%
%     % Day
%     subplot(3, 2, (f - 1)*2+1)
%     if ~isempty(dayProfiles)
%         imagesc(dayProfiles)
%         set(gca, 'YDir', 'reverse')
%         yticks(10:10:config.maxDepth)
%         yticklabels(string(10:10:config.maxDepth))
%         xlabel('Dive Profile')
%         ylabel('Depth (m)')
%         title(['Daytime - ', freqLabels{f}])
%         caxis(clim)
%         colormap(hot)
%         colorbar
%     else
%         text(0.5, 0.5, 'No day data', 'HorizontalAlignment', 'center')
%         axis off
%     end
%
%     % Night
%     subplot(3, 2, (f - 1)*2+2)
%     if ~isempty(nightProfiles)
%         imagesc(nightProfiles)
%         set(gca, 'YDir', 'reverse')
%         yticks(10:10:config.maxDepth)
%         yticklabels(string(10:10:config.maxDepth))
%         xlabel('Dive Profile')
%         ylabel('Depth (m)')
%         title(['Nighttime - ', freqLabels{f}])
%         caxis(clim)
%         colormap(hot)
%         colorbar
%     else
%         text(0.5, 0.5, 'No night data', 'HorizontalAlignment', 'center')
%         axis off
%     end
% end

getFigureDayVNight1(Dives, config, P_indices, freqLabels)

getFigureDayVNight2(Dives, config, P_indices, freqLabels)

% figure;
% tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact')
%
% for f = 1:3
%     j = P_indices(f);
%     dayMeans = [];
%     nightMeans = [];
%
%     for i = 1:length(Dives)
%         profile = Dives(i).P(j).mN(:);
%         meanVal = mean(profile, 'omitnan');
%         if strcmpi(Dives(i).tod, 'Day')
%             dayMeans(end+1) = meanVal;
%         elseif strcmpi(Dives(i).tod, 'Night')
%             nightMeans(end+1) = meanVal;
%         end
%     end
%
%     nexttile
%     hold on
%
%     rng default
%     x1 = normrnd(5, 1, 100, 1);
%     x2 = normrnd(6, 1, 100, 1);
%     figure
%     boxplot([x1, x2], 'Notch', 'on', 'Labels', {'mu = 5', 'mu = 6'})
%     title('Compare Random Data from Different Distributions')
%
%     x = [dayMeans, nightMeans]; % numeric column
%     g = [repmat({'Day'}, numel(dayMeans), 1); repmat({'Night'}, numel(nightMeans), 1)];
%     boxplot(x, g, 'Colors', 'k', 'Symbol', '');
%
%     % boxplot([dayMeans, nightMeans], ...
%     %     [repmat({'Day'}, 1, length(dayMeans)), repmat({'Night'}, 1, length(nightMeans))], ...
%     %     'Colors', 'k', 'Symbol', '')
%
%     xDay = ones(size(dayMeans)) + 0.1 * randn(size(dayMeans));
%     xNight = 2 * ones(size(nightMeans)) + 0.1 * randn(size(nightMeans));
%     scatter(xDay, dayMeans, 30, [0.6, 0.8, 1], 'filled', 'MarkerFaceAlpha', 0.6)
%     scatter(xNight, nightMeans, 30, [0, 0, 0.5], 'filled', 'MarkerFaceAlpha', 0.6)
%
%     ylim([0, 10000])
%
%     set(gca, 'XTickLabel', {'Day', 'Night'})
%     ylabel('Mean mN value')
%     title(freqLabels{f})
%     hold off
% end
%
% nexttile
% axis off
% legend({'Day points', 'Night points'}, 'Location', 'best')

% save the figure
%filename = strcat("C:/Users/Andrea/Documents/AMesquita2025/UNB/dataAnalysis/preyData/processedAZFP/updatedCode/2024/mMasked2024_1_4.png");
%print(gcf,'-dpng',filename,'-r0')
toc
end
