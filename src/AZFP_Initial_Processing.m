%Version 1 created in 2020 by Kim Davies and Delphine Mossman

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update log:

% Jan 7, 2021 AZFP Noise Floor Estimate
% Scott Loranger

% 31 August 2021 Far Field Noise Floor removal code
% Modified Nov 19, 2021
% Delphine Mossman

% Version 2 created by Kim Davies
% Oct 4, 2024
% Code clean up and made some pieces more efficient for processing longer
% deployments.

% Halyna's version
% Feb 19, 2026
%Code clean up; include the code to save output from each day

%% Prepare your workspace and file directories
addpath(genpath(pwd));
clc
clear variables
close all

%% User-defined config variables
config.dateOfData = '25-07-29';
config.xmlFileName = '25072317.XML';
config.gliderVariableName = 'cabot_20250723_213_delayed';
config.gliderFileName = 'cabot_20250723_213_delayed_0660_8583_ae5f.mat';
config.calibrationOffsets = [-0.72, -4.60, -4.63, -0.20];
% Far field per frequency cut-off range, see filterFarField function
config.farFieldCutOffRange = [30, 10, 10, 5];
config.dataCutOffDepth = 67.5;
config.surfaceNoiseDepthMin = 10;
config.surfaceNoiseDepthMax = 14;
config.maxDepth = 105;

%% System config variables
config.azfpDataCacheFilename = "data-" + config.dateOfData + '.mat';
config.divesDataCacheFilename = "dives-" + config.dateOfData + ".mat";
config.sourceFolder = fullfile(pwd, '..', 'data');
config.outputFolder = fullfile(pwd, '..', 'output');
config.sourceFileNames = getFilenamesByDate(config.dateOfData, config.sourceFolder);
fprintf('For date %s files count: %d.\n', config.dateOfData, length(config.sourceFileNames));
config.azfpDataCachePath = fullfile(config.outputFolder, config.azfpDataCacheFilename);
config.divesDataCachePath = fullfile(config.outputFolder, config.divesDataCacheFilename);
config.gliderFileFullName = fullfile(pwd, '..', 'gliderData', config.gliderFileName);

%% Getting output of AZFP raw data processing
% if current date is not processed yet it is processed
% and cache files for single data generated
if ~(exist(config.divesDataCachePath, 'file') == 2)
    config.saveWithDateName = 1;
    % before starting, you may want to increase the amount of memory that
    % MATLAB can use.  Select Home - Preferences - General - Java Heap Memory
    % and use the scale bar to increase memory.
    %
    % If output file alteady exist AZFP processing will not happen - data
    % will be loaded from existing output
    Output = getAZFPProcessResult(config);

    % Make a test figure to make sure the data look right
    % getTestFigure(Output)

    % generating dives data
    [Dives, bottomDepth] = aggregateDivesData(Output, config);
% if current date is already processed then dive data is loaded and merged
% from all available dive files in "output" dir
else
    config.saveWithDateName = 0;
    % merging all output dives files
    filenames = getDivesFilenames(config.outputFolder);
    [Dives, bottomDepth] = mergeDives(filenames);
end
%% make matrices of each Sv variable
chartData = getChartsData(Dives);
%% Draw and save figures
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

depths = 1:config.maxDepth;
P_indices = [2, 3, 4];
freqLabels = {'200 kHz', '455 kHz', '769 kHz'};

% figure('Position', [100, 100, 1200, 900])
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

figure('Position', [100, 100, 1200, 900])

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
    clim = [ min(allVals,[],'omitnan'), min(max(allVals,[],'omitnan'), 1e4) ];


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

figure('Position', [100, 100, 1200, 900]);

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

