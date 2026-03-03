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
%% User-defined variables
dateOfData = '25-07-24';
xmlFileName = '25072317.XML';
gliderVariableName = 'cabot_20250723_213_delayed';
gliderFileName = 'cabot_20250723_213_delayed_0660_8583_ae5f.mat';
calibrationOffsets = [-0.72, -4.60, -4.63, -0.20];
% Far field per frequency cut-off range, see filterFarField function
farFieldCutOffRange = [30, 10, 10, 5];
dataCutOffDepth = 67.5;
surfaceNoiseDepthMin = 10; % m
surfaceNoiseDepthMax = 14; % m
maxDepth = 105;
%% Global variables
azfpDataCacheFilename = "data-" + dateOfData + '.mat';
divesDataCacheFilename = "dives-" + dateOfData + '.mat';
sourceFolder = fullfile(pwd, '..', 'data');
sourceFileNames = getFilenamesByDate(dateOfData, sourceFolder);
fprintf('For data %s files count: %d.\n', dateOfData, length(sourceFileNames));
azfpDataCachePath = fullfile(pwd, '..', 'output', azfpDataCacheFilename);
divesDataCachePath = fullfile(pwd, '..', 'output', divesDataCacheFilename);
gliderFileFullName = fullfile(pwd, '..', 'gliderData', gliderFileName);
%% Getting output of AZFP raw data processing
% before starting, you may want to increase the amount of memory that
% MATLAB can use.  Select Home - Preferences - General - Java Heap Memory
% and use the scale bar to increase memory.
%
% If output file alteady exist AZFP processing will not happen - data
% will be loaded from existing output
Output = getAZFPProcessResult(azfpDataCachePath, sourceFolder, sourceFileNames, xmlFileName);
%% Make a test figure to make sure the data look right
% getTestFigure(Output)
%% Standard Sphere Calibration application - Halyna used 2025 callibration data: 130 kHz (-0.72), 200 kHZ (-4.60), 455 kHZ (-4.63),769 kHZ (-0.20),

% Echosounder 59016 (Davies) calibration offset from standard sphere
% calibration conducted June 2024.  Davies and Mesquita have
% calibration files.
Output = standardSphereCallibration(calibrationOffsets, Output);
%% Trim Transmit Pulse and Near Field from Sv data (Step 5 in tutorial)
Output = filterNearField(Output);
%% Load in glider data
[gliderDepth, gliderTime, gliderData] = loadGliderData(gliderFileFullName, gliderVariableName);
%% glider data filtering for extra pings during inactivity
% Output = gliderDataFilter(gliderData, Output);
%% Time - align the glider and AZFP pressure data (Step 2 in tutorial)
[Output, fixedDepth] = timeAlignAZFPToGliderData(Output, gliderTime, gliderDepth);
fixedDepth(end) = [];

%% Far field noise cut off  (Steps 3 and 4 in the tutorial)
Output = filterFarField(Output, fixedDepth, farFieldCutOffRange);

%% Remove pings at the surface of the ocean because these often have bubbles in them
[Output, Dives] = filterSurfacePings(Output);
%% Histogram of Sv data for all frequencies
% getFigureAllFreqSvHystogram(Output)
%% Make a histogram of the seafloor data decibel strengths
% getFigureSeafloorDecibelStrength(Output)
%% Remove all data below a defined depth
Output = filterDepthDataCutOff(Output, Dives, dataCutOffDepth);

%% Remove seafloor echoes (Step 1 in the tutorial)
Output = filterSeafloorEchoes(Output, Dives);

%% Remove surface noise
Output = filterSurfaceNoise(Output, Dives, surfaceNoiseDepthMin, surfaceNoiseDepthMax);

%% Average 10 cm vertical resolution into 1 m depth bins to make the matrices smaller for better storage space
% Step 6: Organized Dive Structure (Stable Version)- Halyna's version
clear Dives; % Start with a blank slate
[azfpAggregateData, Dives] = aggregateVerticalResolution(Output, maxDepth);

%% Moving average to determine noise floor for each bin (Steps 3 and 4) and noise removal
Dives = filterNoiseFloor(maxDepth, Output, azfpAggregateData, Dives);
%% AZFP_Unmasked_Masked_Comparison routine (D. Mossman) - lines 626–750
Dives = unmaskedMaskedComparison(Dives);
%% Use this code to save Dives data to file
saveDiveData(divesDataCachePath, Dives)
%% Use this code if you need to merge multiple 'Dive' matrices from multiple files

% cd('C:/Users/Andrea/Documents/AMesquita2025/UNB/dataAnalysis/preyData/processedAZFP/updatedCode/2024/entireMission')
% filenames = {
%     "2024_1_Dive.mat", ...
%     "2024_2_Dive.mat", ...
%     "2024_3_Dive.mat", ...
%     "2024_4_Dive.mat"
% };
% nFiles = numel(filenames);
% DivesData = cell(1, nFiles);   % Dives{k} contains the Dive from file k% merge multiple structures
% for k = 1:nFiles
%     S = load(filenames{k}, 'Dives');
%     if isfield(S, 'Dives')
%         DivesData{k} = S.Dives;
%     else
%         DivesData{k} = [];  % or handle missing variable
%     end
% end

% % Variable: Dive
% % Convert structures to tables
% aa_t = struct2table( DivesData{1} );
% bb_t = struct2table( DivesData{2} );
% cc_t = struct2table( DivesData{3} );
% dd_t = struct2table( DivesData{4} );
% % Concatonate tables
% merge_t = [ aa_t ; bb_t ; cc_t ; dd_t ];
% % Convert table to structure
% Dive = table2struct( merge_t );
%% make matrices of each Sv variable

for i = 1:length(Dives)
    Sv_130(i, :) = real(10*log10(Dives(i).P(1).msv));
    Sv_200(i, :) = real(10*log10(Dives(i).P(2).msv));
    Sv_455(i, :) = real(10*log10(Dives(i).P(3).msv));
    Sv_769(i, :) = real(10*log10(Dives(i).P(4).msv));
    Sv_200_130(i, :) = Dives(i).P(1).mDiff;
    Sv_455_200(i, :) = Dives(i).P(2).mDiff;
    Sv_769_400(i, :) = Dives(i).P(3).mDiff;
    masked_130(i, :) = real(10*log10(Dives(i).P(1).mMasked));
    masked_200(i, :) = real(10*log10(Dives(i).P(2).mMasked));
    masked_455(i, :) = real(10*log10(Dives(i).P(3).mMasked));
    masked_769(i, :) = real(10*log10(Dives(i).P(4).mMasked));
    StartDiveTime(i, :) = Dives(i).starttime;
    EndDiveTime(i, :) = Dives(i).endtime;
end
Depth = 1:80;
%% mask echoes below the seafloor using 769 kHz matrix (or masked 455 kHz) - ask Andrea why she disabled lines 883-910

idx = NaN(size(masked_455, 1), 1);

for i = 1:size(masked_455, 1)
    rowData = masked_455(i, :);

    for j = size(rowData, 2):-1:1
        if ~isnan(rowData(j))
            idx(i) = j;
            break;
        end
    end
end

% Mask values below the seafloor (after the last valid bin)
for i = 1:size(Sv_130, 1)
    colStart = idx(i) + 1;
    if ~isnan(colStart) && colStart <= size(Sv_130, 2)
        Sv_130(i, colStart:end) = NaN;
        Sv_200(i, colStart:end) = NaN;
        Sv_455(i, colStart:end) = NaN;
        Sv_200_130(i, colStart:end) = NaN;
        Sv_455_200(i, colStart:end) = NaN;
        masked_130(i, colStart:end) = NaN;
        masked_200(i, colStart:end) = NaN;
        masked_455(i, colStart:end) = NaN;

    end
end

% obtain bottom depth by dive
nDives = length(Dives);
endDive = zeros(nDives, 1);
bottomDepth = zeros(nDives, 1);

for i = 1:nDives
    endDive(i) = Dives(i).Index(2) - 1;
    bottomDepth(i) = Output(1).Depth(endDive(i));
end
%% Plot median Sv for all frequencies; Save a plot
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
    ylim([0, maxDepth])
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
figure1MedianSVFileName = "figure1-median-sv-" + dateOfData + '.png';
figure1MedianSVPath = fullfile(pwd, '..', 'output', figure1MedianSVFileName);
print(gcf, '-dpng', figure1MedianSVPath, '-r0')

clear filename;
%close
%% Plot dB differences for all frequencies; Save a plot
figure(2)
clf

diffLabels = {'200-130 kHz', '455-200 kHz', '769-455 kHz'};
dbDiff_all = {Sv_200_130, Sv_455_200, Sv_769_400};
caxisVals = [3, 9];

for i = 1:3
    subplot(2, 2, i)
    imagesc(dbDiff_all{i}', 'AlphaData', ~isnan(dbDiff_all{i}'))
    colormap('jet')
    caxis(caxisVals)
    xlabel('Ping Number')
    ylabel('Depth (m)')
    title(diffLabels{i})

    %     hold on
    %     plot(1:length(bottomDepth), bottomDepth, 'k-', 'LineWidth', 2)
    %     hold off
end

% Shared colorbar
h = colorbar;
ylabel(h, 'dB Difference')
h.Position(4) = 0.65;
h.Position(1) = .94 - h.Position(3);
h.Position(2) = 0.5 - h.Position(4) / 2;

% save the figure
figure2FileName = "figure2-db-diff-" + dateOfData + '.png';
figure2Path = fullfile(pwd, '..', 'output', figure2FileName);
print(gcf, '-dpng', figure2Path, '-r0')

clear filename;
%close
%% Plot masked Sv for all frequencies; save a plot
figure(3)
clf

freqLabels = {'130 kHz', '200 kHz', '455 kHz', '769 kHz'};
Sv_all = {masked_130, masked_200, masked_455, masked_769};
%caxisVals = [-110 -60];

for i = 1:4
    subplot(2, 2, i)
    imagesc(Sv_all{i}', 'AlphaData', ~isnan(Sv_all{i}'))
    colormap('jet')
    %caxis(caxisVals)
    ylim([0, maxDepth])
    xlabel('Ping Number')
    ylabel('Depth (m)')
    title(freqLabels{i})

    %     hold on
    %     plot(1:length(bottomDepth), bottomDepth, 'k-', 'LineWidth', 2)
    %     hold off
end

% Shared colorbar
h = colorbar;
ylabel(h, 'Volume Backscatter (dB re 1 m^{-1})')
h.Position(4) = 0.65;
h.Position(1) = .94 - h.Position(3);
h.Position(2) = 0.5 - h.Position(4) / 2;

% save the figure
figure3FileName = "figure3-m-masked-" + dateOfData + '.png';
figure3Path = fullfile(pwd, '..', 'output', figure3FileName);
print(gcf, '-dpng', figure3Path, '-r0')
%close
%% save variables

% filename = strcat("d:/AZFP_GSL_2025/AZFP_processed_data/", "July25_Dive.mat");
% save(filename, 'Dive', 'Sv_130', 'Sv_200', 'Sv_455', 'Sv_769', 'Sv_200_130', 'Sv_455_200', 'Sv_769_400', 'masked_130', 'masked_200', 'masked_455', 'masked_769', 'StartDiveTime', 'EndDiveTime', 'bottomDepth', '-v7.3');
%% calculate numerical density (ind. per m^3) from merged Dive structures

sigma_bs = [1.3e-11, 7.6e-11, 7.9e-11]; % for 200kHz, 455kHz, 769kHz

n = length(Dives);

for ii = 1:n
    Dives(ii).P(2).N = Dives(ii).P(2).masked / sigma_bs(1); % 200 kHz
    Dives(ii).P(3).N = Dives(ii).P(3).masked / sigma_bs(2); % 455 kHz
    Dives(ii).P(4).N = Dives(ii).P(4).masked / sigma_bs(3); % 769 kHz

    Dives(ii).P(2).mN = nanmedian(Dives(ii).P(2).N);
    Dives(ii).P(3).mN = nanmedian(Dives(ii).P(3).N);
    Dives(ii).P(4).mN = nanmedian(Dives(ii).P(4).N);
end
% scale density values by volume to estimate abundance
% estimate biomass by incorporating mean organism weight
% estimate energy density by incorporating kJ per gram of lipid

%% Plot numeric abundances for all frequencies
%Approximate civil twilight in July 2024: 7AM-8PM (day), 11PM-4AM (night)

%Numeric abundances were calculated from linearized masked Sv values divided
%by corresponding sigma_bs and averaged by dive profile
%Plot day and night profiles by depth

% for i = 1:length(Dives)
%     dt = datetime(Dives(i).starttime, 'ConvertFrom', 'datenum'); %convert time
%     Dives(i).localtime = dt - hours(3); %convert to local time
% end
% 
% for i = 1:length(Dives)
%     hr = hour(Dives(i).localtime);
% 
%     if hr >= 7 && hr <= 20
%         Dives(i).tod = "Day";
%     elseif hr >= 23 || hr <= 4
%         Dives(i).tod = "Night";
%     else
%         Dives(i).tod = "Other";
%     end
% end
% 
% depths = 1:maxDepth;
% freqLabels = {'200 kHz', '455 kHz', '769 kHz'};
% P_indices = [2, 3, 4];
% 
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
%         yticks(10:10:maxDepth)
%         yticklabels(string(10:10:maxDepth))
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
%         yticks(10:10:maxDepth)
%         yticklabels(string(10:10:maxDepth))
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
% 
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
%     clim = [ min(allVals,[],'omitnan'), min(max(allVals,[],'omitnan'), 1e4) ];
% 
% 
%     % Day
%     subplot(3, 2, (f - 1)*2+1)
%     if ~isempty(dayProfiles)
%         imagesc(dayProfiles)
%         set(gca, 'YDir', 'reverse')
%         yticks(10:10:maxDepth)
%         yticklabels(string(10:10:maxDepth))
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
%         yticks(10:10:maxDepth)
%         yticklabels(string(10:10:maxDepth))
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
% 
% figure('Position', [100, 100, 1200, 900]);
% 
% for f = 1:3
%     j = P_indices(f);
%     allProfiles = [];
%     todLabels = [];
%     diveIDs = [];
% 
%     for i = 1:length(Dives)
%         profile = Dives(i).P(j).mN(:);
%         allProfiles = [allProfiles, profile];
%         todLabels = [todLabels, string(Dives(i).tod)];
%         diveIDs = [diveIDs, i];
%     end
% 
%     clim = [0, 10000];
%     normVals = (allProfiles - clim(1)) / (clim(2) - clim(1));
%     normVals = min(max(normVals, 0), 1);
% 
%     hotMap = hot(256);
%     boneMap = bone(256);
% 
%     idxVals = round(normVals*255) + 1;
%     [rows, cols] = size(allProfiles);
%     RGB = ones(rows, cols, 3);
% 
%     for col = 1:cols
%         if strcmpi(todLabels(col), 'Day')
%             cmap = hotMap;
%         else
%             cmap = boneMap;
%         end
%         for ch = 1:3
%             RGB(:, col, ch) = cmap(idxVals(:, col), ch);
%         end
%     end
% 
%     ax = subplot(3, 1, f);
%     image(RGB)
%     set(gca, 'YDir', 'reverse')
%     yticks(10:10:maxDepth)
%     yticklabels(string(10:10:maxDepth))
%     ylabel('Depth (m)')
%     xlabel('Dive Profile')
%     title(freqLabels{f})
% 
%     tickStep = max(1, floor(cols/10));
%     xtickIdx = 1:tickStep:cols;
%     xtickLabels = diveIDs(xtickIdx);
%     xticks(xtickIdx)
%     xticklabels(xtickLabels)
% 
%     axPos = get(ax, 'Position');
%     cbHot = colorbar('Position', [axPos(1) - 0.05, axPos(2), 0.01, axPos(4)]);
%     colormap(cbHot, hotMap);
%     cbHot.Ticks = linspace(0, 1, 5);
%     cbHot.TickLabels = arrayfun(@(x) sprintf('%.0f', clim(1)+x*(clim(2) - clim(1))), cbHot.Ticks, 'UniformOutput', false);
%     cbHot.Label.String = 'Daytime';
% 
%     cbBone = colorbar('Position', [axPos(1) + axPos(3) + 0.01, axPos(2), 0.01, axPos(4)]);
%     colormap(cbBone, boneMap);
%     cbBone.Ticks = linspace(0, 1, 5);
%     cbBone.TickLabels = arrayfun(@(x) sprintf('%.0f', clim(1)+x*(clim(2) - clim(1))), cbBone.Ticks, 'UniformOutput', false);
%     cbBone.Label.String = 'Nighttime';
% end
% 
figure;
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact')

for f = 1:3
    j = P_indices(f);
    dayMeans = [];
    nightMeans = [];

    for i = 1:length(Dives)
        profile = Dives(i).P(j).mN(:);
        meanVal = mean(profile, 'omitnan');
        if strcmpi(Dives(i).tod, 'Day')
            dayMeans(end+1) = meanVal;
        elseif strcmpi(Dives(i).tod, 'Night')
            nightMeans(end+1) = meanVal;
        end
    end

    nexttile
    hold on

    rng default  
    x1 = normrnd(5,1,100,1);
    x2 = normrnd(6,1,100,1);
    figure
    boxplot([x1,x2],'Notch','on','Labels',{'mu = 5','mu = 6'})
    title('Compare Random Data from Different Distributions')

    x = [dayMeans, nightMeans];   % numeric column
    g = [repmat({'Day'}, numel(dayMeans), 1); repmat({'Night'}, numel(nightMeans), 1)];
    boxplot(x, g, 'Colors', 'k', 'Symbol', '');    

    % boxplot([dayMeans, nightMeans], ...
    %     [repmat({'Day'}, 1, length(dayMeans)), repmat({'Night'}, 1, length(nightMeans))], ...
    %     'Colors', 'k', 'Symbol', '')

    xDay = ones(size(dayMeans)) + 0.1 * randn(size(dayMeans));
    xNight = 2 * ones(size(nightMeans)) + 0.1 * randn(size(nightMeans));
    scatter(xDay, dayMeans, 30, [0.6, 0.8, 1], 'filled', 'MarkerFaceAlpha', 0.6)
    scatter(xNight, nightMeans, 30, [0, 0, 0.5], 'filled', 'MarkerFaceAlpha', 0.6)

    ylim([0, 10000])

    set(gca, 'XTickLabel', {'Day', 'Night'})
    ylabel('Mean mN value')
    title(freqLabels{f})
    hold off
end

nexttile
axis off
legend({'Day points', 'Night points'}, 'Location', 'best')

% save the figure
%filename = strcat("C:/Users/Andrea/Documents/AMesquita2025/UNB/dataAnalysis/preyData/processedAZFP/updatedCode/2024/mMasked2024_1_4.png");
%print(gcf,'-dpng',filename,'-r0')

%clear filename;
%close
%%  Test Plots
% Test plot data by profile
% figure(1)
%
% subplot(2,1,1)
% imagesc([NaN * ones(5, size(P(1).avg_sv, 1)); 10*log10(abs(P(1).avg_sv'))],'AlphaData',~isnan([NaN * ones(5, size(P(1).avg_sv, 1)); 10*log10(abs(P(1).avg_sv'))])) % conversion back to decibels + dealing with the fact that we cut off the first 5 m of data
% colormap('jet');
% caxis([-100 -70]);
% % set(gca, 'Xdir', 'reverse');
% xlabel('Time')
% ylabel('Depth (m)')
% title('130 kHz')
% xt = get(gca,'XTick');
% xtlbl = [];
% for i = 1:numel(xt)
%     temp = Output(1).Date(xt(i));
%     temp = datetime(temp, 'ConvertFrom','datenum','Format','HH:mm');
%     temp = char(temp);
%     xtlbl = [xtlbl;temp];
% end
%
% set(gca, 'XTick',xt, 'XTickLabel',xtlbl, 'XTickLabelRotation',30)
% hold on
% %line(1:length(bott_dep2),bott_dep2(1,1:length(bott_dep2)),'Color','r','LineWidth',1)
% hold off
% colorbar
%
% subplot(2,1,2)
% imagesc([NaN * ones(5, size(P(2).avg_sv, 1)); 10*log10(abs(P(2).avg_sv'))],'AlphaData',~isnan([NaN * ones(5, size(P(2).avg_sv, 1)); 10*log10(abs(P(2).avg_sv'))]))
% colormap('jet');
% caxis([-110 -40]);
% % set(gca, 'Xdir', 'reverse');
% xlabel('Ping Number')
% ylabel('Depth (m)')
% ylim([0,100])
% title('200 kHz')
% set(gca, 'XTick',xt, 'XTickLabel',xtlbl, 'XTickLabelRotation',30)
% hold on
% %line(1:length(bott_dep2),bott_dep2(2,1:length(bott_dep2)),'Color','r','LineWidth',1)
% hold off
% colorbar
%
% subplot(2,2,3)
% imagesc([NaN * ones(5, size(P(3).avg_sv, 1)); 10*log10(abs(P(3).avg_sv'))],'AlphaData',~isnan([NaN * ones(5, size(P(3).avg_sv, 1)); 10*log10(abs(P(3).avg_sv'))]))
% colormap('jet');
% caxis([-80 -40]);
% % set(gca, 'Xdir', 'reverse');
% xlabel('Ping Number')
% ylabel('Depth (m)')
% ylim([0,100])
% title('455 kHz')
% set(gca, 'XTick',xt, 'XTickLabel',xtlbl, 'XTickLabelRotation',30)
% hold on
% %line(1:length(bott_dep2),bott_dep2(3,1:length(bott_dep2)),'Color','r','LineWidth',1)
% hold off
%
% subplot(2,2,4)
% imagesc([NaN * ones(5, size(P(4).avg_sv, 1)); 10*log10(abs(P(4).avg_sv'))],'AlphaData',~isnan([NaN * ones(5, size(P(4).avg_sv, 1)); 10*log10(abs(P(4).avg_sv'))]))
% caxis([-80 -40]);
% % set(gca, 'Xdir', 'reverse');
% xlabel('Ping Number')
% ylabel('Depth (m)')
% ylim([0,100])
% title('769 kHz')
% set(gca, 'XTick',xt, 'XTickLabel',xtlbl, 'XTickLabelRotation',30)
% hold on
% %line(1:length(bott_dep2),bott_dep2(4,1:length(bott_dep2)),'Color','r','LineWidth',1)
% hold off
%
% % single colorbar for all four plots
% h = colorbar;
% set(get(h,'label'),'string','Sv (dB scattering per unit volume)');
%
% % h.Position(4) = 0.65;
% % h.Position(1) = .94-h.Position(3);
% % h.Position(2) = 0.5-h.Position(4)/2;
% %
% % AddLetters2Plots(gcf,'VShift',-0.04)
% % % save the figure
% % filename = strcat("/Users/dmossman/Box/2022 MSc Thesis Work/Visuals/MATLAB Echosounder Figures/",date,"Sept/Frequencies_All_",date,"Sept.png");
% % print(gcf,'-dpng',filename,'-r0')
% %
% % clear filename;
% % %close
%
% %Plot
% figure
% imagesc(Sv_455')
% p=colorbar
% ylabel(p,'Volume Backscatter at 769 kHz P2')
% caxis([-110 -60])
% ylim([0 100])
% colormap(jet)
%% Plot a subset of the data - Andréa added this section to QAQC AZFP data

% startTimes = datetime(StartDiveTime, 'ConvertFrom', 'datenum');
%
% % Define periods as [month startDay]
% periods = [
%     7, 21;
%     8, 14;
%     9,  9
% ];
%
% for p = 1:size(periods,1)
%     month = periods(p, 1);
%     startDay = periods(p, 2);
%
%     figure
%     for d = 0:2
%         thisDay = datetime(2018, month, startDay + d);
%         dayMask = isbetween(startTimes, thisDay, thisDay + days(1));
%         Sv_day = Sv_455(dayMask, :);
%
%         subplot(2,2,d+1)  % arrange as 2 rows x 2 columns
%         imagesc(Sv_day')
%         colorbar
%         ylabel('Depth (m)')
%         xlabel('Dive index (Ping #)')
%         title(['Sv at 455 kHz on ', datestr(thisDay, 'mmm dd')])
%         caxis([-110 -60])
%         ylim([0 100])
%         colormap(jet)
%     end
%
%     % Add a title for the whole figure
%     sgtitle(['Sv at 455 kHz — ', datestr(datetime(2018, month, startDay), 'mmmm'), ...
%              ' ', num2str(startDay), '–', num2str(startDay+2)])
% end
%% extra plots
% %test plot the Sv by dive for the 130 kHz
% divedata=[Dive(1).P(1).sv];
% for ii=2:length(Dive)
%     divedata=[divedata;Dive(ii).P(1).sv];
% end;
%
% figure(1)
% imagesc([NaN * ones(5, size(divedata, 1)); 10*log10(abs(divedata'))]) % conversion back to decibels + dealing with the fact that we cut off the first 5 m of data
% colormap('jet');
% caxis([-100 -40]);
% xlabel('Ping Number')
% ylabel('Depth (m)')
% title('130 kHz')
% colorbar
%
% % plot single dives
% figure
% imagesc(10*log10(abs(Dive(6).P(1).sv))')
% colormap('jet')
% caxis([-80 -60]);
% ylim([0 150])
%% Subset 2023 data

% rowIdx = 142:365;
%
% for n = 1:numel(Output)
%     fields = fieldnames(Output(n));
%
%     for k = 1:numel(fields)
%         thisField = Output(n).(fields{k});
%
%         if isnumeric(thisField) && size(thisField,1) >= max(rowIdx)
%             Output(n).(fields{k}) = thisField(rowIdx,:);
%         end
%     end
% end
%
% removeIdx = 144:364;
%
% for n = 1:numel(Output)
%     fields = fieldnames(Output(n));
%
%     for k = 1:numel(fields)
%         thisField = Output(n).(fields{k});
%
%         if isnumeric(thisField) || isdatetime(thisField) || iscell(thisField)
%             if size(thisField,1) >= max(removeIdx)
%                 keepIdx = setdiff(1:size(thisField,1), removeIdx);
%                 Output(n).(fields{k}) = thisField(keepIdx, :);
%             end
%         end
%     end
% end
