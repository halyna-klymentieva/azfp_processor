%% Use this code if you need to merge multiple 'Dive' matrices from multiple files

%filename = strcat("C:/Users/Andrea/Documents/AMesquita2025/UNB/dataAnalysis/preyData/processedAZFP/updatedCode/2024/","2024_day2_Dive.mat");
%save(filename, 'Dive','-v7.3');
%clear filename;

clc
clear variables
close all

cd('C:/Users/Andrea/Documents/AMesquita2025/UNB/dataAnalysis/preyData/processedAZFP/updatedCode/2023/entireMission')

% merge multiple structures
load('2023_11_Dive.mat')
Dive1=Dive;
bdepth1=bottomDepth
clear Dive Output
% load('2023_2_Dive.mat')
% Dive2=Dive;
% bdepth2=bottomDepth
% clear Dive Output
% load('2023_3_Dive.mat')
% Dive3=Dive;
% bdepth3=bottomDepth
% clear Dive Output
% load('2023_4_Dive.mat')
% Dive4=Dive;
% bdepth4=bottomDepth
% clear Dive Output
% load('2023_5_Dive.mat')
% Dive5=Dive;
% bdepth5=bottomDepth
% clear Dive Output
% load('2023_6_Dive.mat')
% Dive6=Dive;
% bdepth6=bottomDepth
% clear Dive Output
% load('2023_7_Dive.mat')
% Dive7=Dive;
% bdepth7=bottomDepth
% clear Dive Output
% load('2023_8_Dive.mat')
% Dive8=Dive;
% bdepth8=bottomDepth
% clear Dive Output
% load('2023_9_Dive.mat')
% Dive9=Dive;
% bdepth9=bottomDepth
% clear Dive Output
% load('2023_10_Dive.mat')
% Dive10=Dive;
% bdepth10=bottomDepth
% clear Dive Output
% load('2023_11_Dive.mat')
% Dive11=Dive;
% bdepth11=bottomDepth
% clear Dive Output
% load('2023_12_Dive.mat')
% Dive12=Dive;
% bdepth12=bottomDepth
% clear Dive Output
% load('2023_13_Dive.mat')
% Dive13=Dive;
% bdepth13=bottomDepth
% clear Dive Output
% load('2023_14_Dive.mat')
% Dive14=Dive;
% bdepth14=bottomDepth
% clear Dive Output
% load('2023_15_Dive.mat')
% Dive15=Dive;
% bdepth15=bottomDepth
% clear Dive Output

%% Pad short structures with NaNs

targetCols = 100;

for d = 1:1 %number of dive structures
    diveVarName = ['Dive' num2str(d)];
    dive = eval(diveVarName);
    
    for i = 1:numel(dive)
        for j = 1:numel(dive(i).P)
            fields = fieldnames(dive(i).P(j));
            for f = 1:numel(fields)
                val = dive(i).P(j).(fields{f});
                if isnumeric(val) && ismatrix(val)
                    [rows, cols] = size(val);
                    if cols < targetCols
                        dive(i).P(j).(fields{f}) = [val NaN(rows, targetCols - cols)];
                    end
                end
            end
        end
    end
    
    assignin('base', diveVarName, dive);
end

%% Remove inconsistent variable
for i = 1:numel(Dive1)
    if isfield(Dive1(i).P, 'N')
        Dive1(i).P = rmfield(Dive1(i).P, 'N');
    end
    if isfield(Dive1(i).P, 'mN')
        Dive1(i).P = rmfield(Dive1(i).P, 'mN');
    end
end

for i = 1:numel(Dive6)
    if isfield(Dive6, 'localtime')
        Dive6 = rmfield(Dive6, 'localtime');
    end
    if isfield(Dive6, 'tod')
        Dive6 = rmfield(Dive6, 'tod');
    end
end

%% Convert to cell array of structs
for d = 1:1
    diveVarName = ['Dive' num2str(d)];
    dive = eval(diveVarName);
    for i = 1:numel(dive)
        dive(i).P = {dive(i).P};
    end
    assignin('base', diveVarName, dive);
end

%% Convert structures to tables
aa_t  = struct2table(Dive1);
% bb_t  = struct2table(Dive2);
% cc_t  = struct2table(Dive3);
% dd_t  = struct2table(Dive4);
% ee_t  = struct2table(Dive5);
% ff_t  = struct2table(Dive6);
% gg_t  = struct2table(Dive7);
% hh_t  = struct2table(Dive8);
% ii_t  = struct2table(Dive9);
% jj_t  = struct2table(Dive10);
% kk_t  = struct2table(Dive11);
% ll_t  = struct2table(Dive12);
% mm_t  = struct2table(Dive13);
% nn_t  = struct2table(Dive14);
% oo_t  = struct2table(Dive15);
% Concatenate tables
merge_t = [aa_t;
%     bb_t; cc_t; dd_t; ee_t; ff_t; gg_t; hh_t; ii_t; jj_t; 
%     kk_t; ll_t; mm_t; nn_t; oo_t
    ];
% Convert table to structure
Dive = table2struct( merge_t )
% Merge bottom depth
bottomDepth = [bdepth1; 
%     bdepth2; bdepth3; bdepth4; bdepth5; bdepth6; bdepth7; 
%     bdepth8; bdepth9; bdepth10; bdepth11; bdepth12; bdepth13; bdepth14; bdepth15
    ];
%% make matrices of each Sv variable

nDives = length(Dive);
cols = 100;

Sv_130 = NaN(nDives, cols);
Sv_200 = NaN(nDives, cols);
Sv_455 = NaN(nDives, cols);
Sv_769 = NaN(nDives, cols);

Sv_200_130 = NaN(nDives, cols);
Sv_455_200 = NaN(nDives, cols);
Sv_769_400 = NaN(nDives, cols);

masked_130 = NaN(nDives, cols);
masked_200 = NaN(nDives, cols);
masked_455 = NaN(nDives, cols);
masked_769 = NaN(nDives, cols);

StartDiveTime = NaN(nDives, 1);
EndDiveTime = NaN(nDives, 1);

for i = 1:nDives
    %Pstruct = Dive(i).P{1,1};
    
    Sv_130(i,:)      = real(10*log10(Dive(i).P(1).msv));
    Sv_200(i,:)      = real(10*log10(Dive(i).P(2).msv));
    Sv_455(i,:)      = real(10*log10(Dive(i).P(3).msv));
    Sv_769(i,:)      = real(10*log10(Dive(i).P(4).msv));
    
    Sv_200_130(i,:)  = Dive(i).P(1).mDiff;
    Sv_455_200(i,:)  = Dive(i).P(2).mDiff;
    Sv_769_400(i,:)  = Dive(i).P(3).mDiff;
    
    masked_130(i,:)  = real(10*log10(Dive(i).P(1).mMasked));
    masked_200(i,:)  = real(10*log10(Dive(i).P(2).mMasked));
    masked_455(i,:)  = real(10*log10(Dive(i).P(3).mMasked));
    masked_769(i,:)  = real(10*log10(Dive(i).P(4).mMasked));
    
    StartDiveTime(i,:) = Dive(i).starttime;
    EndDiveTime(i,:)   = Dive(i).endtime;
end

%% Plot median Sv for all frequencies
figure(1)
clf

freqLabels = {'130 kHz', '200 kHz', '455 kHz', '769 kHz'};
Sv_all = {Sv_130, Sv_200, Sv_455, Sv_769};
caxisVals = [-110 -60];

for i = 1:4
    subplot(2,2,i)
    imagesc(Sv_all{i}', 'AlphaData', ~isnan(Sv_all{i}'))
    colormap('jet')
    caxis(caxisVals)
    ylim([0 120])
    xlabel('Ping Number')
    ylabel('Depth (m)')
    title(freqLabels{i})
    
    hold on
    plot(1:length(bottomDepth), bottomDepth, 'k-', 'LineWidth', 1, 'Color', [0 0 0 0.3])
    hold off
end

% Shared colorbar
h = colorbar;
ylabel(h, 'Volume Backscatter (dB re 1 m^{-1})')
h.Position(4) = 0.65;
h.Position(1) = .94 - h.Position(3);
h.Position(2) = 0.5 - h.Position(4)/2;

% save the figure
filename = strcat("C:/Users/Andrea/Documents/AMesquita2025/UNB/dataAnalysis/preyData/processedAZFP/updatedCode/2023/entireMission/medianSv2023.png");
print(gcf,'-dpng',filename,'-r0')

clear filename;
%close

%% Plot dB differences for all frequencies
figure(2)
clf

diffLabels = {'200-130 kHz', '455-200 kHz', '769-455 kHz'};
dbDiff_all = {Sv_200_130, Sv_455_200, Sv_769_400};
caxisVals = [0 11]

for i = 1:3
    subplot(2,2,i)
    imagesc(dbDiff_all{i}', 'AlphaData', ~isnan(dbDiff_all{i}'))
    colormap('jet')
    caxis(caxisVals)
    xlabel('Ping Number')
    ylabel('Depth (m)')
    title(diffLabels{i})
end

% Shared colorbar
h = colorbar;
ylabel(h, 'dB Difference')
h.Position(4) = 0.65;
h.Position(1) = .94 - h.Position(3);
h.Position(2) = 0.5 - h.Position(4)/2;

% save the figure
filename = strcat("C:/Users/Andrea/Documents/AMesquita2025/UNB/dataAnalysis/preyData/processedAZFP/updatedCode/2023/entireMission/dbDiff2023.png");
print(gcf,'-dpng',filename,'-r0')

clear filename;
%close

%% Plot masked Sv for all frequencies
figure(3)
clf

freqLabels = {'130 kHz', '200 kHz', '455 kHz', '769 kHz'};
Sv_all = {masked_130, masked_200, masked_455, masked_769};
caxisVals = [-110 -60];

for i = 1:4
    subplot(2,2,i)
    imagesc(Sv_all{i}', 'AlphaData', ~isnan(Sv_all{i}'))
    colormap('jet')
    caxis(caxisVals)
    ylim([0 120])
    xlabel('Ping Number')
    ylabel('Depth (m)')
    title(freqLabels{i})
end

% Shared colorbar
h = colorbar;
ylabel(h, 'Volume Backscatter (dB re 1 m^{-1})')
h.Position(4) = 0.65;
h.Position(1) = .94 - h.Position(3);
h.Position(2) = 0.5 - h.Position(4)/2;

% save the figure
filename = strcat("C:/Users/Andrea/Documents/AMesquita2025/UNB/dataAnalysis/preyData/processedAZFP/updatedCode/2023/entireMission/mMasked2023.png");
print(gcf,'-dpng',filename,'-r0')

clear filename;
%close

%% save variables

filename = strcat("C:/Users/Andrea/Documents/AMesquita2025/UNB/dataAnalysis/preyData/processedAZFP/updatedCode/2023/entireMission/","2023_Dive.mat");
save(filename,'Dive','Sv_130','Sv_200','Sv_455','Sv_769','Sv_200_130','Sv_455_200','Sv_769_400','masked_130','masked_200','masked_455','masked_769','StartDiveTime','EndDiveTime','bottomDepth','-v7.3');

%% calculate numerical density (ind. per m^3) from merged Dive structures

sigma_bs_jun = [2.1e-11, 9.9e-11, 8.8e-11];  % For June
sigma_bs_jul = [1.3e-11, 7.6e-11, 7.9e-11];  % For July
sigma_bs_aug  = [1.2e-11, 7.3e-11, 7.8e-11];  % For August
sigma_bs_sep = [1.7e-11, 8.7e-11, 8.1e-11];  % For September

n = length(Dive);

for ii = 1:n
    dt = datetime(Dive(ii).starttime, 'ConvertFrom', 'datenum');
    Dive(ii).localtime = dt - hours(3);

    mo = month(Dive(ii).localtime);
    if mo == 5
        sigma_bs = sigma_bs_jun;
    elseif mo == 6
        sigma_bs = sigma_bs_jun;    
    elseif mo == 7
        sigma_bs = sigma_bs_jul;
    elseif mo == 8
        sigma_bs = sigma_bs_aug;
    elseif mo == 9
        sigma_bs = sigma_bs_sep;
    end

    % Compute numerical density
    Dive(ii).P(2).N  = Dive(ii).P(2).masked / sigma_bs(1);  % 200 kHz
    Dive(ii).P(3).N  = Dive(ii).P(3).masked / sigma_bs(2);  % 455 kHz
    Dive(ii).P(4).N  = Dive(ii).P(4).masked / sigma_bs(3);  % 769 kHz

    % Median numerical densities
    Dive(ii).P(2).mN = nanmedian(Dive(ii).P(2).N);
    Dive(ii).P(3).mN = nanmedian(Dive(ii).P(3).N);
    Dive(ii).P(4).mN = nanmedian(Dive(ii).P(4).N);
end

%% scale density values by volume to estimate abundance


%% estimate biomass by incorporating mean organism weight


%% estimate energy density by incorporating kJ per gram of lipid


%% Plot numeric abundances for all frequencies
%Approximate civil twilight in July 2024: 7AM-8PM (day), 11PM-4AM (night)

%Numeric abundances were calculated from linearized masked Sv values divided
%by corresponding sigma_bs and averaged by dive profile
%Plot day and night profiles by depth

for i = 1:length(Dive)
    hr = hour(Dive(i).localtime);

    if hr >= 6 && hr <= 17
        Dive(i).tod = "Day";
    elseif hr >= 18 || hr <= 5
        Dive(i).tod = "Night";
    else
        Dive(i).tod = "Other";
    end
end

depths = 1:100;
freqLabels = {'200 kHz', '455 kHz', '769 kHz'};
P_indices = [2, 3, 4];

figure('Position', [100 100 1200 900])

for f = 1:3
    j = P_indices(f);

    dayProfiles = [];
    nightProfiles = [];

    for i = 1:length(Dive)
            profile = Dive(i).P(j).mN(:);
            if strcmpi(Dive(i).tod, 'Day')
                dayProfiles = [dayProfiles, profile];
            elseif strcmpi(Dive(i).tod, 'Night')
                nightProfiles = [nightProfiles, profile];
            end
    end

    allVals = [dayProfiles(:); nightProfiles(:)];
    clim = [nanmin(allVals), nanmax(allVals)];

    % Day
    subplot(3, 2, (f-1)*2 + 1)
    if ~isempty(dayProfiles)
        imagesc(dayProfiles)
        set(gca, 'YDir', 'reverse')
        yticks(10:10:100)
        yticklabels(string(10:10:100))
        xlabel('Dive Profile')
        ylabel('Depth (m)')
        title(['Daytime - ' freqLabels{f}])
        caxis(clim)
        colormap(hot)
        colorbar
    else
        text(0.5, 0.5, 'No day data', 'HorizontalAlignment', 'center')
        axis off
    end

    % Night
    subplot(3, 2, (f-1)*2 + 2)
    if ~isempty(nightProfiles)
        imagesc(nightProfiles)
        set(gca, 'YDir', 'reverse')
        yticks(10:10:100)
        yticklabels(string(10:10:100))
        xlabel('Dive Profile')
        ylabel('Depth (m)')
        title(['Nighttime - ' freqLabels{f}])
        caxis(clim)
        colormap(hot)
        colorbar
    else
        text(0.5, 0.5, 'No night data', 'HorizontalAlignment', 'center')
        axis off
    end
end

figure('Position', [100 100 1200 900])

for f = 1:3
    j = P_indices(f);

    dayProfiles = [];
    nightProfiles = [];

    for i = 1:length(Dive)
        profile = Dive(i).P(j).mN(:);
        if strcmpi(Dive(i).tod, 'Day')
            dayProfiles = [dayProfiles, profile];
        elseif strcmpi(Dive(i).tod, 'Night')
            nightProfiles = [nightProfiles, profile];
        end
    end

    allVals = [dayProfiles(:); nightProfiles(:)];
    clim = [nanmin(allVals), min(nanmax(allVals), 1e4)];

    % Day
    subplot(3, 2, (f-1)*2 + 1)
    if ~isempty(dayProfiles)
        imagesc(dayProfiles)
        set(gca, 'YDir', 'reverse')
        yticks(10:10:100)
        yticklabels(string(10:10:100))
        xlabel('Dive Profile')
        ylabel('Depth (m)')
        title(['Daytime - ' freqLabels{f}])
        caxis(clim)
        colormap(hot)
        colorbar
    else
        text(0.5, 0.5, 'No day data', 'HorizontalAlignment', 'center')
        axis off
    end

    % Night
    subplot(3, 2, (f-1)*2 + 2)
    if ~isempty(nightProfiles)
        imagesc(nightProfiles)
        set(gca, 'YDir', 'reverse')
        yticks(10:10:90)
        yticklabels(string(10:10:100))
        xlabel('Dive Profile')
        ylabel('Depth (m)')
        title(['Nighttime - ' freqLabels{f}])
        caxis(clim)
        colormap(hot)
        colorbar
    else
        text(0.5, 0.5, 'No night data', 'HorizontalAlignment', 'center')
        axis off
    end
end

figure('Position', [100 100 1200 900]);

for f = 1:3
    j = P_indices(f);
    allProfiles = [];
    todLabels = [];
    diveIDs = [];

    for i = 1:length(Dive)
        profile = Dive(i).P(j).mN(:);
        allProfiles = [allProfiles, profile];
        todLabels = [todLabels, string(Dive(i).tod)];
        diveIDs = [diveIDs, i];
    end

    clim = [0, 10000];
    normVals = (allProfiles - clim(1)) / (clim(2) - clim(1));
    normVals = min(max(normVals, 0), 1);

    hotMap = hot(256);
    boneMap = bone(256);

    idxVals = round(normVals * 255) + 1;
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

    ax = subplot(3,1,f);
    image(RGB)
    set(gca, 'YDir', 'reverse')
    yticks(10:10:90)
    yticklabels(string(10:10:100))
    ylabel('Depth (m)')
    xlabel('Dive Profile')
    title(freqLabels{f})

    tickStep = max(1, floor(cols / 10));
    xtickIdx = 1:tickStep:cols;
    xtickLabels = diveIDs(xtickIdx);
    xticks(xtickIdx)
    xticklabels(xtickLabels)

    axPos = get(ax, 'Position');
    cbHot = colorbar('Position', [axPos(1) - 0.05, axPos(2), 0.01, axPos(4)]);
    colormap(cbHot, hotMap);
    cbHot.Ticks = linspace(0, 1, 5);
    cbHot.TickLabels = arrayfun(@(x) sprintf('%.0f', clim(1) + x * (clim(2) - clim(1))), cbHot.Ticks, 'UniformOutput', false);
    cbHot.Label.String = 'Daytime';

    cbBone = colorbar('Position', [axPos(1) + axPos(3) + 0.01, axPos(2), 0.01, axPos(4)]);
    colormap(cbBone, boneMap);
    cbBone.Ticks = linspace(0, 1, 5);
    cbBone.TickLabels = arrayfun(@(x) sprintf('%.0f', clim(1) + x * (clim(2) - clim(1))), cbBone.Ticks, 'UniformOutput', false);
    cbBone.Label.String = 'Nighttime';
end

figure;
tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact')

for f = 1:3
    j = P_indices(f);
    dayMeans = [];
    nightMeans = [];
    
    for i = 1:length(Dive)
            profile = Dive(i).P(j).mN(:);
            meanVal = mean(profile, 'omitnan');
            if strcmpi(Dive(i).tod, 'Day')
                dayMeans(end+1) = meanVal;
            elseif strcmpi(Dive(i).tod, 'Night')
                nightMeans(end+1) = meanVal;
            end
    end
    
    nexttile
    hold on
    
    boxplot([dayMeans, nightMeans], ...
        [repmat({'Day'}, 1, length(dayMeans)), repmat({'Night'}, 1, length(nightMeans))], ...
        'Colors', 'k', 'Symbol', '')
    
    xDay = ones(size(dayMeans)) + 0.1*randn(size(dayMeans));
    xNight = 2*ones(size(nightMeans)) + 0.1*randn(size(nightMeans));
    scatter(xDay, dayMeans, 30, [0.6, 0.8, 1], 'filled', 'MarkerFaceAlpha', 0.6)
    scatter(xNight, nightMeans, 30, [0, 0, 0.5], 'filled', 'MarkerFaceAlpha', 0.6)
    
    ylim([0 50000])
    
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

%% make matrices of numerical densities

for i = 1:nDives
    %Pstruct = Dive(i).P{1,1};
    
    mN_200(i,:)  = Dive(i).P(2).mN;
    mN_455(i,:)  = Dive(i).P(3).mN;
    mN_769(i,:)  = Dive(i).P(4).mN;    
end

%% save variables

filename = strcat("C:/Users/Andrea/Documents/AMesquita2025/UNB/dataAnalysis/preyData/processedAZFP/updatedCode/2023/entireMission/","2023_Dive_densities.mat");
save(filename,'Dive','Sv_130','Sv_200','Sv_455','Sv_769','Sv_200_130','Sv_455_200','Sv_769_400','masked_130','masked_200','masked_455','masked_769','mN_200','mN_455','mN_769','StartDiveTime','EndDiveTime','bottomDepth','-v7.3');
