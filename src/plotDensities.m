%% Plot prey densities concurrently with glider flight
clc
clear variables
close all

addpath(genpath('C:/Users/Andrea/Documents/AMesquita2025/UNB/dataAnalysis/preyData/matlab/november2023Code/AZFP Code - 17 Nov 2023/Code'))

% load processed AZFP data
cd('C:/Users/Andrea/Documents/AMesquita2025/UNB/dataAnalysis/preyData/processedAZFP/updatedCode/2023/entireMission')

load('2023_Dive_densities.mat')
Dive=Dive;

%% Load in and subset glider data
load('/Users/Andrea/Documents/AMesquita2025/UNB/dataAnalysis/preyData/matlab/missions/shad_20240703_196_delayed_8a8c_3290_a513.mat');

% put that data into the variable gliderdata
gliderdata = shad_20240703_196_delayed;
clear shad_20240703_196_delayed

% change glider unix time format to same format as in Output file (matlab
% time format)
unix_epoch = datenum(1970,1,1,0,0,0);
gliderdata.time = gliderdata.time./86400 + unix_epoch;

% convert datetime and time zone
gliderdata.time = datetime(gliderdata.time, 'ConvertFrom', 'datenum');
gliderdata.localtime = gliderdata.time - hours(3);

% subset hourly data
gliderh = dateshift(gliderdata.time, 'start', 'hour');
[~, ia] = unique(gliderh);

fn = fieldnames(gliderdata);
for k = 1:numel(fn)
    if isvector(gliderdata.(fn{k})) && length(gliderdata.(fn{k})) == length(gliderdata.time)
        gliderdataH.(fn{k}) = gliderdata.(fn{k})(ia);
    else
        gliderdataH.(fn{k}) = gliderdata.(fn{k});
    end
end

clear gliderdata

% export glider lat-lon data for checking projection in ArcGIS
lat = gliderdataH.latitude(:);
lon = gliderdataH.longitude(:);
ltime = gliderdataH.localtime(:);
T = table(lat, lon, ltime);
writetable(T, 'gliderdataH.csv');

%% Load GEBCO bathymetry data
tif_path = 'C:\Users\Andrea\Documents\AMesquita2025\UNB\dataAnalysis\preyData\matlab\november2023Code\AZFP Code - 17 Nov 2023\Code\bathymetry1.tif';
[bathy, R] = readgeoraster(tif_path);

% Create latitude and longitude grids from spatial reference
[lat, lon] = meshgrid(R.LatitudeLimits(1):R.CellExtentInLatitude:R.LatitudeLimits(2)-R.CellExtentInLatitude, ...
                      R.LongitudeLimits(1):R.CellExtentInLongitude:R.LongitudeLimits(2)-R.CellExtentInLongitude);
lat = lat';  % Match raster orientation
lon = lon';

% Define bounding box from gliderdataH
lat_min = min(gliderdataH.latitude) -1.4;
lat_max = max(gliderdataH.latitude) -0.96;
lon_min = min(gliderdataH.longitude) -0.1;
lon_max = max(gliderdataH.longitude) +0.05;

% Create mask based on bounding box
lat_mask = lat >= lat_min & lat <= lat_max;
lon_mask = lon >= lon_min & lon <= lon_max;
mask = lat_mask & lon_mask;

% Subset bathymetry and coordinates
bathy_subset = bathy(mask);
lat_subset = lat(mask);
lon_subset = lon(mask);

% invert latitude to correct projection
lat_subset = flipud(lat_subset);

% read -9999 as NaNs
bathy_subset(bathy_subset == -9999) = NaN;

% figure;
% scatter(lon_subset, lat_subset, 10, -double(bathy_subset), 'filled');
% xlabel('Longitude');
% ylabel('Latitude');
% title('Shediac Valley Bathymetry');
% colormap(parula);
% colorbar;
% axis equal tight;

%% Project data
% mstruct = defaultm('tranmerc');   % Transverse Mercator
% mstruct.origin = [0 -63 0];       % Central meridian at -63°
% mstruct.falseeasting  = 304800;
% mstruct.falsenorthing = 0;
% mstruct.scalefactor   = 0.9999;
% mstruct.geoid = wgs84Ellipsoid('meters');
% mstruct = defaultm(mstruct);
% 
% [x_glider, y_glider] = mfwdtran(mstruct, gliderdataH.latitude, gliderdataH.longitude);
% [x_bathy, y_bathy] = mfwdtran(mstruct, lat_subset, lon_subset);
% 
% figure
% scatter(x_bathy(:), y_bathy(:), 10, -double(bathy_subset(:)), 'filled');
% hold on
% plot(x_glider, y_glider, 'k-', 'LineWidth', 1.2);
% scatter(x_glider, y_glider, 20, 'r', 'filled');
% axis equal
% xlabel('Easting (m)'); ylabel('Northing (m)');
% title('Glider Track over Bathymetry (EPSG:32105)');
% colormap(parula); colorbar;

%manually correct offset
offset = 48.1 - 46.8;
lat_subset = lat_subset + offset;

%% Create Animated Map of Glider Flight
figure('Position', [100 100 900 700], 'Color', 'w');

% Plot bathymetry
scatter(lon_subset(:), lat_subset(:), 10, -double(bathy_subset(:)), 'filled');
xlabel('Longitude'); ylabel('Latitude');
title('Glider Mission with Bathymetry');
colormap(parula); colorbar;
axis equal tight; hold on;

% Preplot full track for context
plot(gliderdataH.longitude, gliderdataH.latitude, 'k-', 'LineWidth', 0.5);

% Initialize animated path
h_progress = plot(gliderdataH.longitude(1), gliderdataH.latitude(1), ...
    'r-', 'LineWidth', 2);
h_marker = plot(gliderdataH.longitude(1), gliderdataH.latitude(1), ...
    'ro', 'MarkerFaceColor', 'y', 'MarkerSize', 6);

% Timestamp text
time_text = text(min(gliderdataH.longitude), max(gliderdataH.latitude)+0.2, ...
    ['Time: ' datestr(gliderdataH.localtime(1))], 'FontSize', 12, ...
    'BackgroundColor', 'w', 'EdgeColor', 'k');

% Animation loop
frame_delay = 0.3;
for k = 2:length(gliderdataH.longitude)
    h_progress.XData = gliderdataH.longitude(1:k);
    h_progress.YData = gliderdataH.latitude(1:k);
    h_marker.XData = gliderdataH.longitude(k);
    h_marker.YData = gliderdataH.latitude(k);
    time_text.String = ['Time: ' datestr(gliderdataH.localtime(k))];
    pause(frame_delay)
    drawnow
end

%% Plot prey data (sliding window)
for i = 1:length(Dive)
    hr = hour(Dive(i).localtime);

    if hr >= 7 && hr <= 20
        Dive(i).tod = "Day";
    elseif hr >= 23 || hr <= 4
        Dive(i).tod = "Night";
    else
        Dive(i).tod = "Other";
    end
end

figure('Position', [100 100 1200 900]);

% Extract all datetimes
diveTimes = [Dive.localtime];  % Array of datetime objects

% Verify matrix orientation
if size(mN_455, 1) ~= numel(Dive)
    % Transpose if needed (rows should be dives, columns depth)
    mN_455 = mN_455';
end

% Final size check
if size(mN_455, 1) ~= numel(Dive)
    error('Matrix orientation mismatch: mN_455 has %d rows, but there are %d dives', ...
          size(mN_455, 1), numel(Dive));
end

minTime = min(diveTimes);
maxTime = max(diveTimes);

windowDays = 4;  % 3-day window
stepDays = 1;    % 1-day step

currentStart = minTime;

while currentStart + days(windowDays - 1) <= maxTime
    currentEnd = currentStart + days(windowDays - 1);
    
    % Find dives in current window
    inWindow = (diveTimes >= currentStart) & (diveTimes <= currentEnd);
    idxInWindow = find(inWindow);
    
    if isempty(idxInWindow)
        currentStart = currentStart + days(stepDays);
        continue;
    end
    
    clf;
    
    % Get data for window (now rows are dives)
    windowProfiles = mN_455(idxInWindow, :)';  % Transpose for visualization
    
    % Get time-of-day labels
    windowTOD = {Dive(idxInWindow).tod};
    
    % Normalization (adjust clim as needed)
    clim = [0, 10000];
    normVals = (windowProfiles - clim(1)) / diff(clim);
    normVals = min(max(normVals, 0), 1);
    
    % Create colored image
    [depthBins, numDives] = size(windowProfiles);
    RGB = zeros(depthBins, numDives, 3);
    
    % Colormaps
    dayMap = hot(256);
    nightMap = bone(256);
    
    for col = 1:numDives
        if strcmpi(windowTOD{col}, 'Day')
            cmap = dayMap;
        else
            cmap = nightMap;
        end
        idx = round(normVals(:, col) * 255) + 1;
        RGB(:, col, :) = cmap(idx, :);
    end
    
    % Plot
    ax = subplot(1,1,1);
    image(RGB)
    set(gca, 'YDir', 'reverse')
    yticks(10:10:90)
    yticklabels(string(10:10:90))
    ylabel('Depth (m)')
    xlabel('Dive Number')
    title(sprintf('455 kHz Profiles: %s to %s', ...
        datestr(currentStart, 'yyyy-mm-dd'), ...
        datestr(currentEnd, 'yyyy-mm-dd')))
    
    % Set x-axis ticks
    tickStep = max(1, floor(numDives / 10));
    xticks(1:tickStep:numDives)
    xticklabels(idxInWindow(1:tickStep:end))
    
    % Colorbars
    axPos = get(ax, 'Position');
    
    % Day colorbar (left)
    cLeft = colorbar('Position', [axPos(1)-0.04, axPos(2), 0.02, axPos(4)]);
    colormap(cLeft, dayMap)
    cLeft.Ticks = linspace(0,1,5);
    cLeft.TickLabels = compose('%.0f', linspace(clim(1), clim(2), 5));
    cLeft.Label.String = 'Daytime';
    
    % Night colorbar (right)
    cRight = colorbar('Position', [axPos(1)+axPos(3)+0.02, axPos(2), 0.02, axPos(4)]);
    colormap(cRight, nightMap)
    cRight.Ticks = linspace(0,1,5);
    cRight.TickLabels = compose('%.0f', linspace(clim(1), clim(2), 5));
    cRight.Label.String = 'Nighttime';
    
    drawnow;
    pause(1);  % Adjust pause time as needed
    
    currentStart = currentStart + days(stepDays);
end

%% Unified Animation
figure('Position', [100 100 1800 900], 'Color', 'w');

% 1. Setup Map Subplot (Left)
subplot(1,2,1);
% Bathymetry background
scatter(lon_subset(:), lat_subset(:), 10, -double(bathy_subset(:)), 'filled');
hold on;
h_fulltrack = plot(gliderdataH.longitude, gliderdataH.latitude, 'k-', 'LineWidth', 0.5);
h_glider = plot(gliderdataH.longitude(1), gliderdataH.latitude(1), 'r-', 'LineWidth', 2);
h_pos = plot(gliderdataH.longitude(1), gliderdataH.latitude(1), 'ro', 'MarkerFaceColor', 'y', 'MarkerSize', 8);
h_time = text(min(gliderdataH.longitude), max(gliderdataH.latitude)+0.1, ...
    ['Time: ' datestr(gliderdataH.localtime(1))], 'FontSize', 12, 'BackgroundColor', 'w');
colormap(gca, parula);
h_cb1 = colorbar;
h_cb1.Label.String = 'Bathymetry (m)';
axis equal tight;
title('Glider Position');

% 2. Setup Prey Subplot (Right)
subplot(1,2,2);
% Initial prey data
diveTimes = [Dive.localtime];
currentStart = min(diveTimes);
currentEnd = currentStart + days(4);
inWindow = (diveTimes >= currentStart) & (diveTimes <= currentEnd);
idx = find(inWindow);

% Matrix orientation check
if size(mN_455,1) ~= numel(Dive), mN_455 = mN_455'; end

% Create initial prey plot
h_preyImg = imagesc(1:sum(inWindow), 1:size(mN_455,2), mN_455(idx,:)');
set(gca, 'YDir', 'reverse');
colormap(gca, hot);
caxis([0 10000]);
h_preyTitle = title(sprintf('Prey Data: %s to %s', datestr(currentStart), datestr(currentEnd)));
xlabel('Dive Number'); 
ylabel('Depth (m)');
h_cb2 = colorbar;
h_cb2.Label.String = 'Prey Density';

% 3. Animation Loop with Error Handling
try
    for k = 2:length(gliderdataH.longitude)
        % Update glider position
        subplot(1,2,1);
        set(h_glider, 'XData', gliderdataH.longitude(1:k), 'YData', gliderdataH.latitude(1:k));
        set(h_pos, 'XData', gliderdataH.longitude(k), 'YData', gliderdataH.latitude(k));
        set(h_time, 'String', ['Time: ' datestr(gliderdataH.localtime(k))]);
        
        % Update prey plot every 4 days
        if gliderdataH.localtime(k) >= currentEnd
            currentStart = currentStart + days(4);
            currentEnd = currentEnd + days(4);
            inWindow = (diveTimes >= currentStart) & (diveTimes <= currentEnd);
            idx = find(inWindow);
            
            subplot(1,2,2);
            set(h_preyImg, 'CData', mN_455(idx,:)');
            set(h_preyTitle, 'String', sprintf('Prey Data: %s to %s', ...
                datestr(currentStart), datestr(currentEnd)));
            drawnow limitrate; % Smoother animation
        end
        
        pause(0.1); % Control speed
    end
catch ME
    disp('Animation stopped:');
    disp(ME.message);
end

%% Remove test mission data at the beginning

% --- Remove Dive data before 17 July 2020 ---
cutoffDate = datetime(2020,7,17);

% Filter Dive struct
keepIdx = [Dive.localtime] >= cutoffDate;
Dive = Dive(keepIdx);

% Filter gliderdataH (assuming it's struct with vectors)
keepIdxH = gliderdataH.localtime >= cutoffDate;
gliderdataH.longitude = gliderdataH.longitude(keepIdxH);
gliderdataH.latitude  = gliderdataH.latitude(keepIdxH);
gliderdataH.localtime = gliderdataH.localtime(keepIdxH);

%% Enhanced Unified Animation with Time-Lagged Prey Updates
figure('Position', [100 100 2000 800], 'Color', 'w');
set(groot, 'defaultAxesFontSize', 14);
set(groot, 'defaultTextFontSize', 14);
set(groot, 'defaultAxesLabelFontSizeMultiplier', 1.2);

% 1. Setup Map Subplot (Left)
subplot(1,2,1);
% Bathymetry background
scatter(lon_subset(:), lat_subset(:), 10, -double(bathy_subset(:)), 'filled');
hold on;
h_fulltrack = plot(gliderdataH.longitude, gliderdataH.latitude, 'k-', 'LineWidth', 0.5);
h_glider = plot(gliderdataH.longitude(1), gliderdataH.latitude(1), 'r-', 'LineWidth', 2);
h_pos = plot(gliderdataH.longitude(1), gliderdataH.latitude(1), 'ro', 'MarkerFaceColor', 'y', 'MarkerSize', 8);
ax_pos = get(gca, 'Position');
h_time = annotation('textbox',...
    [ax_pos(1), ax_pos(2)+ax_pos(4)-0.11, 0.11, 0.035],...
    'String', ['Date: ' datestr(gliderdataH.localtime(1), 'dd-mmm-yyyy')],...
    'FontSize', 14,...
    'BackgroundColor', 'w',...
    'EdgeColor', 'k',...
    'Margin', 2,...
    'VerticalAlignment', 'top',...
    'HorizontalAlignment', 'left');
blue_colormap = [
    linspace(0.6, 0.0, 256)', ...
    linspace(0.9, 0.2, 256)', ...
    linspace(1.0, 0.6, 256)'
];
colormap(gca, blue_colormap);
h_cb1 = colorbar;
h_cb1.FontSize = 14;
h_cb1.Label.String = 'Bathymetry (m)';
axis equal tight;
title('Glider Track in the swGSL');

% Initialize segment tracking
segment_colors = {'red', [0.6 0 0]}; % Red and dark red
current_segment = 1;
segment_start_idx = 1;
segment_plots = [];

% 2. Setup Prey Subplot (Right)
subplot(1,2,2);
% Initialize with empty prey plot
h_preyImg = imagesc([], [], []);
set(gca, 'YDir', 'reverse', 'Box', 'off');
colormap(gca, hot);
caxis([0 10000]);
ylim([1 100]);
h_preyTitle = title('Echosounder Data: Waiting for first 4-day window...');
xlabel('Glider Profile Number'); 
ylabel('Depth (m)');
h_cb2 = colorbar;
h_cb2.FontSize = 14;
h_cb2.Label.String = 'Zooplankton Density (ind/m3)';

% Initialize prey data tracking
prey_windows = {};
window_complete = false;
current_window_start = min([Dive.localtime]);

% 3. Enhanced Animation Loop
try
    for k = 2:length(gliderdataH.longitude)
        % Update glider position
        subplot(1,2,1);
        set(h_glider, 'XData', gliderdataH.longitude(1:k), 'YData', gliderdataH.latitude(1:k));
        set(h_pos, 'XData', gliderdataH.longitude(k), 'YData', gliderdataH.latitude(k));
        set(h_time, 'String', ['Date: ' datestr(gliderdataH.localtime(k), 'dd-mmm-yyyy')]);
        
        % Check for new 4-day segment
        if gliderdataH.localtime(k) >= current_window_start + days(4)
            % Finalize current segment
            segment_plots(end+1) = plot(gliderdataH.longitude(segment_start_idx:k), ...
                                      gliderdataH.latitude(segment_start_idx:k), ...
                                      'Color', segment_colors{mod(current_segment-1,2)+1}, ...
                                      'LineWidth', 2);
            
            % Store completed window
            window_end = current_window_start + days(4);
            prey_windows{end+1} = {current_window_start, window_end};
            
            % Prepare next segment
            current_segment = current_segment + 1;
            segment_start_idx = k;
            current_window_start = window_end;
            window_complete = true;
        end
        
        % Update prey plot only after window completes
        if window_complete && ~isempty(prey_windows)
            % Get the oldest unshown window
            window_to_show = prey_windows{1};
            prey_windows(1) = [];
            
            % Find prey data in this window
            inWindow = ([Dive.localtime] >= window_to_show{1}) & ...
                      ([Dive.localtime] <= window_to_show{2});
            idx = find(inWindow);
            
            if any(inWindow)
                subplot(1,2,2);
                set(h_preyImg, 'XData', 1:sum(inWindow), ...
                              'YData', 1:size(mN_455,2), ...
                              'CData', mN_455(idx,:)');
                set(h_preyTitle, 'String', sprintf('Echosounder Data: %s to %s', ...
                            datestr(window_to_show{1}, 'dd-mmm-yyyy'), ...
                            datestr(window_to_show{2}, 'dd-mmm-yyyy')));
            else
                set(h_preyTitle, 'String', sprintf('No data: %s to %s', ...
                            datestr(window_to_show{1}, 'dd-mmm-yyyy'), ...
                            datestr(window_to_show{2}, 'dd-mmm-yyyy')));
            end
            window_complete = false;
        end
        
        pause(0.05); % Control speed
        drawnow limitrate;
    end
    
    % Final segment
    if segment_start_idx < length(gliderdataH.longitude)
        segment_plots(end+1) = plot(gliderdataH.longitude(segment_start_idx:end), ...
                                  gliderdataH.latitude(segment_start_idx:end), ...
                                  'Color', segment_colors{mod(current_segment-1,2)+1}, ...
                                  'LineWidth', 2);
    end
    
catch ME
    disp('Animation stopped:');
    disp(ME.message);
end

%% Save video
% Enhanced Unified Animation with Video Recording
figure('Position', [100 100 2000 800], 'Color', 'w');
set(groot, 'defaultAxesFontSize', 14);
set(groot, 'defaultTextFontSize', 14);
set(groot, 'defaultAxesLabelFontSizeMultiplier', 1.2);

% Initialize video writer
videoFile = 'glider_prey_animation.mp4';
v = VideoWriter(videoFile, 'MPEG-4');
v.Quality = 100;
open(v);

% ====================== 1. MAP SUBPLOT ======================
subplot(1,2,1);

% Vibrant blue colormap (teal to deep blue)
blue_map = [
    linspace(0.6, 0.0, 256)', ...
    linspace(0.9, 0.2, 256)', ...
    linspace(1.0, 0.6, 256)'
];

% Bathymetry plot
scatter(lon_subset(:), lat_subset(:), 10, -double(bathy_subset(:)), 'filled');
colormap(gca, blue_map);
hold on;

% Glider track elements
h_fulltrack = plot(gliderdataH.longitude, gliderdataH.latitude, 'k-', 'LineWidth', 0.5);
h_glider = plot(gliderdataH.longitude(1), gliderdataH.latitude(1), 'r-', 'LineWidth', 2);
h_pos = plot(gliderdataH.longitude(1), gliderdataH.latitude(1), 'ro', 'MarkerFaceColor', 'y', 'MarkerSize', 8);

% Date annotation
ax_pos = get(gca, 'Position');
h_time = annotation('textbox',...
    [ax_pos(1), ax_pos(2)+ax_pos(4)-0.11, 0.11, 0.035],...
    'String', ['Date: ' datestr(gliderdataH.localtime(1), 'dd-mmm-yyyy')],...
    'FontSize', 14,...
    'BackgroundColor', 'w',...
    'EdgeColor', 'k',...
    'Margin', 2,...
    'VerticalAlignment', 'top');

h_cb1 = colorbar;
h_cb1.FontSize = 14;
h_cb1.Label.String = 'Bathymetry (m)';
axis equal tight;
title('Glider Track in the swGSL');

% ====================== 2. PREY SUBPLOT ======================
subplot(1,2,2);
h_preyImg = imagesc([], [], []);
set(gca, 'YDir', 'reverse', 'Box', 'off');
colormap(gca, hot);
caxis([0 10000]);
ylim([1 100]);
h_preyTitle = title('Echosounder Data: Waiting for first 4-day window...');
xlabel('Glider Profile Number'); 
ylabel('Depth (m)');
h_cb2 = colorbar;
h_cb2.FontSize = 14;
h_cb2.Label.String = 'Zooplankton Density (ind/m3)';

% ====================== ANIMATION LOOP ======================
segment_colors = {'red', [0.6 0 0]};
current_segment = 1;
segment_start_idx = 1;
prey_windows = {};
window_complete = false;
current_window_start = min([Dive.localtime]);

try
    for k = 2:length(gliderdataH.longitude)
        % Update glider position
        subplot(1,2,1);
        set(h_glider, 'XData', gliderdataH.longitude(1:k), 'YData', gliderdataH.latitude(1:k));
        set(h_pos, 'XData', gliderdataH.longitude(k), 'YData', gliderdataH.latitude(k));
        set(h_time, 'String', ['Date: ' datestr(gliderdataH.localtime(k), 'dd-mmm-yyyy')]);
        
        % Handle 4-day segments
        if gliderdataH.localtime(k) >= current_window_start + days(4)
            segment_plots(end+1) = plot(gliderdataH.longitude(segment_start_idx:k), ...
                                      gliderdataH.latitude(segment_start_idx:k), ...
                                      'Color', segment_colors{mod(current_segment-1,2)+1}, ...
                                      'LineWidth', 2);
            
            prey_windows{end+1} = {current_window_start, current_window_start + days(4)};
            current_segment = current_segment + 1;
            segment_start_idx = k;
            current_window_start = current_window_start + days(4);
            window_complete = true;
        end
        
        % Update prey plot
        if window_complete && ~isempty(prey_windows)
            window_to_show = prey_windows{1};
            prey_windows(1) = [];
            
            inWindow = ([Dive.localtime] >= window_to_show{1}) & ...
                      ([Dive.localtime] <= window_to_show{2});
            idx = find(inWindow);
            
            subplot(1,2,2);
            if any(inWindow)
                set(h_preyImg, 'XData', 1:sum(inWindow), 'YData', 1:size(mN_455,2), 'CData', mN_455(idx,:)');
                set(h_preyTitle, 'String', sprintf('Echosounder Data: %s to %s', ...
                                datestr(window_to_show{1}, 'dd-mmm-yyyy'), ...
                                datestr(window_to_show{2}, 'dd-mmm-yyyy')));
            else
                set(h_preyTitle, 'String', sprintf('No data: %s to %s', ...
                                datestr(window_to_show{1}, 'dd-mmm-yyyy'), ...
                                datestr(window_to_show{2}, 'dd-mmm-yyyy')));
            end
            window_complete = false;
        end
        
        % Capture frame
        writeVideo(v, getframe(gcf));
        drawnow;
        pause(0.05); % Your preferred animation speed
    end
    
    % Final cleanup
    if segment_start_idx < length(gliderdataH.longitude)
        segment_plots(end+1) = plot(gliderdataH.longitude(segment_start_idx:end), ...
                                  gliderdataH.latitude(segment_start_idx:end), ...
                                  'Color', segment_colors{mod(current_segment-1,2)+1}, ...
                                  'LineWidth', 2);
        writeVideo(v, getframe(gcf));
    end
    
    close(v);
    disp(['Video saved: ' fullfile(pwd, videoFile)]);
    
catch ME
    close(v);
    disp(['Animation stopped: ' ME.message]);
end