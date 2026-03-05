function output = filterSeafloorEchoes(azfpData, dives)
%FILTERSEAFLOORECHOES Remove seafloor echoes (Step 1 in the tutorial)

%% Find the seafloor indices
fprintf('Filtering seafloor echoes...\n');
tic
for k = 1:length(azfpData) % for each frequency
    %for j = 1:cc % for each dive
    % first need to pull out sv and depths for each dive separately, so
    % remove bottom does not delete too much data
    sv = azfpData(k).Sv';
    depth = azfpData(1).Depth';
    % depth = azfpData(1).PingDepth';
    % next run the find and remove bottom functions to delete all data
    % at or below the seafloor
    [bott_sv, ~, bott_ind] = find_bottom_AM(sv, depth); %June 2025: Andréa updated
    sv_nb = remove_bottom_AM(sv, bott_sv, bott_ind); %with new code
    % finally, reinsert the data with the seafloor removed into the
    % original Output file
    % azfpData(k).Sv(dives(j).Index(1):dives(j).Index(2),:) = sv_nb';
    azfpData(k).Sv = sv_nb';
    %end
end
toc

%% Remove residual echoes

for k = 1:length(azfpData)
    sv_data = azfpData(k).Sv';
    depth_data = azfpData(k).Depth';

    % Frequency-specific parameters
    switch k
        case 4 % 769 kHz needs different handling
            depth_range = [50, 98]; % Detection zone (m)
            buffer = 16; % Bins above peak to keep
        otherwise
            depth_range = [50, 98];
            buffer = 10;
    end

    % Create ping-level mask only for dives 1–227
    dive_mask = false(1, size(sv_data, 2));
    for d = 1:length(dives)
        idx = dives(d).Index;

        if isnumeric(idx) && numel(idx) == 2
            ping_start = idx(1);
            ping_end = idx(2);
        elseif isstruct(idx) && isfield(idx, 'StartPing') && isfield(idx, 'EndPing')
            ping_start = idx.StartPing(1);
            ping_end = idx.EndPing(1);
        else
            continue
        end

        % Apply bounds check
        ping_end = min(ping_end, size(sv_data, 2));
        dive_mask(ping_start:ping_end) = true;
    end

    for ping = find(dive_mask)
        sv_col = sv_data(:, ping);
        depth_col = depth_data(:, ping);

        if all(isnan(sv_col)) || max(depth_col) < depth_range(1)
            continue
        end

        % --- STAGE 1: Peak detection within target range ---
        zone_mask = depth_col > depth_range(1) & depth_col < depth_range(2);
        [~, peak_idx] = max(sv_col.*zone_mask);

        % --- STAGE 2: Masking logic ---
        if ~isnan(peak_idx) && sv_col(peak_idx) > -70
            cutoff = min(length(sv_col), peak_idx+buffer);
            sv_col(peak_idx:cutoff) = NaN;
            sv_col(cutoff+1:end) = NaN;
        elseif max(depth_col) > 105
            sv_col(depth_col > 105) = NaN;
        end

        sv_data(:, ping) = sv_col;
    end

    azfpData(k).Sv = sv_data';
end


output = azfpData;
end