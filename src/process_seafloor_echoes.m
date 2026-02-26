%% Tentative streamlined workflor for removing seafloor echoes

function [Output] = process_seafloor_echoes(Output)
% PROCESS_SEAFLOOR_ECHOES - Handles mismatched ping counts between Sv and Range
% Input/Output: Output structure with:
%   Sv: [n_pings × n_range_bins] 
%   Range: [n_depth_cells × n_range_bins] (typically smaller)

for f = 1:4 % For each frequency
    % Get data dimensions
    sv_data = Output(f).Sv; % [15085 × 742]
    range_data = Output(f).Range; % [194 × 742]
    
    % Initialize results
    n_pings = size(sv_data, 1);
    n_range_bins = size(sv_data, 2);
    bott_ind = NaN(n_pings, 1);
    sv_clean = sv_data;
    
    % Frequency parameters
    switch f
        case 1 % 130 kHz
            sv_thresh = -60; min_depth = 30; buffer = 5;
        case 2 % 200 kHz
            sv_thresh = -58; min_depth = 35; buffer = 5;
        case 3 % 455 kHz
            sv_thresh = -55; min_depth = 40; buffer = 6;
        case 4 % 769 kHz
            sv_thresh = -50; min_depth = 40; buffer = 8;
    end
    
    % Process each ping using nearest range profile
    for p = 1:n_pings
        % Select range profile (cycling through available range profiles)
        range_idx = mod(p-1, size(range_data,1)) + 1;
        range_profile = range_data(range_idx, :);
        
        % Bottom detection
        [~, ~, bi] = extrema_bottom(sv_data(p,:)', range_profile', sv_thresh, min_depth);
        
        if ~isnan(bi)
            % Apply bottom mask
            cutoff = bi - buffer;
            if cutoff > 0
                sv_clean(p, cutoff+1:end) = NaN;
                
                % Additional cleaning for 769 kHz
                if f == 4
                    strong_echoes = (1:n_range_bins) > (cutoff-10) & sv_data(p,:) > -45;
                    sv_clean(p, strong_echoes) = NaN;
                end
            end
            bott_ind(p) = bi;
        end
    end
    
    % Store results
    Output(f).bott_ind = bott_ind;
    Output(f).Sv = sv_clean;
end
end

function [bott_sv, bott_dep, bott_ind] = extrema_bottom(sv, depth, sv_thresh, min_depth)
% EXTREMA_BOTTOM - Safe bottom detection
    bott_sv = NaN;
    bott_dep = NaN;
    bott_ind = NaN;
    
    valid = ~isnan(sv) & ~isnan(depth);
    if sum(valid) < 10 % Need at least 10 valid samples
        return
    end
    
    [xmax, imax] = extrema(sv(valid));
    imax_orig = find(valid);
    imax = imax_orig(imax); % Map back to original indices
    
    valid_peaks = xmax >= sv_thresh & depth(imax) > min_depth;
    if any(valid_peaks)
        [~, idx] = max(xmax(valid_peaks)); % Strongest peak
        bott_ind = imax(valid_peaks(idx));
        bott_sv = sv(bott_ind);
        bott_dep = depth(bott_ind);
    end
end