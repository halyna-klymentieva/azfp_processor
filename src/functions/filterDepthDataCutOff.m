function output = filterDepthDataCutOff(azfpData, dives, depth)
%FILTERDEPTHDATACUTOFF Remove all data below a defined depth

arguments (Input)
    azfpData
    dives
    depth
end

arguments (Output)
    output
end

% The purpose of this count in the second cycle is unclear.
divesCount = length(dives);

for k = 1:length(azfpData)
    sv_data = azfpData(k).Sv';
    depth_data = azfpData(k).Depth';

    dive_mask = false(size(depth_data));

    % Define dives
    for d = divesCount:divesCount
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

        if ping_end <= size(dive_mask, 2)
            dive_mask(:, ping_start:ping_end) = true;
        end
    end

    % Define depth
    deep_mask = depth_data > depth;

    % Combine dive + depth masks
    final_mask = dive_mask & deep_mask;

    % Apply mask
    sv_data(final_mask) = NaN;

    % Store Sv data back in Output
    azfpData(k).Sv = sv_data';
end


output = azfpData;
end
