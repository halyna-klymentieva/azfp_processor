function output = filterSurfaceNoise(azfpData, dives, depthMin, depthMax)
%FILTERSURFACENOISE Remove surface noise

divesLength = length(dives);
for k = 1:4
    sv_data = azfpData(k).Sv';
    depth_data = azfpData(k).Depth';

    for d = 95:divesLength
        idx = dives(d).Index;

        % Get ping range for this dive
        ping_start = idx(1);
        ping_end = idx(2);

        % Mask depth band
        for ping = ping_start:ping_end
            depth_col = depth_data(:, ping);
            mask_band = depth_col >= depthMin & depth_col <= depthMax;
            sv_data(mask_band, ping) = NaN;
        end
    end

    % Assign back to Output structure
    azfpData(k).Sv = sv_data';
end

output = azfpData;
end