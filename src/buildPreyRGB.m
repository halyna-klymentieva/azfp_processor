%% Plot densities as a video

function RGB = buildPreyRGB(DiveSubset)
    % Parameters
    clim = [0, 10000];   % color scale
    hotMap  = hot(256);
    boneMap = bone(256);

    % Collect all profiles and TOD
    allProfiles = [];
    todLabels   = [];
    for i = 1:length(DiveSubset)
        % Select mN_455 (assuming index 1 is 455 kHz, adjust if needed)
        profile = DiveSubset(i).P(1).mN(:);
        allProfiles = [allProfiles, profile];
        todLabels   = [todLabels, string(DiveSubset(i).tod)];
    end

    % Normalize values to 0-1
    normVals = (allProfiles - clim(1)) / (clim(2)-clim(1));
    normVals = min(max(normVals,0),1);

    idxVals = round(normVals * 255) + 1;
    [rows, cols] = size(allProfiles);
    RGB = ones(rows, cols, 3);

    % Assign colors column by column
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
end