%% Tentative code for removing seafloor echoes and values below the seafloor Part 1/2

function [bott_sv, bott_dep, bott_ind] = find_bottom_AM(sv, depth)

bott_sv = NaN*ones(1, size(sv, 2));
bott_dep = NaN*ones(1, size(sv, 2)); 
bott_ind = NaN*ones(1, size(sv, 2));

for i = 1:size(sv, 2)

    sv_col = sv(:, i);
    dep_col = depth(:, i);

    if all(isnan(sv_col)) || all(isnan(dep_col))
        continue
    end

    % Identify peaks using extrema
    [xmax, imax] = extrema(sv_col);

    % Filter peaks above a backscatter threshold
    strong_peaks = xmax >= -55;
    strong_inds = imax(strong_peaks);

    % Further filter:
    if ~isempty(strong_inds)
        deep_inds = strong_inds(dep_col(strong_inds) > 45);
    else
        deep_inds = [];
    end

    if isempty(deep_inds)
        continue
    end

    % Choose the shallowest strong & deep peak as seafloor
    ind = min(deep_inds);

    bott_sv(i) = sv_col(ind);
    bott_dep(i) = dep_col(ind);
    bott_ind(i) = ind;

end

end