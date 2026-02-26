%% Tentative code for removing seafloor echoes and values below the seafloor Part 2/2

function [sv_nb] = remove_bottom_AM(sv, bott_sv, bott_ind, buffer_bins)

% Default buffer: bins above the bottom
if nargin < 4
    buffer_bins = 10;
end

sv_nb = NaN(size(sv));

for i = 1:length(bott_ind)
    this_ind = bott_ind(i);

    if isnumeric(this_ind) && ~isnan(this_ind) && this_ind > buffer_bins
        cutoff = this_ind - buffer_bins;
        sv_nb(1:cutoff, i) = sv(1:cutoff, i);
    else       
        sv_nb(:, i) = sv(:, i);
    end
end

end