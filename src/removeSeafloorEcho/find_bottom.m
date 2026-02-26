% Pull the bottom out of the echosounder data

function [bott_sv, bott_dep, bott_ind] = find_bottom(sv,depth)

%sv=varmat;
% depth=echo.depth;
% Allocate space for Sv and depth values associated with bottom (1 value
% per dive)
bott_sv = NaN*ones(1,size(sv,2));
bott_dep = NaN*ones(1,size(sv,2)); 
bott_ind = NaN*ones(1,size(sv,2));

for i = 1:size(sv,2) % 1:number of dives
    
    [xmax, imax] = extrema(sv(:,i)); % Local peak values of Sv and their row indices for the i'th column (dive), extrema ignores NaNs (unlike findpeaks, which ignores peaks if a NaN is within 2 indices of it)
    
    % [SEPARATE THE STATEMENTS BELOW SO THEY YIELD THE INDEX IN THE SV FOR
    % THE BOTTOM]
    
    % [JULY 22,2016 NOTE: THIS SCRIPT WAS WRITTEN BEFORE THE APPLICATION OF
    % THE 2015 CALIBRATION CORRECTION OF +2 DB RELATIVE TO THE -4 DB
    % APPLIED IN THE SCRIPT THAT DOWNLOADS AND PRE-PROCESSES THE DATA; THE
    % VALUES BELOW NOW REFLECT THAT CHANGE = +2 DB ADDED TO THE THRESHOLD
    % VALUES TO MATCH THE OFFSET CREATED BY THE CORRECTION IN THE SV DATA]
    
    % Get value(s) at or above -22 dB intensity and deeper than 45 m
    max_vals = find(xmax >= -20); %Within vector indices for values that are above -22 dB
    max_deps = imax(max_vals); % Depth index values that have intensites greater than -22 dB
    max_depinds = find(max_deps >= 90); % Indices in the max deps that are above 90 bins
    max_deps = max_deps(max_depinds); % Reassign the variable so that it only retains the depths that are above 45 m and that are above -22 dB
    
    if length(max_deps) == 1
        max_ind = max_deps; % If there is only one possible value for the bottom, take the first value
        bott_sv(i) = sv(max_ind,i);
        bott_dep(i) = depth(max_ind,i);
        bott_ind(i) = max_ind; 
        
    elseif isempty(max_deps) % If it is an early dive that does not contain a bottom signal...
        bott_sv(i) = NaN;
        bott_dep(i) = NaN;
        bott_ind(i) = NaN;  
    
    elseif length(max_deps) > 1 % If there are more than one value that could be the bottom... 
        max_ind = min(max_deps); % The shallowest value from among those that are greater than -22 dB intensity and that is deeper than 45 m
        bott_sv(i) = sv(max_ind,i);
        bott_dep(i) = depth(max_ind,i);
        bott_ind(i) = max_ind;
    end
    
end
    
disp('Warning: Plot the bottom depths obtained and compare to the dataset, several may be outliers than need to be changed to NaN or a user specified value.')

end


%     [xmax imax] = extrema(sv(:,i)); % Local peak values of Sv and their row indices for the i'th column (dive), extrema ignores NaNs (unlike findpeaks, which ignores peaks if a NaN is within 2 indices of it)
%     pk_ind = find(xmax >= -22); % Get xmax above -22 dB from the peaks
%     if ~isempty(pk_ind)
%         pk_ind = min(pk_ind); % Select the shallowest peak
%         pk_loc = imax(pk_ind); % Assign the variable index associated with the bottom peak 
%         bott_sv(i)=sv(locs(pk_ind),i);
%         bott_dep(i)=depth(locs(pk_ind),i);
%         bott_ind(i)=locs(pk_ind);
%     end
% % Archived in case the extrema test doesn't work ... just put findpeaks,
% % locs, pks back into the code.
% end



% Need to reduce spikiness:
% Take out the dives that have an incorrect (too shallow or too deep) bottom depth
% (these are early dives where the echosounder didn't see the bottom or where noise ...
%  below the bottom is identified as the bottom)
% Averaging? --> Compare to the i-1 and i+1 depth and if it is more than a
% threshold value above that average, the new value IS the average.

% [If statement that contains a threshold value (from correct bottom avg) that the bottom is
% consistently AND is more than 4-5 bins (not fish).]
