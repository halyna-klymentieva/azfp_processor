%% Separate the  bathymetry datasets into transits of the basin and apply
% them to remove the bottom signal from the sv data.

% Beginning of transit 1 = dive #286
% End of transit 1 = dive #1490
% Beginning of transit 2 = dive #1794
% End of transit 2 = dive #2705

% Taking only the dives in the above ranges removes the sections where the
% glider was reaching the first waypoint at the beginning of the survey,
% while the glider repositioned itself from the last waypoint back to the
% first for the second transit, and from the last waypoint while the glider
% began its return flight to Halifax.

% These ranges can also be used to pull out the cut data, if that becomes
% of interest.
% 
% load bathydataFull
% 
% Separate the bathy data into transits
% 
% bott_sv1 = bott_sv(~isnan(bott_sv));
% bott_dep1 = bott_dep(~isnan(bott_dep));
% bott_ind1 = bott_ind(~isnan(bott_dep)); % These are the row indices for each dive
% 
% bott_sv2 = bott_sv(1794:2705);
% bott_dep2 = bott_dep(1794:2705);
% bott_ind2 = bott_ind(1794:2705);

% Use the bathy data to remove the bottom signal from the Sv data for both
% transits

% Need to use the ind values to put in the rows that each column goes down
% to, and leave the ones without a bottom value alone

% sv1_nb = nan*ones(size(sv1));
% dep1_nb = nan*ones(size(sv1));
% lat1_nb = nan*ones(size(sv1));
% lon1_nb = nan*ones(size(sv1));
% time1_nb = nan*ones(size(sv1));

%
function [sv_nb] = removeBottom(sv,bott_sv,bott_ind)

for i = 1:length(bott_sv)
    if ~isnan(bott_ind(i))
    tempdata = sv(1:bott_ind(i)-7,i); % Take all data from 3 m (= 6 bins) above the bottom bin (potentially increase to 5 m?)
    nanpad = nan*ones(length(sv(:,i))-length(tempdata),1); % Number of nans that need to be added to make the matrix the same size at the original
    tempdata = [tempdata;nanpad];
    sv_nb(:,i) = tempdata;
    else % If there is not bottom echo in the dive, still need to remove noise?
%         revert_ct = i-1; % Want to get the bottom depth from the previous dive and use it to cut the early dive
%         tempdata = sv1(1:bott_ind1(revert_ct)-7,i);
%         nanpad = nan*ones(length(sv1(:,i))-length(tempdata),1); % Number of nans that need to be added to make the matrix the same size at the original
%         tempdata = [tempdata;nanpad];
         sv_nb(:,i) = sv(:,i); % Try leaving the data as is in the new matrix
    end
    
end
end


% % [NEED TO DEBUG THE TRANSIT 1 CODE BEFORE USING THIS]
% % And for transit 2:
% for i = 1:length(bott_sv2)
%     if ~isnan(bott_ind2(i))
%     tempdata = sv2(1:bott_ind2(i)-7,i);
%     nanpad = nan*ones(length(sv2(:,i))-length(tempdata),1);
%     tempdata = [tempdata;nanpad];
%     sv2_nb(:,i) = tempdata;
%     else
%         sv2_nb(:,i) = sv2(:,i);
%     end
% end

% save transit1_nb sv1_nb dep1_nb lat1_nb lon1_nb time1_nb
% save transit2_nb sv2_nb dep2_nb lat2_nb lon2_nb time2_nb