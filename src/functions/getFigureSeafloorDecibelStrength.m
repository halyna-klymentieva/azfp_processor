function getFigureSeafloorDecibelStrength(azfpData)
%GETFIGURESEAFLOORDECIBELSTRENGTH Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    azfpData
end

% % Find the seafloor indices
for k=1:length(azfpData) % for each frequency
    % for j = 1:cc % for each dive
        % first need to pull out sv and depths for each dive separately, so
        % remove bottom does not delete too much data
        sv = azfpData(k).Sv';
        depth = azfpData(k).Depth';
        % depth = azfpData(k).PingDepth';

        % next run the find bottom function
        [bott_sv2(k,:), bott_dep2(k,:), bott_ind2(k,:)] = find_bottom_AM(sv, depth);
    % end
end

figure(1)

% subplot(2,2,1)
histogram(bott_sv2(1,:),'BinWidth',1,'FaceAlpha',0.5);
hold on
[values, edges] = histcounts(bott_sv2(1,:),'BinWidth', 1);
centers = (edges(1:end-1)+edges(2:end))/2;
plot(centers, values,'LineWidth',2)
% hold off
% title('130kHz')
%
% subplot(2,2,2)
histogram(bott_sv2(2,:), 'BinWidth',1,'FaceAlpha',0.5);
hold on
[values, edges] = histcounts(bott_sv2(2,:),'BinWidth', 1);
centers = (edges(1:end-1)+edges(2:end))/2;
plot(centers, values,'LineWidth',2)
% hold off
% title('200kHz')
%
% subplot(2,2,3)
histogram(bott_sv2(3,:),'BinWidth',1,'FaceAlpha',0.5);
hold on
[values, edges] = histcounts(bott_sv2(3,:),'BinWidth', 1);
centers = (edges(1:end-1)+edges(2:end))/2;
plot(centers, values,'LineWidth',2)
% hold off
% title('455kHz')
% %
% % subplot(2,2,4)
histogram(bott_sv2(4,:),'BinWidth',1,'FaceAlpha',0.5);
hold on
[values, edges] = histcounts(bott_sv2(4,:),'BinWidth', 1);
centers = (edges(1:end-1)+edges(2:end))/2;
plot(centers, values,'LineWidth',2)
% % hold off
% % title('769kHz')
%
% sgtitle(strcat('Bottom Depth Sv Values (dB) for ',date, ' July 2022'));
legend({'130 kHz','','200 kHz','','455 kHz','','769 kHz'});

% % % filename = strcat("/Users/dmossman/Box/2022 MSc Thesis Work/Visuals/MATLAB Echosounder Figures/",date,"Sept/Bott_Sv_Hist_Subplots",date,"Sept.png");
% % % print(gcf,'-dpng',filename,'-r0')
% % % clear filename;

end