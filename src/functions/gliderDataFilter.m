function output = gliderDataFilter(gliderdata, azfpData)
%GLIDERDATAFILTER Filtering glider data for pings during inactivity
%
% During the Baffin Bay 2023 and 2024 missions, it appears that the echosounder was
% pinging even when the glider recorded the AZFP as being powered off.
% This means there will be random pings on the upcast or when the glider is
% transmitting data at the surface.  I need to remove these pings.  Keep
% only pings and time stamps where the glider has recorded that the AZFP is
% on.  Talk to Jude and ASL about this issue.
arguments (Input)
    gliderdata
    azfpData
end

arguments (Output)
    output
end

% First test whether there is a time alignment issue using the indicator that the
% AZFP was on (c_azfp_on).  -1 = off, 0 = on
nanindex2=find(~isnan(gliderdata.c_azfp_on));
azfp_on=gliderdata.c_azfp_on(nanindex2);
azfp_time=gliderdata.time(nanindex2);

% figure
% scatter(azfpData(1).Date,zeros(1,length(azfpData(1).Date)))
% hold on
% plot(azfp_time,azfp_on)
% set(gca,'YLim',[-2 1])
% legend('AZFP time stamp per ping','Glider Record of AZFP on/off (c azfp on)','FontSize',20)

% time align
timeindex2 = zeros(1, length(azfpData(1).Date));
for ii = 1:length(azfpData(1).Date)
    [~,timeindex2(ii)] = min(abs(azfpData(1).Date(ii) - azfp_time));
end
azfpData(1).azfp_on=azfp_on(timeindex2);
azfpidx=azfp_on(timeindex2);

% remove the pings where the glider says the AZFP should not be on
azfpData(1).Date = azfpData(1).Date(azfpidx==0);
azfpData(1).BatteryMain = azfpData(1).BatteryMain(azfpidx==0);
azfpData(1).BatteryTx = azfpData(1).BatteryTx(azfpidx==0);
azfpData(1).Depth = azfpData(1).Depth(azfpidx==0);

azfpData(1).N = azfpData(1).N(azfpidx==0,:);
azfpData(2).N = azfpData(2).N(azfpidx==0,:);
azfpData(3).N = azfpData(3).N(azfpidx==0,:);
azfpData(4).N = azfpData(4).N(azfpidx==0,:);

azfpData(1).Sv = azfpData(1).Sv(azfpidx==0,:);
azfpData(2).Sv = azfpData(2).Sv(azfpidx==0,:);
azfpData(3).Sv = azfpData(3).Sv(azfpidx==0,:);
azfpData(4).Sv = azfpData(4).Sv(azfpidx==0,:);

azfpData(1).TS = azfpData(1).TS(azfpidx==0,:);
azfpData(2).TS = azfpData(2).TS(azfpidx==0,:);
azfpData(3).TS = azfpData(3).TS(azfpidx==0,:);
azfpData(4).TS = azfpData(4).TS(azfpidx==0,:);

output = azfpData;
end