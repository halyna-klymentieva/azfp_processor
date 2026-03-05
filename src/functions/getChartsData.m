function chartData = getChartsData(Dives)
%GETCHARTSDATA make matrices of each Sv variable

%% make matrices of each Sv variable
chartData = struct();
for i = 1:length(Dives)
    chartData.Sv_130(i, :) = real(10*log10(Dives(i).P(1).msv));
    chartData.Sv_200(i, :) = real(10*log10(Dives(i).P(2).msv));
    chartData.Sv_455(i, :) = real(10*log10(Dives(i).P(3).msv));
    chartData.Sv_769(i, :) = real(10*log10(Dives(i).P(4).msv));
    chartData.Sv_200_130(i, :) = Dives(i).P(1).mDiff;
    chartData.Sv_455_200(i, :) = Dives(i).P(2).mDiff;
    chartData.Sv_769_400(i, :) = Dives(i).P(3).mDiff;
    chartData.masked_130(i, :) = real(10*log10(Dives(i).P(1).mMasked));
    chartData.masked_200(i, :) = real(10*log10(Dives(i).P(2).mMasked));
    chartData.masked_455(i, :) = real(10*log10(Dives(i).P(3).mMasked));
    chartData.masked_769(i, :) = real(10*log10(Dives(i).P(4).mMasked));
    chartData.StartDiveTime(i, :) = Dives(i).starttime;
    chartData.EndDiveTime(i, :) = Dives(i).endtime;
end
%% mask echoes below the seafloor using 769 kHz matrix (or masked 455 kHz) - ask Andrea why she disabled lines 883-910

idx = NaN(size(chartData.masked_455, 1), 1);

for i = 1:size(chartData.masked_455, 1)
    rowData = chartData.masked_455(i, :);

    for j = size(rowData, 2):-1:1
        if ~isnan(rowData(j))
            idx(i) = j;
            break;
        end
    end
end

% Mask values below the seafloor (after the last valid bin)
for i = 1:size(chartData.Sv_130, 1)
    colStart = idx(i) + 1;
    if ~isnan(colStart) && colStart <= size(chartData.Sv_130, 2)
        chartData.Sv_130(i, colStart:end) = NaN;
        chartData.Sv_200(i, colStart:end) = NaN;
        chartData.Sv_455(i, colStart:end) = NaN;
        chartData.Sv_200_130(i, colStart:end) = NaN;
        chartData.Sv_455_200(i, colStart:end) = NaN;
        chartData.masked_130(i, colStart:end) = NaN;
        chartData.masked_200(i, colStart:end) = NaN;
        chartData.masked_455(i, colStart:end) = NaN;

    end
end
end
