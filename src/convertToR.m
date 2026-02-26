%% Load processed, differenced, and masked data
% A. Mesquita
% 1 February 2024

clc
clear variables
close all
addpath( genpath('/Users/amesquit/Documents/AMesquita/UNB/dataAnalysis/processedAZFP/novCode/2022'))

date = input('Enter the numerical day of the data: ','s');

filename = strcat("/Users/amesquit/Documents/AMesquita/UNB/dataAnalysis/processedAZFP/novCode/2022/",date,"Jun_Differencing_Data.mat");
load(filename);

%% Convert dates to DateTime Format
dateTime = datetime(Output(1).Date, 'ConvertFrom', 'datenum');

%% Save variables to read in R

filename = strcat("/Users/amesquit/Documents/AMesquita/UNB/dataAnalysis/processedAZFP/novCode/2022/",date,"Aug_Differencing_Data_copy.mat");
fields = {'Tx','Ty','T','filename','HourlyAvgTemp','SoundSpeed','N','Range','TiltCorrRange','Sv','TS','seaAbs','Freq','Bins2Avg','Time2Avg','BurstInt','PingPerProfile','NumAcqPings','DataType'};
Output = rmfield(Output, fields);
Dive = rmfield(Dive, 'P');
save(filename, 'Output', 'P', 'Dive', 'Diff', 'DiveDiff', 'PDiff', '-V6');

clear filename;

%% Create csv to export for further analysis in R
% Gina L. Lonati
% 15 February 2023

% concatenate all dives by frequency
PSum = struct('P1', [], 'P2', [], 'P3', [], 'P4', []); % create field for each frequency

for cc = 1:length(Dive)
    PSum.P1 = vertcat(PSum.P1, P(1).avg_sv(Dive(cc).Index(1):Dive(cc).Index(2), :));
    PSum.P2 = vertcat(PSum.P2, P(2).avg_sv(Dive(cc).Index(1):Dive(cc).Index(2), :));
    PSum.P3 = vertcat(PSum.P3, P(3).avg_sv(Dive(cc).Index(1):Dive(cc).Index(2), :));
    PSum.P4 = vertcat(PSum.P4, P(4).avg_sv(Dive(cc).Index(1):Dive(cc).Index(2), :));
end

export = struct('msv', [], 'freq', [], 'profile', [], 'deployment', [], 'depth', []);

l = 1; % to iterate through multiple dives/profiles

for i = 1:length(Dive) % for the number of dives/profiles
    for jj = 1:4 % frequency index
        export.msv = horzcat(export.msv, Dive(i).P(jj).msv);
    end
    export.freq = repmat([repmat({'130 kHz'}, 1, width(PSum.P1)) repmat({'200 kHz'}, 1, width(PSum.P2)) repmat({'455 kHz'}, 1, width(PSum.P3)) repmat({'769 kHz'}, 1, width(PSum.P4))], 1, i); % assumes each dive/profile is same length
    export.profile(l:length(export.msv)) = repelem(i, length(export.msv)-(l-1));
    export.depth(l:length(export.msv)) = repmat(6:1:width(PSum.P1)+5, 1, 4); % start at 5, because first 5 m trimmed; assumes all frequencies same length
    export.deployment = repmat({date}, 1, length(export.msv));
    l = 1 + length(export.msv);
end

export = [num2cell(export.msv'), export.freq', num2cell(export.profile'), export.deployment', num2cell(export.depth')];

export = cell2table(export, 'VariableNames', {'med_sv', 'freq', 'profile', 'deployment', 'depth'});

% write table
writetable(export, strcat("/Users/amesquit/Documents/AMesquita/UNB/dataAnalysis/processedAZFP/novCode/2022/csv/",date,"Jun_MedianSvByProfile.txt"))