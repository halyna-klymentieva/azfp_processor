function output = getAZFPProcessResult(savePath, sourceFolder)
%GETAZFPPROCESSRESULT Summary of this function goes here
% AZFP Code from ASL for converting from engineering to real units
% Parameter description and the default value (if the value is omitted):
% >> ParametersAZFP;[Output,Par] = ProcessAZFP(Parameters);
%
% Ver 1.3 September 2017
% written by Dave Billenness
% ASL Environmental Sciences Inc.
% 1-6703 Rajpur Place, Victoria, B.C., V8M 1Z5, Canada
% T: +1 (250) 656-0177 ext. 126
% E: dbillenness@aslenv.com
% w: http://www.aslenv.com/
% For any suggestions, comments, questions or collaboration, please contact me.

arguments (Input)
    savePath string
    sourceFolder string = ""
end

arguments (Output)
    output
end

% If savePath already exists data is loaded from file
if isfile(savePath)
    tic
    fprintf('Loading data from cache. Will not be re-generated if not deleted.\n');
    load(savePath, 'Output');
    toc
else
    %% Params
    % FILE LOADING AND AVERAGING:
    % Parameters.ProcDir = 0; 1 will prompt for an entire directory to
    % process, = 0 will prompt to load individual files in a directory
    Parameters.ProcDir = 1;

    % FILE LOADING AND AVERAGING:
    % Parameters.ProcDir = 0; 1 will prompt for an entire directory to
    % process, = 0 will prompt to load individual files in a directory
    Parameters.folderName = sourceFolder;

    % Parameters.datafilename = ''; % '' will prompt for hourly AZFP
    % file(s) to load, example '16010100.01A'
    Parameters.datafilename = '';

    % Parameters.xmlfilename = ''; % prompt for XML filename if no XML file exists
    % in the directory, example '15101614.XML'
    Parameters.xmlfilename = '';

    % Parameters.Salinity = 35; % Salinity in psu
    Parameters.Salinity = 33;

    % Parameters.Bins2Avg = 10; % number of range bins to average
    Parameters.Bins2Avg = 1; % 0.45m ->15

    % Parameters.Time2Avg = 60; % number of time values to average
    Parameters.Time2Avg = 1; %1

    % Parameters.Pressure = 50; % in dbars (~ depth of instrument in meters).
    % This can be approximate and is used in the soundspeed and absorption calc
    Parameters.Pressure = 100;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOTTING:
    % Parameters.Plot = 1; % show an echogram plot for each channel
    Parameters.Plot = 1;
    % Parameters.Channel: freq to plot #1-4, default 1
    Parameters.Channel = 1;
    % Parameters.Value2Plot = 2; 1,2,3,4 = Counts, Sv, TS, Temperature/Tilts, default 2
    Parameters.Value2Plot = 2;

    % Parameters.NoiseFloor = 10000; % for Sv and Ts plotting only, values
    % with counts < NoiseFloor will be set to -150, can use individual values
    % for each frequency, ex. 'NoiseFloor',[10000; 11000; 10500; 12500]
    % Default = 10000.
    %Parameters.NoiseFloor = [11014;9966;13206;12760]; % AZFP SN 59015
    Parameters.NoiseFloor = [12564; 9660; 13604; 10780]; % AZFP SN 59016
    % from the manufacturer; 59015/59016 are serial numbers

    % Parameters.Orientation = 0 instrument on bottom looking up (range bins), 1 at surface
    % looking down (depth bins). This changes the ydir on the echogram plots only. Default is 1.
    Parameters.Orientation = 1;

    % Parameters.UseTiltCorr = 0; Use the tilt corrected ranges for the echogram plots,
    % default 0. Will give a warning if the tilt magnitudes are unreasonable (> 20 deg)
    Parameters.UseTiltCorr = 0;
    %% Load AZFP Data
    % If no sourceFolder - will pull up your file explorer on the current path;
    % first select the *folder* with the AZFP data, then select the XML file
    % within that folder. If there is only one XML file, the code will find
    % it automatically, and you will not need to select it
    tic
    [Output, ~] = ProcessAZFP(Parameters);
    toc
    %% Davies lab code: sort Output by date
    [~, order] = sort([Output(1).Date], 'ascend');

    Output(1).Date = Output(1).Date(order);
    Output(1).BatteryMain = Output(1).BatteryMain(order);
    Output(1).BatteryTx = Output(1).BatteryTx(order);
    Output(1).Depth = Output(1).Depth(order);

    Output(1).N = Output(1).N(order, :);
    Output(2).N = Output(2).N(order, :);
    Output(3).N = Output(3).N(order, :);
    Output(4).N = Output(4).N(order, :);

    Output(1).Sv = Output(1).Sv(order, :);
    Output(2).Sv = Output(2).Sv(order, :);
    Output(3).Sv = Output(3).Sv(order, :);
    Output(4).Sv = Output(4).Sv(order, :);

    Output(1).TS = Output(1).TS(order, :);
    Output(2).TS = Output(2).TS(order, :);
    Output(3).TS = Output(3).TS(order, :);
    Output(4).TS = Output(4).TS(order, :);

    % Save the Processed Data for later use
    cacheOutput(savePath, Output)
end

output = Output;
end
