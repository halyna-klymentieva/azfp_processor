%Version 1 created in 2020 by Halyna Klymentieva

% sv - "volume backscattering strength", intencity of sound reflected
%  by biological organism within the volume of water

%% Prepare your workspace and file directories
addpath(genpath(pwd));
clc
clear variables
close all

%% Merge data
config.saveWithDateName = 0;
config.outputFolder = fullfile(pwd, '..', 'output');
config.maxDepth = 105;
% merging all output dives files
filenames = getDivesFilenames(config.outputFolder);
[Dives, bottomDepth] = mergeDives(filenames);
saveDiveData(fullfile(config.outputFolder, '0-merged-dives.mat'), Dives, bottomDepth)

%% Merge data
drawAndSaveFigures(Dives, bottomDepth, config)
