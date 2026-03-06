%Version 1 created in 2020 by Halyna Klymentieva

%% Prepare your workspace and file directories
addpath(genpath(pwd));
clc
clear variables
close all

%% Merge data
config.saveWithDateName = 0;
config.outputFolder = fullfile(pwd, '..', 'output');
% merging all output dives files
filenames = getDivesFilenames(config.outputFolder);
[Dives, bottomDepth] = mergeDives(filenames);
saveDiveData(fullfile(config.outputFolder, '0-merged-dives.mat'), Dives, bottomDepth)

%% Merge data
drawAndSaveFigures(Dives, bottomDepth, config)
