function [gliderDepth, gliderTime, gliderData] = loadGliderData(gliderFileName, gliderVariableName)
%LOADGLIDERDATA Loads depths and times from glider raw data file
arguments (Input)
    gliderFileName
    gliderVariableName
end

arguments (Output)
    gliderDepth
    gliderTime
    gliderData
end

S = load(gliderFileName, gliderVariableName);
gliderdata = S.(gliderVariableName);
clear S

% change glider unix time format to same format as in Output file (matlab
% time format)
unix_epoch = datenum(1970, 1, 1, 0, 0, 0);
gliderdata.time = gliderdata.time ./ 86400 + unix_epoch;

% find the non-NaN indices of glider depth and get their values + the time
% at which they were recorded
nanindex = find(~isnan(gliderdata.depth));

gliderDepth = gliderdata.depth(nanindex);
gliderTime = gliderdata.time(nanindex);
gliderData = gliderdata;
end
