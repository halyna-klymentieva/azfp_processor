%% Load in glider data

% load up the raw datafile downloaded from CEOTR
load('/Users/Andrea/Documents/AMesquita2025/UNB/dataAnalysis/preyData/matlab/missions/shad_20230601_168_delayed_43c7_19bc_dd9c.mat');

% put that data into the variable gliderdata
gliderdata1 = shad_20230601_168_delayed;
clear shad_20230601_168_delayed

%% Load in glider data

% load up the raw datafile downloaded from CEOTR
load('/Users/Andrea/Documents/AMesquita2025/UNB/dataAnalysis/preyData/matlab/missions/shad_20230601_168_delayed_4392_d8fb_7b88.mat');

% put that data into the variable gliderdata
gliderdata2 = shad_20230601_168_delayed;
clear shad_20230601_168_delayed

%% Load in glider data

% load up the raw datafile downloaded from CEOTR
load('/Users/Andrea/Documents/AMesquita2025/UNB/dataAnalysis/preyData/matlab/missions/shad_20230601_168_delayed_bfa7_8c28_ee8f.mat');

% put that data into the variable gliderdata
gliderdata3 = shad_20230601_168_delayed;
clear shad_20230601_168_delayed

%% Remove variables and bind structures

%variables to keep
keywords = ["time", "depth", "lat", "lon", "pressure", "temperature", ...
            "salinity", "density", "echosounder"];

%filter fields for gliderdata1
fn = fieldnames(gliderdata1);
keep = false(size(fn));
for k = 1:numel(fn)
    for kw = 1:numel(keywords)
        if contains(fn{k}, keywords(kw), 'IgnoreCase', true)
            keep(k) = true;
            break
        end
    end
end
gliderdata1 = rmfield(gliderdata1, fn(~keep));

%filter fields for gliderdata2
fn = fieldnames(gliderdata2);
keep = false(size(fn));
for k = 1:numel(fn)
    for kw = 1:numel(keywords)
        if contains(fn{k}, keywords(kw), 'IgnoreCase', true)
            keep(k) = true;
            break
        end
    end
end
gliderdata2 = rmfield(gliderdata2, fn(~keep));

%filter fields for gliderdata3
fn = fieldnames(gliderdata3);
keep = false(size(fn));
for k = 1:numel(fn)
    for kw = 1:numel(keywords)
        if contains(fn{k}, keywords(kw), 'IgnoreCase', true)
            keep(k) = true;
            break
        end
    end
end
gliderdata3 = rmfield(gliderdata3, fn(~keep));

%merge
merged = gliderdata1;
fields = fieldnames(merged);

for i = 1:numel(fields)
    merged.(fields{i}) = [gliderdata1.(fields{i}); ...
                          gliderdata2.(fields{i}); ...
                          gliderdata3.(fields{i})];
end

%% Save file
shad_20230601_168_delayed = merged;
save('/Users/Andrea/Documents/AMesquita2025/UNB/dataAnalysis/preyData/matlab/missions/shad_20230601_168_delayed_merged.mat', 'shad_20230601_168_delayed');