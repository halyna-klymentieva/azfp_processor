function [Dives, bottomDepth] = mergeDives(filenames)
%MERGEDIVES Merges a list of dives from files

nFiles = numel(filenames);

fprintf('Merging all %d dives...\n', nFiles)
tic

DivesData = cell(1, nFiles);
DepthData = cell(1, nFiles);
for k = 1:nFiles
    fullName = fullfile(pwd, '..', 'output', filenames(k,:));
    S = load(fullName);
    DivesData{k} = S.Dives;
    DepthData{k} = S.bottomDepth;
end

tbls = cell(size(DivesData));
for k = 1:numel(DivesData)
    s = DivesData{k};
    tbls{k} = struct2table(s);
end
merge_t = vertcat(tbls{:});
Dives = table2struct( merge_t );
bottomDepth = vertcat(DepthData{:});
toc
end