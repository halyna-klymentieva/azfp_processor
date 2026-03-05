function filenames = getDivesFilenames(sourceFolder)
%GETDIVESFILENAMES Scans output directory for dive files
cwd = pwd;
cd(sourceFolder)

filelist = dir("dives-*.mat");
[~, I] = sort({filelist.name});
f = char({filelist.name});
filelist = f(I, :);
cd(cwd)

filenames = string(filelist);
end