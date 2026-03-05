function filenames = getFilenamesByDate(dateOfData, sourceFolder)
%GETFILENAMESBYDATE Scans source folder for files with specific date in names

cwd = pwd;
cd(sourceFolder)

noHyphen = replace(dateOfData, "-", "");
filelist = dir(noHyphen+"*.01*.");
filelist = filelist(~endsWith({filelist.name}, {'.evi'}));
[~, I] = sort(datenum({filelist.date}));
f = char({filelist.name});
filelist = f(I, :);
cd(cwd)

filenames = filelist;
end
