function cam = getCameraDives(Dives, config)
%GETCAMERADIVES True for dives recorded while the shadowgraph camera was on
%   Each dive is matched, by its end time, to the last glider descent that
%   started before it in config.scheduleFile (instrument_schedule_by_cast.csv
%   from CTD_processor), and takes that descent's sg_on value.
%   Returns [] (no marks) if the schedule file is not set or not found.
cam = [];
if ~isfield(config, 'scheduleFile') || isempty(config.scheduleFile) || ~isfile(config.scheduleFile)
    warning('Camera schedule file not found - camera dives are not marked on the figures.');
    return
end
% columns: profile_id, start, end, ..., direction (9), ..., sg_on (14)
% ('end' is a MATLAB keyword, so columns are taken by position)
opts = detectImportOptions(config.scheduleFile);
vn = opts.VariableNames;
opts = setvartype(opts, vn([2 9]), 'char');
sch = readtable(config.scheduleFile, opts);
isDown = strcmp(sch.(vn{9}), 'down');
tStart = datenum(datetime(sch.(vn{2})(isDown), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z'''));
sg = sch.(vn{14})(isDown);
if iscell(sg), sg = strcmpi(sg, 'TRUE'); else, sg = logical(sg); end
[tStart, o] = sort(tStart);
sg = sg(o);

cam = false(1, numel(Dives));
for i = 1:numel(Dives)
    te = Dives(i).endtime;
    if isempty(te) || isnan(te), continue; end
    j = find(tStart <= te, 1, 'last');
    if ~isempty(j), cam(i) = sg(j); end
end
fprintf('Camera on for %d of %d dives.\n', sum(cam), numel(cam));
end
