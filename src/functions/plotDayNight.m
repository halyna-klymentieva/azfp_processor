function Dives = plotDayNight(Dives)
%PLOTDAYNIGHT Summary of this function goes here
%Plot day and night profiles by depth

for i = 1:length(Dives)
    dt = datetime(Dives(i).starttime, 'ConvertFrom', 'datenum'); %convert time
    Dives(i).localtime = dt - hours(3); %convert to local time
end

for i = 1:length(Dives)
    hr = hour(Dives(i).localtime);

    if hr >= 7 && hr <= 20
        Dives(i).tod = "Day";
    elseif hr >= 23 || hr <= 4
        Dives(i).tod = "Night";
    else
        Dives(i).tod = "Other";
    end
end
end