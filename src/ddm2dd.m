function [dd] = ddm2dd(ddm)
%Degrees Decimal Minutes to Decimal Degrees Converter
%   Matlab version of Hansen's R function that converts latitudes and
%   longitudes from degrees decimal minutes to decimal degrees
dd = NaN*ones(length(ddm),1);

for i = 1:length(ddm)
    ddm = split(ddm, ' ');
    
    deg = str2num(ddm(i,1));
    dec = str2num(ddm(i,2))/60;
    
    if deg < 1
        dd(i) = deg - dec;
    else
        dd(i) = deg + dec;
    end
end
end

