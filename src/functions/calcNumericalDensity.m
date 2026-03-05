function Dives = calcNumericalDensity(Dives)
%CALCNUMERICALDENSITY Calculate dives numerical density
% calculate numerical density (ind. per m^3) from merged Dive structures
sigma_bs = [1.3e-11, 7.6e-11, 7.9e-11]; % for 200kHz, 455kHz, 769kHz

n = length(Dives);

for ii = 1:n
    Dives(ii).P(2).N = Dives(ii).P(2).masked / sigma_bs(1); % 200 kHz
    Dives(ii).P(3).N = Dives(ii).P(3).masked / sigma_bs(2); % 455 kHz
    Dives(ii).P(4).N = Dives(ii).P(4).masked / sigma_bs(3); % 769 kHz

    Dives(ii).P(2).mN = median(Dives(ii).P(2).N, 'omitnan');
    Dives(ii).P(3).mN = median(Dives(ii).P(3).N, 'omitnan');
    Dives(ii).P(4).mN = median(Dives(ii).P(4).N, 'omitnan');
end
% scale density values by volume to estimate abundance
% estimate biomass by incorporating mean organism weight
% estimate energy density by incorporating kJ per gram of lipid
end