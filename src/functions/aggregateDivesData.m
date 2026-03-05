function [Dives, bottomDepth] = aggregateDivesData(Output, config)
%AGGREGATEDIVESDATA Aggregate dives data

%% Average 10 cm vertical resolution into 1 m depth bins to make the matrices smaller for better storage space
% Step 6: Organized Dive Structure (Stable Version)- Halyna's version
[azfpAggregateData, Dives] = aggregateVerticalResolution(Output, config.maxDepth);

%% Moving average to determine noise floor for each bin (Steps 3 and 4) and noise removal
Dives = filterNoiseFloor(config.maxDepth, Output, azfpAggregateData, Dives);
%% AZFP_Unmasked_Masked_Comparison routine (D. Mossman) - lines 626–750
Dives = unmaskedMaskedComparison(Dives);

%% Calculate numerical density
Dives = calcNumericalDensity(Dives);

%% Plot day/night
Dives = plotDayNight(Dives);

%% Calculate bottom depth and end dives
bottomDepth = calculateDiveBottomDepth(Output, Dives);

%% Use this code to save Dives data to file
saveDiveData(config.divesDataCachePath, Dives, bottomDepth)
end