function output = standardSphereCallibration(calibrationOffsets,azfpData)
%STANDARDSPHERECALLIBRATION Standard Sphere Calibration application

% Echosounder 59016 (Davies) calibration offsets from the post-deployment
% standard sphere calibration, Aquatron tower tank, 28-29 Aug 2025
% (Mesquita & Davies, 59016CalReport2025.pdf, Table 1).
% Offset ~ observed - theoretical TS, incl. ASL term (negative = AZFP reads too low).
% Corrected Sv = observed Sv - offset (e.g. Sv - (-4.60) = Sv + 4.60).
arguments (Input)
    calibrationOffsets
    azfpData
end

arguments (Output)
    output
end


fprintf('Applying standard sphere calibration...\n');
tic
azfpData(1).Sv(:, :) = azfpData(1).Sv(:, :) - calibrationOffsets(1);
azfpData(2).Sv(:, :) = azfpData(2).Sv(:, :) - calibrationOffsets(2);
azfpData(3).Sv(:, :) = azfpData(3).Sv(:, :) - calibrationOffsets(3);
azfpData(4).Sv(:, :) = azfpData(4).Sv(:, :) - calibrationOffsets(4);
toc

output = azfpData;
end