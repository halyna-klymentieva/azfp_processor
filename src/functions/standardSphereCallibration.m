function output = standardSphereCallibration(calibrationOffsets,azfpData)
%STANDARDSPHERECALLIBRATION Standard Sphere Calibration application

% Echosounder 59016 (Davies) calibration offset from standard sphere
% calibration conducted June 2024.  Davies and Mesquita have
% calibration files.
arguments (Input)
    calibrationOffsets
    azfpData
end

arguments (Output)
    output
end


fprintf('Applying standard sphere calibration...\n');
tic
azfpData(1).Sv(:, :) = azfpData(1).Sv(:, :) + calibrationOffsets(1);
azfpData(2).Sv(:, :) = azfpData(2).Sv(:, :) + calibrationOffsets(2);
azfpData(3).Sv(:, :) = azfpData(3).Sv(:, :) + calibrationOffsets(3);
azfpData(4).Sv(:, :) = azfpData(4).Sv(:, :) + calibrationOffsets(4);
toc

output = azfpData;
end