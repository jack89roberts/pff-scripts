function [ calibrationConstants, bprNames ] = loadBPRCalibrationConstants( filePath )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    calFile = fopen(filePath);
    importedData = textscan(calFile,'%s%s%s%s%s');
    fclose(calFile);
    
    bprNames = strrep(importedData{1},',','');
    calAmplitudes = importedData{2};
    calFrequency = importedData{3};
    calPhaseShift = importedData{4};
    calOffset = importedData{5};
    
    calibrationConstants = NaN*ones(length(calAmplitudes),4);
    for i=1:length(calAmplitudes)
        calibrationConstants(i,1) = str2double(calAmplitudes{i});
        calibrationConstants(i,2) = str2double(calFrequency{i});
        calibrationConstants(i,3) = str2double(calPhaseShift{i});
        calibrationConstants(i,4) = str2double(calOffset{i});
    end
end

