function [ calibrationConstants, useMixerOverSqrtDiode ] = loadFONTCalibrationConstants( filePath )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    calFile = fopen(filePath);
    importedData = textscan(calFile,'%s%s%s%s%s');
    fclose(calFile);
    
    calAmplitudes = importedData{2}(1:3);
    calFrequency = importedData{3}(1:3);
    calPhaseShift = importedData{4}(1:3);
    calOffset = importedData{5}(1:3);
    try
        useMixerOverSqrtDiode = str2double(importedData{2}(4));
    catch
        disp('Searching for mixerOverSqrtDiode in .mat file');
        calMatFilePath = [strrep(filePath,'Constants','') '.mat'];
        calMatFilePath = strrep(calMatFilePath,'Processed','Extracted');
        try
            calMatFile = load(calMatFilePath);
            useMixerOverSqrtDiode = calMatFile.useMixerOverSqrtDiode;
            
            calFile = fopen(filePath,'a');
            fwrite(calFile,sprintf('useMixerOverSqrtDiode, %d',useMixerOverSqrtDiode));
            fclose(calFile);
            
        catch
            disp('Could not find useMixerOverSqrtDiode in %s. Set to NaN.',calMatFile);
            useMixerOverSqrtDiode = NaN;
        end
    end
    
    calibrationConstants = NaN*ones(3,4);
    for i=1:3
        calibrationConstants(i,1) = str2double(calAmplitudes{i});
        calibrationConstants(i,2) = str2double(calFrequency{i});
        calibrationConstants(i,3) = str2double(calPhaseShift{i});
        calibrationConstants(i,4) = str2double(calOffset{i});
    end
end

