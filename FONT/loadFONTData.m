function [ FONTData, calStruct, dataDir ] = loadFONTData( dataSetName, calStamp, dataDir, calDir )
%loadFONTData Loads all data from the data set dataSetName in to the
%array of structs FONTData. Also return calibration information for most relevant calibration.
%   dataSetName and calStamp are names not file paths (e.g.
%   20150421_0951_test).
%   dataDir and calDir (optional): Directory where dataSetName and calTimeStamp are found.
%   Saves extracted data for this data set if it is not already saved.

    if (nargin<3 || isempty(dataDir))
            dataDir = ['/home/jack/PhaseFeedforward/FONTData/' dataSetName(1:6)];
            calDir = [dataDir '/Calibration'];
    elseif (isempty(calDir))
        warning('If dataDir given manually, calDir should also be given. calStruct will be empty.');
        calDir = '';
    end
    
    extractedDir = [dataDir '/Extracted'];
    extractedName = [extractedDir '/' dataSetName '.mat'];
    if (~exist(extractedDir,'dir'))
        mkdir(extractedDir);
    end
    
    if (exist(extractedName,'file')) % If extracted data already exists, load it.
        fprintf('Loading data in %s...\n',extractedName);
        FONTData = loadExtractedFONTData(extractedName);
        
    else % If not create and save the merged data.
        fprintf('Creating and saving extracted data for %s...\n',dataSetName);
        FONTData = saveExtractedFONTData(dataSetName,dataDir);
    end
    
    if (nargin < 2 || isempty(calStamp)) % If the calibration to use wasn't given, load it from the initSettings file in the data set.
        [calStamp,calDir] = getBestCalStamp(dataSetName);  
    end
    
    
    calStruct = struct();
    if (~isempty(calDir) && ~isempty(calStamp))
        fprintf('Loading calibration constants from %s...\n',calStamp);
        calPath = [calDir '/' calStamp '/Processed/frascatiCalibrationConstants_' calStamp];
        [calibrationConstants, useMixerOverSqrtDiode] = loadFONTCalibrationConstants(calPath);

        calStruct.frascatiCalTimeStamp = calStamp;
        calStruct.calDir = calDir;
        calStruct.calibrationConstants = calibrationConstants;
        calStruct.useMixerOverSqrtDiode = useMixerOverSqrtDiode;
    end


end

