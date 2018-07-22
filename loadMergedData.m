function [ CTFData, calStruct, dataDir ] = loadMergedData( dataSetName, frascatiCalTimeStamp, dataBaseDir, calDir )
%loadMergedData Loads all data from the data set dataSetName in to the
%array of structs CTFData. Also return frascati calibration information.
%   dataSetName and frascatiCalStamp are names not file paths (e.g.
%   20150421_0951_test).
%   dataDir and calDir (optional): Directory where dataSetName and frascatiCalTimeStamp are found.
%   Saves merged data for this data set if it is not already saved.

    if (nargin<3 || isempty(dataBaseDir))
        if (strcmp(getenv('USER'),'jack') == 1) % If username is 'jack' assume environment is Jack's laptop
            dataBaseDir = ['/home/jack/PhaseFeedforward/CTFData/' dataSetName(1:6)];
            dataDir =  [dataBaseDir '/' dataSetName];
            %calDir = [dataBaseDir '/FrascatiCalibrations'];
        else % otherwise assume control room environment
            %dataDir = ['/nfs/cs-ccr-nfs6/vol29/share/ctf/data/ctfmod/MatLab/PhaseFeedforward/data/' dataSetName];
            dataDir = ['/user/ctf3op/PhaseFeedforward/Data/' dataSetName];
            %calDir = '/nfs/cs-ccr-nfs6/vol29/share/ctf/data/ctfmod/MatLab/PhaseFeedforward/FrascatiCalibrations/data';
            %calDir = '/user/ctf3op/PhaseFeedforward/FrascatiCalibrations';
        end
    elseif (isempty(calDir))
        dataDir = [dataBaseDir '/' dataSetName];
        warning('If dataDir given manually, calDir should also be given. calStruct will be empty.');
        calDir = '';
    end
    
    if (exist([dataDir '/Merged'],'dir')) % If merged data already exists, load it.
        mergedName = [dataDir '/Merged/' dataSetName '.mat'];
        fprintf('Loading merged data in %s...\n',mergedName);
        load(mergedName);
    else % If not create and save the merged data.
        fprintf('Creating and saving merged data for %s...\n',dataSetName);
        CTFData = saveMergedData(dataDir);
    end
    
    if (nargin < 2 || isempty(frascatiCalTimeStamp)) % If the calibration to use wasn't given, load it from the initSettings file in the data set.
        [frascatiCalTimeStamp,~,~] = loadInitSettings([dataDir '/initSettings.dat']);
    end
    
    if (~isempty(frascatiCalTimeStamp))
        if (strcmp(getenv('USER'),'jack') == 1) % If username is 'jack' assume environment is Jack's laptop
            calDir = ['/home/jack/PhaseFeedforward/CTFData/' frascatiCalTimeStamp(1:6) '/FrascatiCalibrations'];
        else % otherwise assume control room environment
            calDir = '/user/ctf3op/PhaseFeedforward/FrascatiCalibrations';
        end
    else
        calDir = '';
    end

    
    
    calStruct = struct();
    if (~isempty(calDir) && ~isempty(frascatiCalTimeStamp))
        fprintf('Loading calibration constants from %s...\n',frascatiCalTimeStamp);
        [calibrationConstants, useMixerOverSqrtDiode] = loadFrascatiCalibrationConstants(sprintf([calDir '/frascatiCalibrationConstants_' frascatiCalTimeStamp]));

        calStruct.frascatiCalTimeStamp = frascatiCalTimeStamp;
        calStruct.calDir = calDir;
        calStruct.calibrationConstants = calibrationConstants;
        calStruct.useMixerOverSqrtDiode = useMixerOverSqrtDiode;
    end
end

