function [ processedData ] = processData( dataDir, frascatiCalFile, saveIfBPMDataBad, saveIfPhaseDataBad)
% old: function [ processedData ] = processData( dataDir, frascatiCalFile, sampleRanges )
%processData Takes data directory and creates an output struct
%processedData with phase information and statistics.
%   sampleRanges. Struct defining sample ranges to use when
%   calculating statistics (sampleRanges.frascatiStartSamples,
%   sampleRanges.frascatiEndSamples, sampleRanges.petsStartSample,
%   sampleRages.petsEndSample etc.). If sampleRanges is not given, user is
%   asked to manually select them.
%   saveIfBPMDataBad - Will save phaseData to processedData even if BPM
%   data is bad (optional, default = false).
%   saveIfPhaseDataBad - Will save BPM Data to processedData even if Phase
%   data is bad (optional, default = false).
    if nargin < 3
        saveIfBPMDataBad = false;
    end
    if nargin<4
        saveIfPhaseDataBad = false;
    end

    % Get calibration constants
    [frascatiCalibrationFactors,useMixerOverSqrtDiode] = loadFrascatiCalibrationConstants(frascatiCalFile);
    frascatiCalStruct = struct();
    frascatiCalStruct.calibrationFactors = frascatiCalibrationFactors;
    frascatiCalStruct.useMixerOverSqrtDiode = useMixerOverSqrtDiode;
    
    % Get file list
    fileList = dir([dataDir '/*.mat']);
    nFiles = length(fileList);
    
    processedData = struct();
    processedData.phaseData = struct();
    processedData.bpmData = struct();
    
    processedData.dataDir = dataDir;
    
    % split up the information in the name of the data directory (assuming
    % names with format data_time_description)
    [~,tmpName,~] = fileparts(dataDir);
    processedData.dataName = tmpName;
    splitName = regexp(tmpName,'_','split');
    if (length(splitName) >= 2)
        processedData.dataDate = splitName{1};
        processedData.dataTime = splitName{2};
        
        if (length(splitName) > 2)
            processedData.dataDescription = '';
            for i=3:length(splitName)
                if i<length(splitName)
                    processedData.dataDescription = [processedData.dataDescription splitName{i} '_'];
                else
                    processedData.dataDescription = [processedData.dataDescription splitName{i}];
                end  
            end
        end
    end 
    processedData.phaseData.frascatiCalibrations.fileName = frascatiCalFile;
    processedData.phaseData.frascatiCalibrations.calibrationFactors = frascatiCalibrationFactors;
    processedData.phaseData.frascatiCalibrations.useMixerOverSqrtDiode = useMixerOverSqrtDiode;
    
    
    
    % Save init settings to struct
    initSettingsFileList = dir([dataDir '/initSettings_*.dat']);
    processedData.initSettings = struct();
    for i=1:length(initSettingsFileList)
        tmpInitSettings = loadInitSettings([dataDir '/' initSettingsFileList(i).name]);
        
        dateStamp = initSettingsFileList(i).name;
        dateStamp = strrep(dateStamp,'.dat','');
        
        eval(['processedData.initSettings.' dateStamp ' = tmpInitSettings;']);
    end
     
    
    phaseArraysAreInit = false;
    bpmArraysAreInit = false;
    otherArraysAreInit = false;
    nGoodFiles = 0;
    fprintf(1,'Attempting to load %d files in %s...\n',nFiles, dataDir);
    for i=1:length(fileList)
        fprintf(1,'File %d of %d: %s\n',i,length(fileList),fileList(i).name);
        load([dataDir '/' fileList(i).name]); % loads tmpDataStruct
        
        % Extract phase data
        [tmpPhaseData,phaseDataIsGood] = extractPhaseData(tmpDataStruct,frascatiCalStruct);
        
        % Extract BPM data
        [tmpBPMData,bpmDataIsGood] = extractBPMData(tmpDataStruct);

%         fprintf(1,'%d %d\n',phaseDataIsGood,bpmDataIsGood);
        
        % INITIALISE SOME STRUCTS/ARRAYS
        if (~otherArraysAreInit)
            processedData.dataSetComments = cell(1,nFiles);
            otherArraysAreInit = true;
        end
        
        if (~phaseArraysAreInit && phaseDataIsGood)
            frascatiNMonitors = tmpPhaseData.frascatiNMonitors;
            frascatiNSamples = tmpPhaseData.frascatiNSamples;
            frascatiTimePerSample =  tmpPhaseData.frascatiTimePerSample;
            petsNSamples = tmpPhaseData.petsNSamples;
            petsTimePerSample = tmpPhaseData.petsTimePerSample;
            
            processedData.phaseData.frascatiMixers = NaN*ones(frascatiNMonitors,nFiles,frascatiNSamples);
            processedData.phaseData.frascatiDiodes = NaN*ones(frascatiNMonitors,nFiles,frascatiNSamples);
            processedData.phaseData.frascatiPhases = NaN*ones(frascatiNMonitors,nFiles,frascatiNSamples);
            processedData.phaseData.petsI = NaN*ones(nFiles,petsNSamples);
            processedData.phaseData.petsQ = NaN*ones(nFiles,petsNSamples);
            processedData.phaseData.petsPower = NaN*ones(nFiles,petsNSamples);
            processedData.phaseData.petsPhase = NaN*ones(nFiles,petsNSamples);
            processedData.phaseData.frascatiNSamples = frascatiNSamples;
            processedData.phaseData.frascatiNMonitors = frascatiNMonitors;
            processedData.phaseData.petsNSamples = petsNSamples;
            processedData.phaseData.frascatiTimePerSample = frascatiTimePerSample;
            processedData.phaseData.petsTimePerSample = petsTimePerSample;
                        
            phaseArraysAreInit = true;
        end
        
        if (~bpmArraysAreInit && bpmDataIsGood) % initialise bpm arrays
            nBPMs = tmpBPMData.nBPMs;
            nProps=  tmpBPMData.nProps;
            bpmNames = tmpBPMData.bpmNames;
            bpmProps = tmpBPMData.bpmProps;
            
            processedData.bpmData.nBPMs = nBPMs;
            processedData.bpmData.nProps = nProps;
            processedData.bpmData.bpmNames = bpmNames;
            processedData.bpmData.bpmProps = bpmProps;
            
            for bpm=1:nBPMs    
                for prop=1:nProps
                    bpmNSamples = eval(['length(tmpBPMData.' bpmNames{bpm} bpmProps{prop} ');']);
                    
                    eval(['processedData.bpmData.' bpmNames{bpm} bpmProps{prop} '= NaN*ones(nFiles,bpmNSamples);']);
                end  
            end
            
            bpmArraysAreInit = true;
        end

        if ((phaseDataIsGood || saveIfPhaseDataBad) && (bpmDataIsGood || saveIfBPMDataBad)) % save the new data if it's good
            nGoodFiles = nGoodFiles + 1;
            
            processedData.dataSetComments{nGoodFiles} = tmpDataStruct.comment;
            
            if (phaseDataIsGood)
                processedData.phaseData.frascatiMixers(:,nGoodFiles,:) = tmpPhaseData.frascatiMixers;
                processedData.phaseData.frascatiDiodes(:,nGoodFiles,:) = tmpPhaseData.frascatiDiodes;
                processedData.phaseData.frascatiPhases(:,nGoodFiles,:) = tmpPhaseData.frascatiPhases;
                processedData.phaseData.petsI(nGoodFiles,:) = tmpPhaseData.petsI;
                processedData.phaseData.petsQ(nGoodFiles,:) = tmpPhaseData.petsQ;
                processedData.phaseData.petsPower(nGoodFiles,:) = tmpPhaseData.petsPower;
                processedData.phaseData.petsPhase(nGoodFiles,:) = tmpPhaseData.petsPhase;
            end

            if (bpmDataIsGood)
                for bpm=1:nBPMs    
                    for prop=1:nProps
                        eval(['processedData.bpmData.' bpmNames{bpm} bpmProps{prop} '(nGoodFiles,:) = tmpBPMData.' bpmNames{bpm} bpmProps{prop} ';']);
                        %e.g. processedData.bpmData.CC.SVBPM0235.Samples.samples(nGoodFiles,:) = tmpBPMData.CC.SVBPM0235.Samples.samples;
                    end  
                end
            end

        else
            fprintf(1,'Skipped file %s (bad data)\n', fileList(i).name);
            % move file to another directory
            % (badData/badFile/)
        end
    
  
    end
    fprintf(1,'%d out of %d files were loaded.\n', nGoodFiles, nFiles);
    
    % Strip NaNs off the end of arrays (made them as long as number of
    % files but some files may have been skipped)
    if nGoodFiles == 0
        processedData = [];
        return; 
    elseif nGoodFiles < nFiles
        if (exist('processedData.phaseData.frascatiMixers','var'))
            processedData.phaseData.frascatiMixers(:,nGoodFiles+1:nFiles,:) = [];
            processedData.phaseData.frascatiDiodes(:,nGoodFiles+1:nFiles,:) = [];
            processedData.phaseData.frascatiPhases(:,nGoodFiles+1:nFiles,:) = [];
            processedData.phaseData.petsI(nGoodFiles+1:nFiles,:) = [];
            processedData.phaseData.petsQ(nGoodFiles+1:nFiles,:) = [];
            processedData.phaseData.petsPower(nGoodFiles+1:nFiles,:) = [];
            processedData.phaseData.petsPhase(nGoodFiles+1:nFiles,:) = [];
        end

        processedData.dataSetComments(nGoodFiles+1:nFiles) = [];
        
        for bpm=1:nBPMs    
            for prop=1:nProps
                if (exist(sprintf('processedData.bpmData.%s.%s',bpmNames{bpm},bpmProps{prop}),'var'))
                    eval(['processedData.bpmData.' bpmNames{bpm} bpmProps{prop} '(nGoodFiles+1:nFiles,:) = [];']);
                end
            end  
        end
    end
    
    % Align pulses (see Piotr getAligned2)
    
    % Delete pulses where we don't have both bpm and phase data
    % move to badData/overlap?
    % Check BPM vs. phase mon timings (FOR NOW JUST FIXED AT 1 PULSE AS
    % CHECKED IN ONE DATA SET 
%     oasisDelay = 1;
%     removeIndOasis = (nGoodFiles-oasisDelay+1):nGoodFiles;
%     removeIndXeneric = 1:oasisDelay;
%     
%     processedData.phaseData.frascatiMixers(:, removeIndOasis, :) = [];
%     processedData.phaseData.frascatiDiodes(:, removeIndOasis, :) = [];
%     processedData.phaseData.frascatiPhases(:, removeIndOasis, :) = [];
%     processedData.phaseData.petsI(removeIndOasis, :) = [];
%     processedData.phaseData.petsQ(removeIndOasis, :) = [];
%     processedData.phaseData.petsPower(removeIndOasis, :) = [];
%     processedData.phaseData.petsPhase(removeIndOasis, :) = [];
% 
%     processedData.dataSetComments{removeIndXeneric} = [];
% 
%     for bpm=1:nBPMs    
%         for prop=1:nProps
%             eval(['processedData.bpmData.' bpmNames{bpm} bpmProps{prop} '(removeIndXeneric,:) = [];']);
%         end  
%     end

    
    
%     % Start and end samples
%     if nargin<3
%         sampleRanges = selectDataSampleRanges(processedData.phaseData);
%     end
%     processedData.phaseData.frascatiStartSamples = sampleRanges.frascatiStartSamples;
%     processedData.phaseData.frascatiEndSamples = sampleRanges.frascatiEndSamples;
%     processedData.phaseData.petsStartSample = sampleRanges.petsStartSample;
%     processedData.phaseData.petsEndSample = sampleRanges.petsEndSample;
%     
%     % Calculate time axes 
%     tmpFrascatiTimeAxis = linspace(0, (frascatiTimePerSample.*frascatiNSamples)-frascatiTimePerSample, frascatiNSamples);
%     processedData.phaseData.frascatiTimeAxes = {};
%     for mon=1:frascatiNMonitors
%         processedData.phaseData.frascatiTimeAxes{mon} = tmpFrascatiTimeAxis - tmpFrascatiTimeAxis(sampleRanges.frascatiStartSamples(mon));
%     end
% 
%     tmpPetsTimeAxis = linspace(0, (petsTimePerSample.*petsNSamples)-petsTimePerSample, petsNSamples);
%     processedData.phaseData.petsTimeAxis = tmpPetsTimeAxis - tmpPetsTimeAxis(sampleRanges.petsStartSample);
%     
    
    
    % Is good pulse flag   
    % move pulses with bad transmission to another directory
    % (badData/badTransmission/)
%     for i=1:nGoodFiles
%         
%         tmpTransmission = processedData.bpmData.CT.SVBPM0155S.Samples.samples(i,:);
%         
%         
%         
%     end
     nGoodPulses = nGoodFiles;
    
    
    % Subtract init phase
%    
%     subtractInitNPulses = ceil(nGoodPulses./10); % subtract mean of first 10% of pulses
%     
%     tmpMeanFrascatiPhases = zeros(1,frascatiNMonitors);
%     for mon=1:frascatiNMonitors
%         tmpFrascatiPhase = squeeze(processedData.phaseData.frascatiPhases(mon,1:subtractInitNPulses,:));
%         tmpFrascatiPhase = nanmean(nanmean(tmpFrascatiPhase(:,sampleRanges.frascatiStartSamples(mon):sampleRanges.frascatiEndSamples(mon))));
%         tmpMeanFrascatiPhases(mon) = tmpFrascatiPhase;
%     end
% 
%     tmpPetsPhase = processedData.phaseData.petsPhase(1:subtractInitNPulses,:);
%     tmpPetsMeanPhase = nanmean(nanmean(tmpPetsPhase(:,sampleRanges.petsStartSample:sampleRanges.petsEndSample)));
% 
%     for mon=1:frascatiNMonitors
%         processedData.phaseData.frascatiPhases(mon,:,:) = processedData.phaseData.frascatiPhases(mon,:,:) - tmpMeanFrascatiPhases(mon);
%     end
%     processedData.phaseData.petsPhase = processedData.phaseData.petsPhase - tmpPetsMeanPhase;



    
    
    
    %Calculate phase statistics
%     processedData.phaseStatistics = calculatePhaseStatistics(processedData.phaseData);
end

