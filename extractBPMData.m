function [ bpmData, dataIsGood ] = extractBPMData( dataStruct )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    % Have removed:
    % Initial phase subtraction
    % Time axis calculation (based on pulse start sample)
    try
        
        % names of bpms to extract
        devsPath = 'devices/';%'/nfs/cs-ccr-nfs6/vol29/share/ctf/data/ctfmod/MatLab/PhaseFeedforward/devices/';
        bpmDevFiles = {'bpmCL.devs', 'bpmCT.devs', 'bpmCR.devs', 'bpmCC.devs', 'bpmCB.devs'};
        

        bpmNames = cell(0,0);
        for i=1:length(bpmDevFiles)
            tmpBPMFile = fopen([devsPath bpmDevFiles{i}]);
            tmpBPMNames = textscan(tmpBPMFile,'%s');
            fclose(tmpBPMFile);
            
            for j=1:length(tmpBPMNames)
                bpmNames = [bpmNames; tmpBPMNames{j}];
            end
        end
        strrep(bpmNames,'.','_');
        nBPMs = length(bpmNames);
        
        bpmPropsFile = fopen([devsPath 'bpm.props']);
        bpmProps = textscan(bpmPropsFile,'%s');
        fclose(bpmPropsFile);
        bpmProps = strrep(bpmProps{1},'/','.');
        bpmProps = strrep(bpmProps,'#','.');
        nProps = length(bpmProps);

        % create bpmData struct
        bpmData = struct();
        bpmData.nBPMs = nBPMs;
        bpmData.nProps = nProps;
        bpmData.bpmNames = bpmNames;
        bpmData.bpmProps = bpmProps;
        bpmIsPresent =   false(1,nBPMs);
        propIsPresent = false(1,nProps);
        
        for bpm=1:nBPMs
            for prop=1:nProps
                try
                    % extract signals from the data ('_' separated)
                    tmpData = extractCTFSignalFromMergedData([bpmNames{bpm} bpmProps{prop} '.value'],dataStruct);
                    eval(['bpmData.' bpmNames{bpm} bpmProps{prop} ' = tmpData;']);
                    bpmIsPresent(bpm) = 1;
                    propIsPresent(prop) = 1;
                catch
                    try
                        % extract signals from the data ('.' separated)
                        tmpData = extractCTFSignalFromMergedData([strrep(bpmNames{bpm},'_','.') bpmProps{prop} '.value'],dataStruct);
                        eval(['bpmData.' bpmNames{bpm} bpmProps{prop} ' = tmpData;']);
                        bpmIsPresent(bpm) = 1;
                        propIsPresent(prop) = 1;
                    catch
                        % this bpm or this property not present
%                         if (bpmIsPresent)
%                             propIsPresent(prop) = 0;
%                         end
%                         if (propIsPresent)
%                             bpmIsPresent(bpm) = 0;
%                         end
                        
                    end
                end
            end
        end

        origNBPMs = nBPMs;
        origBPMNames = bpmNames;
        origNProps = nProps;
        origBPMProps = bpmProps;
        origBPMISPresent = bpmIsPresent;
        origPropIsPresent = propIsPresent;
        
        nBPMs = sum(bpmIsPresent);
        nProps = sum(propIsPresent);
        
        bpm = 1;
        while (bpm <= length(bpmNames)) % remove bpms that are not in the data from the list of bpm names
            if (~bpmIsPresent(bpm))
                bpmNames(bpm) = [];
                bpmIsPresent(bpm) = [];
            else
                bpm = bpm+1;
            end
        end
        
        prop = 1;
        while (prop <= length(bpmProps)) % remove bpms that are not in the data from the list of bpm names
            if (~propIsPresent(prop))
                bpmProps(prop) = [];
                propIsPresent(prop) = [];
            else
                prop = prop+1;
            end
        end
        
        if (nBPMs < 0.5*origNBPMs || nProps < 2)
            dataIsGood = 0;
            return;
        end
        
        bpmData.nBPMs = nBPMs;
        bpmData.bpmNames = bpmNames;
        bpmData.nProps = nProps;
        bpmData.bpmProps = bpmProps;    

        % perform some basic check on data to see if the file is good
        propIsSamples = strfind(bpmProps,'Samples');
        for prop=1:nProps
            if (isempty(propIsSamples{prop}))
                continue; % only check Samples (not MeanAtCursor)
            else    
                for bpm=1:nBPMs
                    nSamplesProp = eval(['length(bpmData.' bpmNames{bpm} bpmProps{prop} ');']);

                    if (nSamplesProp < 100)
                        dataIsGood = 0;
                        return;
                    end
                end   
            end  
        end
        dataIsGood = 1;

%         if (isempty(testBPM1) || isnan(min(testBPM1(:))) || isempty(testBPM2) || isnan(min(testBPM2(:))))
%             dataIsGood = 0;
%         else
%             dataIsGood = 1;
%         end
        
    catch % any other issue assume the file is bad
        bpmData = struct();
        dataIsGood = 0;
    end

end

