classdef closure_Class < handle
    %frascatiCalibrationClass 
    % Class used to run calibrations of the Frascati phase monitors.
    
    properties
        % data saving
        saveData = true;
        saveDir = '/user/ctf3op/PhaseFeedforward/SlowCorrection/';
        scanTimeStamp;
        
        nAvg = 10; % no. of pulses acquired for each point

        % correctors
        corr465Offsets = [-2.0 -1.0 1.0 2.0]; % offsets to apply to corrector 465
        minCorrRatio = -2.0; % min765 offset = min465 offset * minCorr765Ratio
        maxCorrRatio = -0.5; % max765 offset = max465 offset * maxCorr765Ratio
        n765Offsets = 10; % no. of points in scan (for each setting of 465)
        n465Offsets;
        corr765Offsets;
        scanCorrRatios;
        init465;
        init765;

        % bpms
        bpmNames = {...
            'CC.SVBPM0435',...
            'CC.SVBPI0535',...
            'CC.SVBPI0645',...
            'CC.SVBPI0685',...
            'CC.SVBPI0735',...
            'CC.SVBPM0845',...
            'CC.SVBPM0930',...
            'CB.SVBPM0150',...
            'CB.SVBPS0210',...
            'CB.SVBPS0250',...
            'CB.SVBPS0310'...
        };
        bpmClosureIndices = 6:11; % use these bpms to calculate orbit closure
        nBPMs;        
        bpmNSamples;
        
        % BPM data arrays
        initOrbitSet;
        initBPMH;
        initBPMS;
        scanBPMH;
        scanBPMS;
        
        % Scan results struct
        scanResults;
        
        % monitor object
        dataMonitor; 

    end
    
    properties (Access=private)
        
        cycleName = 'SCT.USER.SETUP';
        beamEnabledSignal = 'CL.STBPM0402S/Samples#samples';
        bpmAcqName = 'Samples#samples';
        pulseIndex;
        index465;
        index765;
        
    end
    
    methods
        % constructor
        function obj=closure_Class()
            addpath('/nfs/cs-ccr-nfs6/vol29/share/ctf/data/ctfmod/MatLab/PhaseFeedforward');
        end
        
        function startScan(obj)
            obj.initialiseArrays();
            obj.makeMonitor();
            obj.dataMonitor.start();
        end
        
        function pauseScan(obj)
            obj.dataMonitor.pauseOn();
        end
        function resumeScan(obj)
            obj.dataMonitor.pauseOff();
        end
        function stopScan(obj)
            obj.dataMonitor.stop();
        end
          
        function changeCorrs(obj)
            fprintf('Changing correctors...\n');

            if (obj.index765 == obj.n765Offsets)
                obj.index465 = obj.index465 + 1;
                obj.index765 = 1;
            else
                obj.index765 = obj.index765 + 1;
            end

            toSet465 = obj.corr465Offsets(obj.index465) + obj.init465;
            toSet765 = obj.corr765Offsets(obj.index765) + obj.init765;
            
            setCorrValues(toSet465,toSet765);
            
            fprintf('done!\n');
        end
        
        function processScanResults(obj)
            fprintf('Processing scan results...');
            obj.scanResults =...
                closure_processScanResults(...
                    obj.corr465Offsets,...
                    obj.corr765Offsets,...
                    obj.bpmClosureIndices,...
                    obj.initBPMH,...
                    obj.initBPMS,...
                    obj.scanBPMH,...
                    obj.scanBPMS...
                );
            fprintf('done!\n');
        end
        
        function saveScanData(obj)
            fprintf('Saving data...');

            saveName = [obj.saveDir 'closure_' obj.scanTimeStamp '.mat'];
            
            varsToSave = {... 
                'scanResults',...  
                'initOrbitSet',...
                'initBPMH',...
                'initBPMS',...
                'scanBPMH',...
                'scanBPMS',...
                'saveData',...
                'saveDir',...
                'scanTimeStamp',...        
                'nAvg',...     
                'corr465Offsets',... 
                'minCorrRatio',...
                'maxCorrRatio',...
                'n765Offsets',...
                'n465Offsets',...
                'corr765Offsets',...
                'scanCorrRatios',...
                'init465',...
                'init765',...
            };            
            eval(sprintf('%s = obj.%s;',varsToSave{1},varsToSave{1}));
            
            save(saveName,varsToSave{1});
            for i=2:length(varsToSave)
                eval(sprintf('%s = obj.%s;',varsToSave{i},varsToSave{i}));
                save(saveName,varsToSave{i},'-append');
            end
            
            fprintf('done!\n');
            
        end
        
        function saveRefFile(obj)
            refFile = fopen('refRatio','w');
            fprintf(refFile,obj.scanResults.bestEstRatio);
            fclose(refFile);
        end
  
    end
    methods(Access=private)
        % called by dataMonitor when new acquisition received
        function dataReceived(obj, dataStruct)
                    
            % skip first 5 acquisitions
            if (obj.dataMonitor.recordedCycles == 2)
                fprintf('Waiting for good acquisition... ');
            end
            if (obj.dataMonitor.recordedCycles < 5)
                return;
            end
            if (obj.dataMonitor.recordedCycles==5)
                fprintf('done!\n');
            end
            
            % check beam enabled
            transmission = matlabDataAndSignalsHelper.simpleExtractSingleSignal(dataStruct,obj.beamEnabledSignal);
            
            % DEBUG - run even if beam disabled
%             transmission = NaN;
%             warning('DEBUG: Transmission set to NaN to force scan to continue. Search for DEBUG in class and comment appropriate lines to stop this.');
            
            if (min(transmission) > -2) % NB: If min(transmission) returns NaN, beam is assumed to be present.
                obj.dataMonitor.pauseOn();
                fprintf('Beam Disabled. Waiting to continue scan...\n');
                pause(20);
                obj.dataMonitor.pauseOff();
                return;
            else
                obj.pulseIndex = obj.pulseIndex+1;
                
                % save data to arrays (different arrays depending on
                % whether initial orbit is being acquired or not)
                if (obj.initOrbitSet)
                    obj.saveDataToScanArrays(dataStruct);
                else
                    obj.saveDataToInitArrays(dataStruct);              
                end
                
                % if acquired nAvg pulses for these corrector settings,
                % change the correctors or stop the scan and process
                % results if it's the end
                if (obj.pulseIndex == obj.nAvg)
                    obj.dataMonitor.pauseOn();
                    
                    obj.pulseIndex = 0;
                    
                    if (obj.initOrbitSet)
                        if (obj.index465==obj.n465Offsets && obj.index765==obj.n765Offsets) % end of the scan
                            fprintf('-------- END OF SCAN --------\n');
                            obj.dataMonitor.stop();
                            fprintf('Resetting initial corrector values...');
                            setCorrValues(obj.init465,obj.init765); % reset initial corrector values
                            fprintf('done!\n');
                            obj.processScanResults();                    
                            if (obj.saveData)
                                obj.saveScanData();
                                obj.saveRefFile();
                            end
                            return;
                        end
                        
                    else
                        obj.initOrbitSet = 1;
                    end
                    
                    obj.changeCorrs();                   
                    
                    if (obj.saveData)
                        obj.saveScanData();
                    end
                    
                    obj.dataMonitor.pauseOff();
                    
                end
                
            end
            
        end
        
        function saveDataToScanArrays(obj,dataStruct)
            [tmpBPMH,tmpBPMS] = obj.extractBPMsFromDataStruct(dataStruct);
            obj.scanBPMH(obj.index465,obj.index765,obj.pulseIndex,:,:) = tmpBPMH;
            obj.scanBPMS(obj.index465,obj.index765,obj.pulseIndex,:,:) = tmpBPMS;        
        end
        
        function saveDataToInitArrays(obj,dataStruct)
            [tmpBPMH,tmpBPMS] = obj.extractBPMsFromDataStruct(dataStruct);
            obj.initBPMH(obj.pulseIndex,:,:) = tmpBPMH;
            obj.initBPMS(obj.pulseIndex,:,:) = tmpBPMS;           
        end
        
        function [bpmH,bpmS] = extractBPMsFromDataStruct(obj,dataStruct)
            bpmH = NaN(obj.nBPMs,obj.bpmNSamples);
            bpmS = NaN(obj.nBPMs,obj.bpmNSamples);
            for i=1:obj.nBPMs
                [tmpBPMH, tmpBPMS, ~] = extractBPMFromCTFData( strrep(obj.bpmNames{i},'.','_'), dataStruct);
                bpmH(i,:) = tmpBPMH;
                bpmS(i,:) = tmpBPMS;
            end
        end
        
        function initialiseArrays(obj) % pre-define arrays etc.
            fprintf('Initialising arrays...');
            obj.scanResults = struct();
            
            tmpSName = [obj.bpmNames{1} 'S/' obj.bpmAcqName];
            tmpBPMS = matlabJapc.staticGetSignal(obj.cycleName,tmpSName);
            obj.bpmNSamples = length(tmpBPMS);
            obj.nBPMs = length(obj.bpmNames);
            
            obj.scanCorrRatios = linspace(obj.minCorrRatio,obj.maxCorrRatio,obj.n765Offsets);
            obj.n465Offsets = length(obj.corr465Offsets);
            
            obj.corr765Offsets = NaN(obj.n465Offsets,obj.n765Offsets);
            for ind465Val = 1:obj.n465Offsets
                tmp465Offset = obj.corr465Offsets(ind465Val);
                obj.corr765Offsets(ind465Val,:) = linspace(tmp465Offset.*obj.minCorrRatio,...
                                                           tmp465Offset.*obj.maxCorrRatio,...
                                                           obj.n765Offsets);
            end
            
            obj.index465 = 1;
            obj.index765 = 0;
            
            [tmp465,tmp765] = getCorrValues();
            obj.init465 = tmp465;
            obj.init765 = tmp765;
    
            obj.initOrbitSet = 0;
            obj.initBPMH = NaN(obj.nAvg,obj.nBPMs,obj.bpmNSamples);
            obj.initBPMS = NaN(obj.nAvg,obj.nBPMs,obj.bpmNSamples);
            obj.scanBPMH = NaN(obj.n465Offsets,obj.n765Offsets,obj.nAvg,obj.nBPMs,obj.bpmNSamples);
            obj.scanBPMS = NaN(obj.n465Offsets,obj.n765Offsets,obj.nAvg,obj.nBPMs,obj.bpmNSamples);
            obj.pulseIndex = 0;
            
            obj.index465 = 1;
            obj.index765 = 1;

            obj.scanTimeStamp = getTimeStampString();
            
            fprintf('done!\n');

        end
        
        function makeMonitor(obj)
            fprintf('Making japc monitor object...\n');
            tmpSignals = cell(1,(2*obj.nBPMs)+1);
            i=1;
            for bpm=1:obj.nBPMs
                tmpSignals{i} = [obj.bpmNames{bpm} 'S/' obj.bpmAcqName];
                i=i+1;
                tmpSignals{i} = [obj.bpmNames{bpm} 'H/' obj.bpmAcqName];
                i=i+1;
            end
            tmpSignals{i} = obj.beamEnabledSignal;
            
            obj.dataMonitor = matlabJapcMonitor(obj.cycleName, tmpSignals, ...
                @(e)dataReceived(obj,e),'');
            
            obj.dataMonitor.useFastStrategy(1);
            
            fprintf('done!\n');

        end
        
    end
    
end

