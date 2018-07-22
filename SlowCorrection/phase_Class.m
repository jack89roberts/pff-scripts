classdef phase_Class < handle
    %frascatiCalibrationClass 
    % Class used to run calibrations of the Frascati phase monitors.
    
    properties
        % data saving
        saveData = true;
        saveDir = '/user/ctf3op/PhaseFeedforward/SlowCorrection/';
        scanTimeStamp;
        
        nAvg = 20; % no. of pulses acquired for each point

        % correctors
        corr465Offsets = [-2.0 -1.0 0.0 1.0 2.0]; % offsets to apply to corrector 465
        corrRatio;
        nCorrOffsets;
        corr765Offsets;
        init465;
        init765;
        
        calFile = '../FrascatiCalibrations/frascatiCalibrationConstants';
        scanDiodes;
        scanMixers;
        
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
        corrIndex;
        nMons = 3;
        nSamples;
        
    end
    
    methods
        % constructor
        function obj=phase_Class(myCorrRatio)
            obj.corrRatio = myCorrRatio;
            addpath('/nfs/cs-ccr-nfs6/vol29/share/ctf/data/ctfmod/MatLab/PhaseFeedforward');
        end
        
        function startScan(obj)
            obj.initialiseArrays();
            obj.makeMonitor();
            obj.changeCorrs();
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
            obj.corrIndex = obj.corrIndex + 1;

            toSet465 = obj.corr465Offsets(obj.corrIndex) + obj.init465;
            toSet765 = obj.corr765Offsets(obj.corrIndex) + obj.init765;
            
            setCorrValues(toSet465,toSet765);
            fprintf('done!\n');
        end
        
        
        function processScanResults(obj)
            fprintf('Processing scan results...');
            obj.scanResults =   phase_processScanResults(...
                                    obj.corr465Offsets,...
                                    obj.calFile,...
                                    obj.scanDiodes,...
                                    obj.scanMixers...
                                );
            fprintf('done!\n');
        end
        
        function saveScanData(obj)
            fprintf('Saving data...');
            
            saveName = [obj.saveDir 'phase_' obj.scanTimeStamp '.mat'];
            
            varsToSave = {... 
                'scanResults',...
                'saveData',...
                'saveDir',...
                'scanTimeStamp',...
                'nAvg',...
                'corr465Offsets',...
                'corrRatio',...
                'nCorrOffsets',...
                'corr765Offsets',...
                'init465',...
                'init765',...
                'calFile',...
                'scanDiodes',...
                'scanMixers',...
                'cycleName',...
                'beamEnabledSignal',...
                'bpmAcqName',...
                'pulseIndex',...
                'corrIndex',...
                'nMons',...
                'nSamples'... 
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
            fprintf('Saving reference file...');
            refFile = fopen('refPhasePerAmp','w');
            fprintf(refFile,'%.5f',obj.scanResults.fitCoeffs(1));
            fclose(refFile);
            fprintf('done!\n');
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
                fprintf('Acquired pulse %d of %d\n',obj.pulseIndex,obj.nAvg);
                % save data to arrays (different arrays depending on
                % whether initial orbit is being acquired or not)
                obj.saveDataToScanArrays(dataStruct);
                
                % if acquired nAvg pulses for these corrector settings,
                % change the correctors or stop the scan and process
                % results if it's the end
                if (obj.pulseIndex == obj.nAvg)
                    obj.dataMonitor.pauseOn();
                    
                    obj.pulseIndex = 0;
                    
                    if (obj.corrIndex==obj.nCorrOffsets) % end of the scan
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
                    
                    else  % change the correctors and continue the scan              
                        obj.changeCorrs();
                        
                        if (obj.saveData)
                            obj.saveScanData();
                        end
                        obj.dataMonitor.pauseOff();  
                        
                    end                    
                    
                end
                
            end
            
        end
        
        function saveDataToScanArrays(obj,dataStruct)
            [tmpMixers,tmpDiodes] = extractMixerDiode(dataStruct);
            obj.scanMixers(obj.corrIndex,obj.pulseIndex,:,:) = tmpMixers;
            obj.scanDiodes(obj.corrIndex,obj.pulseIndex,:,:) = tmpDiodes;        
        end
        
        function initialiseArrays(obj) % pre-define arrays etc.
            fprintf('Initialising arrays...');
            obj.scanResults = struct();
            
            obj.nCorrOffsets = length(obj.corr465Offsets);
            obj.corr765Offsets = obj.corr465Offsets.*obj.corrRatio;
          
            [tmp465,tmp765] = getCorrValues();
            obj.init465 = tmp465;
            obj.init765 = tmp765;
            
            tmpSigName = 'CT.SCOPE01.CH01/Acquisition#value';
            tmpSig = matlabJapc.staticGetSignal(obj.cycleName,tmpSigName);
            obj.nSamples = length(tmpSig);
            obj.scanMixers = NaN(obj.nCorrOffsets,obj.nAvg,obj.nMons,obj.nSamples);
            obj.scanDiodes = NaN(obj.nCorrOffsets,obj.nAvg,obj.nMons,obj.nSamples);
            
            obj.corrIndex = 0;
            obj.pulseIndex = 0;          

            obj.scanTimeStamp = getTimeStampString();
            fprintf('done!\n');
        end
        
        function makeMonitor(obj)
            fprintf('Making japc monitor object...\n');
            tmpSignals =makePropertyListFromFiles('../devices/scope.devs','../devices/scope.props');
            tmpSignals{end+1} = obj.beamEnabledSignal;
            
            obj.dataMonitor = matlabJapcMonitor(obj.cycleName, tmpSignals, ...
                @(e)dataReceived(obj,e),'');
            
            obj.dataMonitor.useFastStrategy(1);
            fprintf('done!\n');

        end
        
    end
    
end

