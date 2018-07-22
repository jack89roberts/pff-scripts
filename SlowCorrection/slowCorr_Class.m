classdef slowCorr_Class < handle
    %frascatiCalibrationClass 
    % Class used to run calibrations of the Frascati phase monitors.
    
    properties
        % data saving
        saveData = true;
        saveDir = '/user/ctf3op/PhaseFeedforward/SlowCorrection/';
        saveEvery = 200;
        scanTimeStamp;
        
        nAvg = 20; % correction updates every nAvg pusles
        nHistory = 1000; % no. of pulses to keep in history

        gain = 1.0; % if calculated phase offset is x degrees, correction attempted is gain*x.

        % correctors        
        maxCorrVal = 5.0; % don't let correctors go beyond this value
        corrRatio = [];
        init465;
        init765;
        resetCorrs = 1;
        history465;
        history765;
        
        phaseShiftPerAmp = [];
        
        calibrationConstants = [];
        useMixerOverSqrtDiode = [];
        
        initPhase;
        phaseHistory;
        meanPhaseHistory;
        phaseSampleRange = [];
        
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
        bpmAxisNames;
        bpmSHistory;
        bpmHHistory;
        initBPMS;
        initBPMH;
        bpmSampleRange = [];
        
        monColours = [0 0 1; 1 0 0; 0 0.8 0];

        % monitor object
        dataMonitor; 

    end
    
    properties (Access=private)
        calFile = '../FrascatiCalibrations/frascatiCalibrationConstants';
        cycleName = 'SCT.USER.SETUP';
        beamEnabledSignal = 'CL.STBPM0402S/Samples#samples';
        bpmAcqName = 'Samples#samples';
        pulseIndex;
        nMons = 3;
        nBPMs;
        bpmNSamples;
        frascatiNSamples;
        corr465Name = 'CC.DHF0465';
        corr765Name = 'CC.DHF0765';
        corrAcqName = 'Acquisition#currentAverage';
        plotWindow;
        
    end
    
    methods
        % constructor
        function obj=slowCorr_Class()            
            addpath('/nfs/cs-ccr-nfs6/vol29/share/ctf/data/ctfmod/MatLab/PhaseFeedforward');
        end
        
        function loadRefValues(obj)
            fprintf('Loading reference files...');
            obj.corrRatio = importdata('refRatio');
            obj.phaseShiftPerAmp = importdata('refPhasePerAmp');
            fprintf('done!\n');
        end
        
        function loadCalibrations(obj)
            fprintf('Loading calibration constants...');
            [obj.calibrationConstants,obj.useMixerOverSqrtDiode] = loadFrascatiCalibrationConstants(obj.calFile);
            fprintf('done!\n');
        end
        
        function startCorrection(obj)
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
        function stopCorrection(obj)
            fprintf('Stopping japc monitor...');
            obj.dataMonitor.stop();
            fprintf('done!\n');
            if (obj.resetCorrs)                
                fprintf('Resetting initial corrector values...');
                setCorrValues(obj.init465,obj.init765); % reset initial corrector values
                fprintf('done!\n');
            end
            
        end
          
        function changeCorrs(obj)
            fprintf('Changing correctors...');
            
            % should move this to the function for updating arrays
            if (obj.pulseIndex <= obj.nHistory)
                lastPhases = obj.phaseHistory(:,(obj.pulseIndex-obj.nAvg+1):obj.pulseIndex);
            else
                lastPhases = obj.phaseHistory(:,(obj.nHistory-obj.nAvg+1):obj.nHistory);
            end
            lastMeanPhases = nanmean(lastPhases,2);
            if (obj.pulseIndex <= obj.nHistory)
                obj.meanPhaseHistory(1,(obj.pulseIndex-obj.nAvg+1):obj.pulseIndex) = lastMeanPhases(1);
                obj.meanPhaseHistory(2,(obj.pulseIndex-obj.nAvg+1):obj.pulseIndex) = lastMeanPhases(2);
                obj.meanPhaseHistory(3,(obj.pulseIndex-obj.nAvg+1):obj.pulseIndex) = lastMeanPhases(3);
            else
                obj.meanPhaseHistory(1,(obj.nHistory-obj.nAvg+1):obj.nHistory) = lastMeanPhases(1);
                obj.meanPhaseHistory(2,(obj.nHistory-obj.nAvg+1):obj.nHistory) = lastMeanPhases(2);
                obj.meanPhaseHistory(3,(obj.nHistory-obj.nAvg+1):obj.nHistory) = lastMeanPhases(3);
            end

            phaseOffset = lastMeanPhases(3);
            diff465 = phaseOffset ./ obj.phaseShiftPerAmp;
            
            if (obj.pulseIndex<= obj.nHistory)
                toSet465 = obj.history465(obj.pulseIndex) + obj.gain*diff465;
                toSet765 = obj.history765(obj.pulseIndex) + obj.gain*diff465*obj.corrRatio;
            else
                toSet465 = obj.history465(obj.nHistory) + obj.gain*diff465;
                toSet765 = obj.history765(obj.nHistory) + obj.gain*diff465*obj.corrRatio;               
            end

            if (toSet465 > obj.maxCorrVal)
                toSet465 = obj.maxCorrVal;
            elseif (toSet465 < -obj.maxCorrVal)
                toSet465 = -obj.maxCorrVal;
            end
            if (toSet765 > obj.maxCorrVal)
                toSet765 = obj.maxCorrVal;
            elseif (toSet765 < -obj.maxCorrVal)
                toSet765 = -obj.maxCorrVal;
            end
            
            setCorrValues(toSet465,toSet765);
            fprintf('done!\n');
        end
        
        function calculateInit(obj)
            fprintf('Calculating initial phase and orbit...');
            tmpPhases = obj.phaseHistory(:,1:obj.nAvg);
            obj.initPhase = nanmean(tmpPhases,2);
            
            tmpBPMS = obj.bpmSHistory(:,1:obj.nAvg);
            obj.initBPMS = nanmean(tmpBPMS,2);
            
            tmpBPMH = obj.bpmHHistory(:,1:obj.nAvg);
            obj.initBPMH = nanmean(tmpBPMH,2);
            
            for mon=1:obj.nMons
                obj.phaseHistory(mon,:) = obj.phaseHistory(mon,:)-obj.initPhase(mon);
            end
            fprintf('done!\n');
        end
        
        function updatePlots(obj)
            set(0,'CurrentFigure',obj.plotWindow);            
            
            subplot(2,3,1:3);
            for mon=1:obj.nMons
                plot(obj.phaseHistory(mon,:),'--','Color',obj.monColours(mon,:),'LineWidth',0.25);
                hold all;
                plot(obj.meanPhaseHistory(mon,:),'Color',obj.monColours(mon,:),'LineWidth',4)
            end
            title('Phase');
            legLabels = {'Mon1','Mon1 mean','Mon2','Mon2 mean','Mon3','Mon3 mean'};
            xlabel('Pulse No.');
            ylabel('Phase [degrees]');
            if (obj.pulseIndex < obj.nHistory)
                xlim([0 obj.pulseIndex]);
                if (obj.pulseIndex > obj.nAvg)
                    for mon=1:obj.nMons
                        legLabels{2*mon} = [legLabels{2*mon} sprintf(' (last: %.2f^o)',obj.meanPhaseHistory(mon,obj.pulseIndex-obj.nAvg))];
                    end
                end
            else
                for mon=1:obj.nMons
                    legLabels{2*mon} = [legLabels{2*mon} sprintf(' (last: %.2f^o)',obj.meanPhaseHistory(mon,obj.nHistory-obj.nAvg))];
                end
            end
            legend(legLabels,'Location','EastOutside');
            hold off;
            
            subplot(2,3,4);
            plot(obj.history465);
            hold all;
            plot(obj.history765);
            title('Corrector Currents')
            xlabel('Pulse No.');
            ylabel('Current [A]');
            if (obj.pulseIndex < obj.nHistory)
                xlim([0 obj.pulseIndex]);
            end
            legend('465','765');
            hold off;
            
            subplot(2,3,5);
            plot(obj.initBPMH,'k','LineWidth',2);
            hold all;            
            if (obj.pulseIndex < obj.nHistory)
                plot(obj.bpmHHistory(:,obj.pulseIndex))
            else
                plot(obj.bpmHHistory(:,obj.nHistory));
            end
            title('BPM H Orbit');
            xlabel('BPM Index');
            ylabel('Position [mm]');
            xlim([1 obj.nBPMs]);
            legend('Initial','Current');
            hold off;
            
            subplot(2,3,6);
            plot(obj.initBPMS,'k','LineWidth',2);
            hold all;            
            if (obj.pulseIndex < obj.nHistory)
                plot(obj.bpmSHistory(:,obj.pulseIndex))
            else
                plot(obj.bpmSHistory(:,obj.nHistory));
            end
            title('BPM Transmission');
            xlabel('BPM Index');
            ylabel('Current [A]');
            xlim([1 obj.nBPMs]);
            set(gca,'XTick',1:3:obj.nBPMs);
            set(gca,'XTickLabel',{obj.bpmAxisNames{1:3:11}});
            legend('Initial','Current');
            hold off;
        end
        
        function saveHistory(obj)
            fprintf('Saving data...');
            obj.scanTimeStamp = getTimeStampString();
            saveName = [obj.saveDir 'slowCorr_' obj.scanTimeStamp '.mat'];
            
            varsToSave = {... 
                'saveDir',...
                'saveEvery',...
                'scanTimeStamp',...        
                'nAvg',...
                'nHistory',...
                'gain',...     
                'maxCorrVal',...
                'corrRatio',...
                'init465',...
                'init765',...
                'resetCorrs',...
                'history465',...
                'history765' ,...       
                'phaseShiftPerAmp',...        
                'calibrationConstants',...
                'useMixerOverSqrtDiode',...        
                'initPhase',...
                'phaseHistory',...
                'meanPhaseHistory',...
                'phaseSampleRange',...        
                'bpmNames',...
                'bpmSHistory',...
                'bpmHHistory',...
                'initBPMS',...
                'initBPMH',...
                'bpmSampleRange',...
                'pulseIndex'...
            };            
            eval(sprintf('%s = obj.%s;',varsToSave{1},varsToSave{1}));
            
            save(saveName,varsToSave{1});
            for i=2:length(varsToSave)
                eval(sprintf('%s = obj.%s;',varsToSave{i},varsToSave{i}));
                save(saveName,varsToSave{i},'-append');
            end
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
                fprintf('Beam Disabled. Waiting to continue...\n');
                pause(20);
                obj.dataMonitor.pauseOff();
                return;
            else
                obj.pulseIndex = obj.pulseIndex+1;
                
                % save data to arrays (different arrays depending on
                % whether initial orbit is being acquired or not)
                obj.saveDataToArrays(dataStruct);
                
                % plot every 2 pulses
                if (mod(obj.pulseIndex,2)==0)
                    obj.updatePlots();
                end
                
                % calculate initial phase/orbit after nAvg pulses
                if (obj.pulseIndex == obj.nAvg)
                    obj.calculateInit();
                end
                
                % update correction every nAvg pulses
                if (mod(obj.pulseIndex,obj.nAvg) == 0)
                    obj.dataMonitor.pauseOn();
                    obj.changeCorrs();
                    obj.dataMonitor.pauseOff();                                              
                end
                
                % save data every saveEvery pulses
                if (obj.saveData && mod(obj.pulseIndex,obj.saveEvery) == 0)
                    obj.dataMonitor.pauseOn();
                    obj.saveHistory();
                    obj.dataMonitor.pauseOff();                                              
                end
                
            end
            
        end
        
        function saveDataToArrays(obj,dataStruct)
            tmpPhases = obj.extractPhasesFromDataStruct(dataStruct);            
            [tmpBPMH,tmpBPMS] = obj.extractBPMsFromDataStruct(dataStruct);
            [tmp465,tmp765] = obj.extractCorrsFromDataStruct(dataStruct);
            
            meanPhases = nanmean(tmpPhases(:,obj.phaseSampleRange),2);
            meanBPMH = nanmean(tmpBPMH(:,obj.bpmSampleRange),2);
            meanBPMS = nanmean(tmpBPMS(:,obj.bpmSampleRange),2);
            
            if (obj.pulseIndex <= obj.nHistory)
                obj.phaseHistory(:,obj.pulseIndex) = meanPhases;
                obj.bpmSHistory(:,obj.pulseIndex) = meanBPMS;
                obj.bpmHHistory(:,obj.pulseIndex) = meanBPMH;
                obj.history465(obj.pulseIndex) = tmp465;
                obj.history765(obj.pulseIndex) = tmp765;
            else
                obj.phaseHistory(:,1:(obj.nHistory-1)) = obj.phaseHistory(:,2:obj.nHistory);
                obj.meanPhaseHistory(:,1:(obj.nHistory-1)) = obj.meanPhaseHistory(:,2:obj.nHistory);
                obj.bpmSHistory(:,1:(obj.nHistory-1)) = obj.bpmSHistory(:,2:obj.nHistory);
                obj.bpmHHistory(:,1:(obj.nHistory-1)) = obj.bpmHHistory(:,2:obj.nHistory);
                obj.history465(1:(obj.nHistory-1)) = obj.history465(2:obj.nHistory);
                obj.history765(1:(obj.nHistory-1)) = obj.history765(2:obj.nHistory);
                obj.phaseHistory(:,obj.nHistory) = meanPhases;
                obj.bpmSHistory(:,obj.nHistory) = meanBPMS;
                obj.bpmHHistory(:,obj.nHistory) = meanBPMH;
                obj.history465(obj.nHistory) = tmp465;
                obj.history765(obj.nHistory) = tmp765;
                obj.meanPhaseHistory(:,obj.nHistory) = NaN;
            end
        end
        
        function [val465,val765] = extractCorrsFromDataStruct(obj,dataStruct)
            name465 = [strrep(obj.corr465Name,'.','_') '.' strrep(obj.corrAcqName, '#','.') '.value'];
            name765 = [strrep(obj.corr765Name,'.','_') '.' strrep(obj.corrAcqName, '#','.') '.value'];           
            val465 = extractCTFSignalFromMergedData(name465,dataStruct);
            val765 = extractCTFSignalFromMergedData(name765,dataStruct);
        end
        
        function phases = extractPhasesFromDataStruct(obj,dataStruct)
            
            [mixers,diodes] = extractMixerDiode(dataStruct);
            if (isempty(obj.frascatiNSamples))
                obj.frascatiNSamples = size(mixers,3);
            end
            
            phases = NaN(obj.nMons,obj.frascatiNSamples);
            for mon=1:obj.nMons
                m = squeeze(mixers(mon,:,:));
                d = squeeze(diodes(mon,:,:));
                [d,m] = getAlignedXCorr(d','end',{m'});
                m = m{1};
               
                calAmp = obj.calibrationConstants(mon,1);
                calOff = obj.calibrationConstants(mon,4);
                
                if (obj.useMixerOverSqrtDiode)
                    phases(mon,:) = getPhaseMixerDiode(m,d,calAmp,calOff)-obj.initPhase(mon);
                else
                    phases(mon,:) = getPhaseMixerDiode(m,[],calAmp,calOff)-obj.initPhase(mon);
                end
            end
 
        end
        
        function [bpmH,bpmS] = extractBPMsFromDataStruct(obj,dataStruct)

            for i=1:obj.nBPMs
                [tmpBPMH, tmpBPMS, ~] = extractBPMFromCTFData( strrep(obj.bpmNames{i},'.','_'), dataStruct);
                
                [tmpBPMS,tmpBPMH] = getAlignedXCorr(tmpBPMS,'end',{tmpBPMH});
                tmpBPMH = tmpBPMH{1};
                
                if (isempty(obj.bpmNSamples))
                    obj.bpmNSamples = length(tmpBPMS);
                    bpmH = NaN(obj.nBPMs,obj.bpmNSamples);
                    bpmS = NaN(obj.nBPMs,obj.bpmNSamples);
                end
                
                bpmH(i,:) = tmpBPMH;
                bpmS(i,:) = tmpBPMS;
            end
        end       
        
        function initialiseArrays(obj) % pre-define arrays etc.
            if (isempty(obj.corrRatio))
                obj.loadRefValues();
            end
            if (isempty(obj.calibrationConstants))
                obj.loadCalibrations();
            end
            if (isempty(obj.phaseSampleRange))
                obj.phaseSampleRange = getAlignedFrascatiSampleRange();
            end
            if (isempty(obj.bpmSampleRange))
                obj.bpmSampleRange = getBPMSampleRange();
            end
            close all;
                                  
            [obj.init465,obj.init765] = getCorrValues();
            
            obj.nBPMs = length(obj.bpmNames);
            obj.bpmAxisNames = cell(1,obj.nBPMs);
            for bpm=1:obj.nBPMs
                obj.bpmAxisNames{bpm} = obj.bpmNames{bpm}(end-2:end);
            end
            obj.phaseHistory = NaN(obj.nMons,obj.nHistory);
            obj.meanPhaseHistory = NaN(obj.nMons,obj.nHistory);
            obj.bpmSHistory = NaN(obj.nBPMs,obj.nHistory);
            obj.bpmHHistory = NaN(obj.nBPMs,obj.nHistory);
          
            obj.history465 = NaN(1,obj.nHistory);
            obj.history765 = NaN(1,obj.nHistory);
            
            obj.initPhase = zeros(1,obj.nMons);
            obj.initBPMS = zeros(1,obj.nBPMs);
            obj.initBPMH = zeros(1,obj.nBPMs);
            
            obj.pulseIndex = 0;
            
            obj.plotWindow = figure();

        end
        
        function makeMonitor(obj)
 
            tmpSignals = cell(1,2*obj.nBPMs);
            i=1;
            for bpm=1:obj.nBPMs
                tmpSignals{i} = [obj.bpmNames{bpm} 'S/' obj.bpmAcqName];
                i=i+1;
                tmpSignals{i} = [obj.bpmNames{bpm} 'H/' obj.bpmAcqName];
                i=i+1;
            end
            tmpSignals = [tmpSignals makePropertyListFromFiles('../devices/scope.devs','../devices/scope.props')];
            tmpSignals{end+1} = obj.beamEnabledSignal;
            tmpSignals{end+1} = [obj.corr465Name '/' obj.corrAcqName];
            tmpSignals{end+1} = [obj.corr765Name '/' obj.corrAcqName];

            obj.dataMonitor = matlabJapcMonitor(obj.cycleName, tmpSignals, ...
                @(e)dataReceived(obj,e),'');
            
            obj.dataMonitor.useFastStrategy(1);

        end
        
    end
    
end

