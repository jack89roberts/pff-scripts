classdef OnlineR56Opt < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        saveData = true;
        saveDir = '/user/ctf3op/PhaseFeedforward/R56/';
        saveName;
        
        beamEnergy = 121;
        cycleName = 'SCT.USER.SETUP';
        
        wiggleGun = false;
        initGunCurrent;
        stepGunCurrent = 0.05;
        maxOffsetGunCurrent = 0.05;
        gunDevName = 'CL.CTF3-GUN/CCV4';
        gunCurrentValues;
        gunCurrentChangeEvery;
        gunCurrentValueIndex = 1;

        nAvg = 80; % no. pulses to average over when calculating stats
        pulseIndex = 0;
        
        ccBPMName = 'CC.STBPM0930S'; % bpm to check transmission
        ctBPMName = 'CT.STBPI0608'; % bpm to check transmission
        bpmSampleRange = [];
      
        initQuadCurrents;
        initR56 = 0.0;
        currentR56;
        maxR56 = 0.6;
        minR56 = -0.1;
        r56Step = 0.05;
        stepDirection = 1;
        
        phaseSampleRange = [];
        subtractPhase = 0;
        subtractUpstreamPhase = 0;
                
        lastCTPhases;
        lastCBPhases;        
        lastCTPositions;
        lastCCTransmissions;
        
        cbPhaseMeanHistory = [];
        cbPhaseJitterHistory = [];
        r56History = [];
        ccTransmissionHistory = [];
        ccTransmissionJitterHistory = [];
        ctPhaseJitterHistory = [];
        ctPhaseMeanHistory  = [];
        ctPositionHistory = [];
        ctPositionJitterHistory = [];
        phaseCorrHistory = [];
        ctPosCTPhaseCorrHistory = [];
        ctPosCBPhaseCorrHistory = [];
        
        calConsts;
        calUseDio;
        
        monitor;
        
        currentFigure;
        historyFigure;
    end
    
    methods
        % constructor
        function obj=OnlineR56Opt()
            addpath('/nfs/cs-ccr-nfs6/vol29/share/ctf/data/ctfmod/MatLab/PhaseFeedforward');
            
            fprintf('Initialising arrays...\n');
            obj.lastCBPhases = NaN(1,obj.nAvg);
            obj.lastCCTransmissions = NaN(1,obj.nAvg);
            %obj.lastGunCurrents = NaN(1,obj.nAvg);
            obj.lastCTPhases = NaN(1,obj.nAvg);
            obj.lastCTPositions = NaN(1,obj.nAvg);
            
            obj.currentR56 = obj.initR56;
            
            if (obj.saveData)
                dV = datevec(now);
                dataTimeStamp = sprintf('%02d%02d%02d_%02d%02d',dV(1),dV(2),dV(3),dV(4),dV(5));
                obj.saveName = [obj.saveDir dataTimeStamp '.mat'];
            end
            
            fprintf('Acquiring gun current...\n');
            obj.initGunCurrent = matlabJapc.staticGetSignal(obj.cycleName,obj.gunDevName);
            
            fprintf('Loading calibration constants...\n');
            [obj.calConsts, obj.calUseDio] = loadFrascatiCalibrationConstants('/nfs/cs-ccr-nfs6/vol29/share/ctf/data/ctfmod/MatLab/PhaseFeedforward/FrascatiCalibrations/frascatiCalibrationConstants');                
    
            fprintf('Making monitor object...\n'); 
            devsFilePath = '/nfs/cs-ccr-nfs6/vol29/share/ctf/data/ctfmod/MatLab/PhaseFeedforward/devices/';
            [signals, ~] = combinePropertyLists(... 
                                            makePropertyListFromFiles([devsFilePath 'frascatiScope.devs'],[devsFilePath 'scope.props'])...
                                         ); 
            signals = [signals 'CL.STBPM0402S/Samples#samples'];            % bpm to check there's beam
            signals = [signals [obj.ccBPMName '/Samples#samples']];            % bpm for transmission check
            signals = [signals [obj.ctBPMName 'S/Samples#samples']];
            signals = [signals [obj.ctBPMName 'H/Samples#samples']];
            obj.monitor = matlabJapcMonitor(obj.cycleName, signals, @(e)obj.dataReceived(e),'');
            obj.monitor.useFastStrategy(1);
            obj.monitor.synchronizedExternalCalls = true; % don't do anything with new data if not finished processing previous data. This will screw up crate alignment if it is needed.
            
            fprintf('Making plot windows...\n');
            obj.currentFigure = figure();
            obj.historyFigure = figure();

        end
        
        function dataReceived(obj,dataStruct)

            %fprintf('Checking beam is present...\n');
            tmpTransmission = extractCTFSignalFromMergedData('CL_STBPM0402S.Samples.samples.value',dataStruct); 
            if min(tmpTransmission) > -2 % no beam in linac if true
                return; % don't update anything
            else
                obj.pulseIndex = obj.pulseIndex + 1;
            end
            
            %fprintf('Saving new data to arrays...\n');
            obj.extractData(dataStruct);
            
            if (obj.pulseIndex == obj.nAvg)
                obj.monitor.pauseOn();

                %fprintf('End of current cycle...\n');
                obj.calculateStats();
                if (length(obj.cbPhaseJitterHistory)>1)
                    obj.changeR56();
                    fprintf('Old phase jitter: %.2f, New phase jitter: %.2f\n',obj.cbPhaseJitterHistory(end-1),obj.cbPhaseJitterHistory(end));
                    fprintf('Old R56: %.2f, New R56 %.2f, Step direction: %.2f\n',obj.r56History(end-1),obj.r56History(end),obj.stepDirection); 
                end
                obj.updateHistoryPlot();
                obj.pulseIndex = 0;
                
                if (obj.saveData)
                    obj.saveDataHistory();
                end
                
                obj.gunCurrentValueIndex = 1;
                obj.setGunCurrent(1);
                
                obj.monitor.pauseOff();
            else
                %fprintf('Changing gun current...\n');
                obj.changeGunCurrent();
         
                obj.updateCurrentPlot();
            end
            
        end
        
        function startOptimisation(obj)
            fprintf('Getting/setting R56 in TL1...\n');
            obj.initQuadCurrents = getCurrentsTL1();
            setR56TL1(obj.initR56,obj.beamEnergy);
            
            if (isempty(obj.phaseSampleRange))
                obj.changeSampleRanges();
            end
            
            fprintf('Calculating gun current values...\n');
            if (obj.wiggleGun)
                minCurrent = obj.initGunCurrent-obj.maxOffsetGunCurrent;
                maxCurrent = obj.initGunCurrent+obj.maxOffsetGunCurrent;
                obj.gunCurrentValues = minCurrent:obj.stepGunCurrent:maxCurrent;
                obj.gunCurrentChangeEvery = round(obj.nAvg./length(obj.gunCurrentValues));
                obj.setGunCurrent(1);
            else
                obj.gunCurrentValues = obj.initGunCurrent;
            end

            fprintf('Starting acquisition...\n');
            obj.monitor.start();
        end
        
        function stopOptimisation(obj)
            fprintf('Stopping acquisition...\n');
            obj.monitor.stop();
            obj.resetInitialGunCurrent();
        end
        
        function resetInitialGunCurrent(obj)
            fprintf('Resetting initial gun current...\n');
            matlabJapc.staticSetSignal(obj.cycleName,obj.gunDevName,obj.initGunCurrent);
        end
        
        function setGunCurrent(obj, valueIndex)
            newCurrent = obj.gunCurrentValues(valueIndex)
            matlabJapc.staticSetSignal(obj.cycleName,obj.gunDevName,newCurrent);
            
            fprintf('New gun current: %.3f\n',newCurrent);
        end
            
        function saveDataHistory(obj)
            dataStruct = struct();
            dataStruct.cbPhaseMeanHistory = obj.cbPhaseMeanHistory;
            dataStruct.cbPhaseJitterHistory = obj.cbPhaseJitterHistory;
            dataStruct.r56History = obj.r56History;
            dataStruct.ccTransmissionHistory = obj.ccTransmissionHistory;
            dataStruct.ccTransmissionJitterHistory = obj.ccTransmissionJitterHistory;
            dataStruct.ctPhaseJitterHistory = obj.ctPhaseJitterHistory;
            dataStruct.ctPhaseMeanHistory  = obj.ctPhaseMeanHistory;
            dataStruct.ctPositionHistory = obj.ctPositionHistory;
            dataStruct.ctPositionJitterHistory = obj.ctPositionJitterHistory;
            dataStruct.phaseCorrHistory = obj.phaseCorrHistory;
            dataStruct.ctPosCTPhaseCorrHistory = obj.ctPosCTPhaseCorrHistory;
            dataStruct.ctPosCBPhaseCorrHistory = obj.ctPosCBPhaseCorrHistory;
            
            dataStruct.beamEnergy = obj.beamEnergy;
            dataStruct.initGunCurrent = obj.initGunCurrent;
            dataStruct.gunCurrentValues = obj.gunCurrentValues;
            dataStruct.nAvg = obj.nAvg; % no. pulses to average over when calculating stats
            dataStruct.ccBPMName = obj.ccBPMName; % bpm to check transmission
            dataStruct.ctBPMName = obj.ctBPMName; % bpm to check transmission
            dataStruct.bpmSampleRange = obj.bpmSampleRange;
            dataStruct.initQuadCurrents = obj.initQuadCurrents;
            dataStruct.initR56 = obj.initR56;
            dataStruct.currentR56 = obj.currentR56;
            dataStruct.r56Step = obj.r56Step;
            dataStruct.phaseSampleRange = obj.phaseSampleRange;
            dataStruct.calConsts = obj.calConsts;
            dataStruct.calUseDio = obj.calUseDio;
            
            save(obj.saveName,'dataStruct');
        end
        
        function changeSampleRanges(obj)
            fprintf('Prompting for sample ranges...\n');
            
            figure;
            tmpDio = matlabJapc.staticGet(obj.cycleName,'CT.SCOPE01.CH05','Acquisition','value');
            tmpMix = matlabJapc.staticGet(obj.cycleName,'CT.SCOPE01.CH06','Acquisition','value');
            [tmpDio,tmpMix] = getAlignedXCorr(double(tmpDio),'end',{double(tmpMix)});
            tmpMix = tmpMix{1};
            subplot(1,2,1)
            plot(tmpDio);
            subplot(1,2,2)
            plot(tmpMix);
            obj.phaseSampleRange = input('Phase: Sample range to use: ');
            
            figure;
            subplot(1,3,1)
            tmpBPM = matlabJapc.staticGet(obj.cycleName,obj.ccBPMName,'Samples','samples');
            tmpBPM = getAlignedXCorr(tmpBPM,'end');
            plot(tmpBPM);
            title(obj.ccBPMName);
            subplot(1,3,2)
            tmpBPMH = matlabJapc.staticGet(obj.cycleName,[obj.ctBPMName 'H'],'Samples','samples');
            tmpBPMS = matlabJapc.staticGet(obj.cycleName,[obj.ctBPMName 'S'],'Samples','samples');
            [tmpBPMS,tmpBPMH] = getAlignedXCorr(tmpBPMS,'end',{tmpBPMH});
            tmpBPMH = tmpBPMH{1};
            plot(tmpBPMS);
            title(obj.ctBPMName);
            subplot(1,3,3)
            plot(tmpBPMH);
            title(obj.ctBPMName);
            obj.bpmSampleRange = input('BPM Sample range to use: ');
            
            close all;
        end
        
        function changeGunCurrent(obj)
            
            if (    obj.wiggleGun &&...
                    obj.pulseIndex>0 &&...
                    mod(obj.pulseIndex,obj.gunCurrentChangeEvery)==0)
                
                obj.monitor.pauseOn();
                
                obj.gunCurrentValueIndex = obj.gunCurrentValueIndex + 1;
                obj.setGunCurrent(obj.gunCurrentValueIndex);
                              
                pause(2);
                obj.monitor.pauseOff();
            end
        end
        
        function extractData(obj,dataStruct)
            [mixers,diodes] = extractMixerDiode(dataStruct);
            if (obj.calUseDio)
                phaseDown = getPhaseMixerDiode(squeeze(mixers(3,:,:)),squeeze(diodes(3,:,:)),obj.calConsts(3,1),obj.calConsts(3,4));
                phaseUp = getPhaseMixerDiode(squeeze(mixers(2,:,:)),squeeze(diodes(2,:,:)),obj.calConsts(2,1),obj.calConsts(2,4));
            else
                phaseDown = getPhaseMixerDiode(squeeze(mixers(3,:,:)),[],obj.calConsts(3,1),obj.calConsts(3,4));
                phaseUp = getPhaseMixerDiode(squeeze(mixers(2,:,:)),[],obj.calConsts(2,1),obj.calConsts(2,4));      
            end
            [~,phaseDown] = getAlignedXCorr(squeeze(diodes(3,:,:))','end',{phaseDown'});
            phaseDown = phaseDown{1};
            
            [~,phaseUp] = getAlignedXCorr(squeeze(diodes(2,:,:))','end',{phaseUp'});
            phaseUp = phaseUp{1};
            
            if (isempty(obj.cbPhaseJitterHistory))
                obj.subtractPhase = nanmean(phaseDown(obj.phaseSampleRange));
                obj.subtractUpstreamPhase = nanmean(phaseUp(obj.phaseSampleRange));
            end
            
            phaseDown = phaseDown - obj.subtractPhase;
            phaseUp = phaseUp - obj.subtractUpstreamPhase;
            
            obj.lastCBPhases(obj.pulseIndex) = nanmean(phaseDown(obj.phaseSampleRange));
            obj.lastCTPhases(obj.pulseIndex) = nanmean(phaseUp(obj.phaseSampleRange));
            
            %obj.lastGunCurrents(obj.pulseIndex) = obj.gunCurrent;
            
            tmpTransmission = extractCTFSignalFromMergedData([ strrep(obj.ccBPMName,'.','_') '.Samples.samples.value'],dataStruct);
            tmpTransmission = getAlignedXCorr(tmpTransmission,'end');
            obj.lastCCTransmissions(obj.pulseIndex) = nanmean(tmpTransmission(obj.bpmSampleRange));
            
            tmpPosition = extractCTFSignalFromMergedData([ strrep(obj.ctBPMName,'.','_') 'H.Samples.samples.value'],dataStruct);
            tmpTransmission = extractCTFSignalFromMergedData([ strrep(obj.ctBPMName,'.','_') 'S.Samples.samples.value'],dataStruct);
            [~,tmpPosition] = getAlignedXCorr(tmpTransmission,'end',{tmpPosition});
            tmpPosition = tmpPosition{1};
            obj.lastCTPositions(obj.pulseIndex) = nanmean(tmpPosition(obj.bpmSampleRange));

        end
        
        function calculateStats(obj)
            %fprintf('Calculating statistics...\n');
            downstreamPhase = nanmean(obj.lastCBPhases);
            jitterDownstreamPhase = nanstd(obj.lastCBPhases);
            
            upstreamPhase = nanmean(obj.lastCTPhases);
            jitterUpstreamPhase = nanstd(obj.lastCTPhases);
           
            transmission = nanmean(obj.lastCCTransmissions);
            transmissionJitter = nanstd(obj.lastCCTransmissions);
            
            position = nanmean(obj.lastCTPositions);
            positionJitter = nanstd(obj.lastCTPositions);

            phaseCorr = nancorrcoef(obj.lastCTPhases,obj.lastCBPhases);
            ctPosCTPhaseCorr = nancorrcoef(obj.lastCTPositions,obj.lastCTPhases);
            ctPosCBPhaseCorr = nancorrcoef(obj.lastCTPositions,obj.lastCBPhases);

            if (isempty(obj.cbPhaseJitterHistory))
                obj.cbPhaseMeanHistory(1) = downstreamPhase;
                obj.cbPhaseJitterHistory(1) = jitterDownstreamPhase;
                obj.ccTransmissionHistory(1) = transmission;
                obj.ccTransmissionJitterHistory(1) = transmissionJitter;
                obj.ctPhaseJitterHistory(1) = jitterUpstreamPhase;
                obj.ctPhaseMeanHistory(1)  = upstreamPhase;
                obj.r56History(1) = obj.currentR56;            
                obj.ctPositionHistory(1) = position;
                obj.ctPositionJitterHistory(1) = positionJitter;
                obj.phaseCorrHistory(1) = phaseCorr;
                obj.ctPosCTPhaseCorrHistory(1) = ctPosCTPhaseCorr;
                obj.ctPosCBPhaseCorrHistory(1) = ctPosCBPhaseCorr;
            else
                obj.cbPhaseMeanHistory(end+1) = downstreamPhase;
                obj.cbPhaseJitterHistory(end+1) = jitterDownstreamPhase;
                obj.ccTransmissionHistory(end+1) = transmission;
                obj.ccTransmissionJitterHistory(end+1) = transmissionJitter;
                obj.ctPhaseJitterHistory(end+1) = jitterUpstreamPhase;
                obj.ctPhaseMeanHistory(end+1)  = upstreamPhase;
                obj.ctPositionHistory(end+1) = position;
                obj.ctPositionJitterHistory(end+1) = positionJitter;
                obj.phaseCorrHistory(end+1) = phaseCorr;
                obj.ctPosCTPhaseCorrHistory(end+1) = ctPosCTPhaseCorr;
                obj.ctPosCBPhaseCorrHistory(end+1) = ctPosCBPhaseCorr;
            end
            
        end
        
        function changeR56(obj)
            fprintf('Setting new R56...\n');
            
            oldDiff = abs(obj.ctPosCBPhaseCorrHistory(end-1)- obj.ctPosCTPhaseCorrHistory(end-1));
            newDiff = abs(obj.ctPosCBPhaseCorrHistory(end)- obj.ctPosCTPhaseCorrHistory(end));
            
            if (newDiff > oldDiff)
                obj.stepDirection = -obj.stepDirection;  
            end
                        
            obj.currentR56 = obj.currentR56 + obj.stepDirection.*obj.r56Step;
            
            if (obj.currentR56 > obj.maxR56)
                obj.currentR56 = 0.6;
            elseif (obj.currentR56 < obj.minR56);
                obj.currentR56 = -0.1;
            end
            
            obj.r56History(end+1) = obj.currentR56;
            
            setR56TL1(obj.currentR56,obj.beamEnergy);
        end
        
        function updateHistoryPlot(obj)
            set(0,'CurrentFigure',obj.historyFigure);  
            nIterations = length(obj.r56History);
            xAxis = 1:nIterations;
            
            subplot(2,2,1);
            plot(obj.r56History);
            title('R56 Value History');
            xlabel('Iteration');
            ylabel('R56 [m]');
            
            subplot(2,2,2)
            plot(obj.phaseCorrHistory);            
            hold all;
            plot(obj.ctPosCTPhaseCorrHistory);
            plot(obj.ctPosCBPhaseCorrHistory);
            legend('CTPhase-CBPhase','CTPos-CTPhase','CTPos-CBPhase');
            ylabel('Correlation');
            ylim([-1 1]);
            xlabel('Iteration');
            title('Correlation History')
            hold off;
            
            subplot(2,2,3);
            [yAxes,H1,H2] = plotyy(xAxis,obj.ctPositionJitterHistory,xAxis,[obj.ctPhaseJitterHistory;obj.cbPhaseJitterHistory]);
            hold all;
            legend([H1;H2],'CT Position','CT Phase','CB Phase');
            ylabel(yAxes(1),'Position Jitter [mm]');
            ylabel(yAxes(2),'Phase Jitter [degrees]');
            xlabel('Iteration');
            hold off;
            title('Jitter History');
            hold off;
            
            subplot(2,2,4);
            [yAxes,H1,H2] = plotyy(xAxis,obj.ctPositionHistory,xAxis,[obj.ctPhaseMeanHistory;obj.cbPhaseMeanHistory]);
            hold all;
            legend([H1;H2],'CT Position','CT Phase','CB Phase');
            ylabel(yAxes(1),'Position [mm]');
            ylabel(yAxes(2),'Phase [degrees]');
            xlabel('Iteration');
            hold off;
            title('Mean History');
            hold off;
    
        end
        
        function updateCurrentPlot(obj)
            set(0,'CurrentFigure',obj.currentFigure);  
            
            subplot(1,3,1);
            plot(obj.lastCTPhases);
            hold all;
            plot(obj.lastCBPhases);
            legend('Upstream','Downstream');
            title('Current iteration: Mean phase');
            xlabel('Pulse No.');
            ylabel('Phase');
            hold off;
            
            subplot(1,3,2);
            plot(obj.lastCTPositions);
            title(sprintf('Current iteration: Mean position %s',obj.ctBPMName));
            xlabel('Pulse No.');
            ylabel('Phase');
            hold off;
            
            subplot(1,3,3);
            plot(obj.lastCCTransmissions);
            title(obj.ccBPMName);
            xlabel('Pulse No.');
            ylabel('Current [A]');
     
            
        end
    end
    
end

