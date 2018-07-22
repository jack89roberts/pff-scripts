classdef onlinePhaseMonitor < handle
    % onlinePhaseMonitor
    %   Small class used to suscribe some signals and plot them online.
    %   Just create an object as
    %     myobj = subscribeAndPlot(<selector>, <signals cell array>);
    %   and start acquiring data
    %     myobj.startAcquisition();
    
    properties
        % cell list with labels for different signals. It is automatically
        % populated, but you can change it later.
        %signalsLabels;
        % colors for the different plots. Automatically generated, but you
        % can change it later.
        %signalsColors;
        
        % cell with all last data
        lastData;
        lastTreatedPhaseData;
        
        noUpdateIfBeamDisabled = true;
        
        subtractInitPhase = true; % if true, will subtract the mean phase of the first subtractInitNPulses pulses 
        subtractInitNPulses = 10;
        
        % no. of pulses to keep in memory
        lengthPulseHistory = 100;
        lengthArrays; % actual arrays longer than pulse history in order to deal with crate delays
        
        % data arrays
        allTreatedPhaseData;
        allPhaseStatistics;
        
        % Calibration factors etc. for the Frascati monitors
        frascatiCalibrationFactors; % loaded from reference file by default
        %frascatiUseMixerOverSqrtDiode = true; % To implement in the future
        frascatiTimePerSample;
        frascatiStartSamples;
        frascatiEndSamples;
        frascatiSubtractPhases = [0,0,0];
        frascatiNMonitors = 3;
        frascatiNames = {'Mon1','Mon2','Mon3'};
        frascatiTimeAxes;
        plotSampleRange;
        
        % Pets info
        usePets = true;
        petsTimePerSample;
        petsStartSample;
        petsEndSample;
        petsSubtractPhase = 0;
        petsName = 'PETS';
        petsTimeAxis;
        
        % BPM info
        ct285StartSample;
        ct285EndSample;
        ct285TimePerSample = 5.2085;%10.417;
        ct285TimeAxis;
        ct285CrateDelay = 0; 
        
        % for data processing
        alignTo = 'end'; % whether to align signals to start or end of pulse
        
        %matlabJacMonitor object
        datamonitor;
  
    end
    properties (Access=private)
        % cell list of signals in use (DON'T CHANGE IT!).
        signals;
        %monitorIdxs;
        
        useMixerOverSqrtDiode;
        
        currentPulseIndex = 0;
        lastDataIndex = 0;
        skipFirstFiveIndex = 0; % just used to skip the first 5 pulses (acquisition can be unstabled for first few pusles)
        acquisitionsSinceLastMessage = 20; % used to control how often beam disabled messages are printed. DON'T CHANGE VALUE HERE WITHOUT CHANGING BELOW!
        
        % figure handler (don't touch)
        lastPhaseFigureHandle;
        meanPhaseFigureHandle;
        correlationFigureHandle;
        
        goodPulses;
    
    end
    methods
        
        % constructor
        function obj=onlinePhaseMonitor()
            fprintf(1,'Creating phase monitor instance...\n');
            
            addpath('/nfs/cs-ccr-nfs6/vol29/share/ctf/data/ctfmod/MatLab/PhaseFeedforward');
            close all;
            
            cycleName = 'SCT.USER.SETUP';
            
            devsFilePath = '/nfs/cs-ccr-nfs6/vol29/share/ctf/data/ctfmod/MatLab/PhaseFeedforward/devices/';
            if (obj.usePets)
                [obj.signals, ~] = combinePropertyLists(... 
                                            ...%makePropertyListFromFiles([devsFilePath 'bpr.devs'],[devsFilePath 'bpr.props']),...
                                            makePropertyListFromFiles([devsFilePath 'scope.devs'],[devsFilePath 'scope.props'])...
                                         ); 
            else
                [obj.signals, ~] = combinePropertyLists(... 
                                            ...%makePropertyListFromFiles([devsFilePath 'bpr.devs'],[devsFilePath 'bpr.props']),...
                                            makePropertyListFromFiles([devsFilePath 'frascatiScope.devs'],[devsFilePath 'scope.props'])...
                                         ); 
            end
            % bpm to check there's beam
            obj.signals = [obj.signals 'CL.STBPM0402S/Samples#samples'];
            
            % CT.285 for energy-phase correlation
            obj.signals = [obj.signals 'CT.STBPI0622H/Samples#samples'];
            obj.signals = [obj.signals 'CT.STBPI0622S/Samples#samples'];
            
            [obj.frascatiCalibrationFactors, obj.useMixerOverSqrtDiode] = loadFrascatiCalibrationConstants('/nfs/cs-ccr-nfs6/vol29/share/ctf/data/ctfmod/MatLab/PhaseFeedforward/FrascatiCalibrations/frascatiCalibrationConstants');                
            
            obj.datamonitor = matlabJapcMonitor(cycleName, obj.signals, @(e)obj.updatePhaseData(e),'');
            obj.datamonitor.useFastStrategy(1);
            obj.datamonitor.synchronizedExternalCalls = true; % don't do anything with new data if not finished processing previous data. This will screw up crate alignment if it is needed.
 
            fprintf('Done!\n');
        end
        
        function stopAcquisition(obj)
            %obj.lastPhaseFigureHandle = [];
            obj.datamonitor.stop();
        end
        function startAcquisition(obj)            
            if (isempty(obj.frascatiStartSamples) || isempty(obj.frascatiEndSamples))
                error('Please define pulse start and end samples.');
            end
            if (isempty(obj.plotSampleRange))
                error('Please define plot sample range.');
            end
   
            obj.skipFirstFiveIndex = 0;
            obj.datamonitor.start(1);
        end
        function pauseAcquisition(obj)
            obj.datamonitor.pauseOn();
        end
        function resumeAcquisition(obj)
            obj.skipFirstFiveIndex = 0;
            obj.datamonitor.pauseOff();
        end
        function delete(obj)
            obj.stopAcquisition();
            obj.datamonitor.delete();
        end

        function plotSignals(obj)
            %disp('plot signals');
 
            tmpCurrentFigure = get(0,'CurrentFigure');
            
            % plot last acquired phases / changed to mean phase plots
            set(0,'CurrentFigure',obj.lastPhaseFigureHandle);            
 
            subplot(1,3,1);
            for mon=1:obj.frascatiNMonitors 
                tmpPhase = obj.allPhaseStatistics.frascatiMeanPulsePhases(mon,:);
                plot(squeeze(tmpPhase));
                hold all
            end
            if (obj.usePets)
                plot(obj.allPhaseStatistics.petsMeanPulsePhase);
                tmpLegStr = cell(1,obj.frascatiNMonitors+1);
            else
                tmpLegStr = cell(1,obj.frascatiNMonitors);
            end
            
            for mon=1:obj.frascatiNMonitors
                tmpLegStr{mon} = [obj.frascatiNames{mon} sprintf(' (std %.1f^{o})', obj.allPhaseStatistics.frascatiStdPulsePhases(mon))];
            end
            if (obj.usePets)
                tmpLegStr{obj.frascatiNMonitors+1} = [obj.petsName sprintf(' (std %.1f^{o})', obj.allPhaseStatistics.petsStdPulsePhase)];
            end
            legend(tmpLegStr);%legend([obj.frascatiNames obj.petsName]); 
            %legend(tmpLegStr{2},tmpLegStr{3},tmpLegStr{4});
            title('Mean Pulse Phase vs. Time');
            xlabel('Time [pulse no.]');
            ylabel('Phase [12GHz degrees]');
            grid
            hold off
                        
            subplot(1,3,2);
            for mon=1:obj.frascatiNMonitors 
            %for mon=2:obj.frascatiNMonitors 
                %tmpPhase = obj.allPhaseStatistics.frascatiMeanSamplePhases(mon,obj.frascatiStartSamples(mon):obj.frascatiEndSamples(mon));
                tmpPhase = obj.allPhaseStatistics.frascatiMeanSamplePhases;
                tmpTimeAxis = obj.frascatiTimeAxes{mon};
                %tmpTimeAxis = tmpTimeAxis(obj.frascatiStartSamples(mon):obj.frascatiEndSamples(mon));
                plot(tmpTimeAxis,squeeze(tmpPhase));
                hold all
            end
            timeStart = obj.frascatiTimeAxes{3}(obj.frascatiStartSamples(3));
            timeEnd = obj.frascatiTimeAxes{3}(obj.frascatiEndSamples(3));
            plot([timeStart timeStart],get(gca,'YLim'),'k','LineWidth',1);
            plot([timeEnd timeEnd],get(gca,'YLim'),'k','LineWidth',1);

            if (obj.usePets)
                tmpTimeAxis = obj.petsTimeAxis(obj.petsStartSample:obj.petsEndSample);
                plot(tmpTimeAxis,obj.allPhaseStatistics.petsMeanSamplePhase(obj.petsStartSample:obj.petsEndSample));
                legend([obj.frascatiNames obj.petsName]);
            else
                legend(obj.frascatiNames);
            end
            title('Mean Phase Along the Pulse');
            xlabel('Time [ns]');
            ylabel('Phase [12GHz degrees]');
            %legend(obj.frascatiNames{2}, obj.frascatiNames{3}, obj.petsName);
            grid
            xlim([tmpTimeAxis(obj.plotSampleRange(1)) tmpTimeAxis(obj.plotSampleRange(end))]);
            hold off
             
            subplot(1,3,3);
            for mon=1:obj.frascatiNMonitors 
            %for mon=2:obj.frascatiNMonitors 
                %tmpPhase = obj.allPhaseStatistics.frascatiStdSamplePhases(mon,obj.frascatiStartSamples(mon):obj.frascatiEndSamples(mon));
                tmpPhase = obj.allPhaseStatistics.frascatiStdSamplePhases(mon,:);
                tmpTimeAxis = obj.frascatiTimeAxes{mon};
                %tmpTimeAxis = tmpTimeAxis(obj.frascatiStartSamples(mon):obj.frascatiEndSamples(mon));
                plot(tmpTimeAxis,squeeze(tmpPhase));
                hold all
            end
            if (obj.usePets)
                tmpTimeAxis = obj.petsTimeAxis(obj.petsStartSample:obj.petsEndSample);
                plot(tmpTimeAxis,obj.allPhaseStatistics.petsStdSamplePhase(obj.petsStartSample:obj.petsEndSample));
                tmpLegStr = cell(1,obj.frascatiNMonitors+1);
            else
                tmpLegStr = cell(1,obj.frascatiNMonitors);
            end
            
            for mon=1:obj.frascatiNMonitors
                tmpLegStr{mon} = [obj.frascatiNames{mon} sprintf(' (mean %.1f^{o})', obj.allPhaseStatistics.frascatiMeanStdSamplePhases(mon))];
            end
            if (obj.usePets)
                tmpLegStr{obj.frascatiNMonitors+1} = [obj.petsName sprintf(' (mean %.1f^{o})', obj.allPhaseStatistics.petsMeanStdSamplePhase)];
            end
            legend(tmpLegStr);%legend([obj.frascatiNames obj.petsName]);
            %legend(tmpLegStr{2},tmpLegStr{3},tmpLegStr{4});%legend([obj.frascatiNames obj.petsName]);
            title('Std Phase Along the Pulse');
            xlabel('Time [ns]');
            ylabel('Phase [12GHz degrees]');
            grid
            plot([timeStart timeStart],get(gca,'YLim'),'k','LineWidth',1);
            plot([timeEnd timeEnd],get(gca,'YLim'),'k','LineWidth',1);
            xlim([tmpTimeAxis(obj.plotSampleRange(1)) tmpTimeAxis(obj.plotSampleRange(end))]);
            hold off
             
            % mean position plots
            set(0,'CurrentFigure',obj.meanPhaseFigureHandle);                    
 
            subplot(1,3,1);
            tmpPhase = obj.allPhaseStatistics.ct285HMeanPulse;
            plot(tmpPhase);
            title('CT.622H Mean Position vs. Time');
            %size(obj.allPhaseStatistics.ct285HStdPulse)
            legend(sprintf('std: %.2fmm',obj.allPhaseStatistics.ct285HStdPulse))
            xlabel('Time [pulse no.]');
            ylabel('Position [mm]');
            grid
                        
            subplot(1,3,2);
            tmpPhase = obj.allPhaseStatistics.ct285HMeanSample;
            tmpTimeAxis = obj.ct285TimeAxis;
            %tmpTimeAxis = tmpTimeAxis(obj.ct285StartSample:obj.ct285EndSample);
            plot(tmpTimeAxis,tmpPhase);
            title('CT.622H Mean Position Along the Pulse');
            xlabel('Time [ns]');
            ylabel('Position [mm]');
            grid
            hold all;
            plot([timeStart timeStart],get(gca,'YLim'),'k','LineWidth',1);
            plot([timeEnd timeEnd],get(gca,'YLim'),'k','LineWidth',1);
            xlim([tmpTimeAxis(obj.plotSampleRange(1)) tmpTimeAxis(obj.plotSampleRange(end))]);
            hold off;
             
            subplot(1,3,3);
            tmpPhase = obj.allPhaseStatistics.ct285HStdSample;
            tmpTimeAxis = obj.ct285TimeAxis;
            tmpTimeAxis = tmpTimeAxis;
            plot(tmpTimeAxis,tmpPhase);
            title('CT.622H Std Position Along the Pulse');
            legend(sprintf('mean: %.2fmm',obj.allPhaseStatistics.ct285HMeanStdSample));
            xlabel('Time [ns]');
            ylabel('Phase [12GHz degrees]');
            grid
            hold all;
            plot([timeStart timeStart],get(gca,'YLim'),'k','LineWidth',1);
            plot([timeEnd timeEnd],get(gca,'YLim'),'k','LineWidth',1);
            xlim([tmpTimeAxis(obj.plotSampleRange(1)) tmpTimeAxis(obj.plotSampleRange(end))]);
%            yLim = [min(tmpPhase(obj.ct285StartSample:obj.ct285EndSample)) max(tmpPhase(obj.ct285StartSample:obj.ct285EndSample))];
%            ylim(yLim);
            hold off;
             
            % correlation plots
            set(0,'CurrentFigure',obj.correlationFigureHandle);            
            
            subplot(2,2,1)
            tmpPhaseA = obj.allPhaseStatistics.goodMeanFrascati1;
            tmpPhaseB = obj.allPhaseStatistics.goodMeanFrascati2;
            plot(squeeze(tmpPhaseA),squeeze(tmpPhaseB),'o');
            tmpCorr = obj.allPhaseStatistics.corr12_meanPulsePhase;
            tmpGrad = obj.allPhaseStatistics.fit12_meanPulsePhase(1);
            titString = sprintf('Mon1 vs. Mon2: corr=%.2f, grad=%.2f',tmpCorr,tmpGrad);
            title(titString);
            xlabel('Mon1 Phase [12GHz degrees]');
            ylabel('Mon2 Phase [12GHz degrees]');
            grid
            
            subplot(2,2,2)
            tmpPhaseA = obj.allPhaseStatistics.goodMeanFrascati2;
            tmpPhaseB = obj.allPhaseStatistics.goodMeanFrascati3;
            plot(squeeze(tmpPhaseA),squeeze(tmpPhaseB),'o');
            tmpCorr = obj.allPhaseStatistics.corr23_meanPulsePhase;
            tmpGrad = obj.allPhaseStatistics.fit23_meanPulsePhase(1);
            titString = sprintf('Mon2 vs. Mon3: corr=%.2f, grad=%.2f',tmpCorr,tmpGrad);
            title(titString);
            xlabel('Mon2 Phase [12GHz degrees]');
            ylabel('Mon3 Phase [12GHz degrees]');
            grid
            
            subplot(2,2,3)
            tmpPhaseA = obj.allPhaseStatistics.goodMeanCT285H;
            tmpPhaseB = obj.allPhaseStatistics.goodMeanFrascati2;
            plot(squeeze(tmpPhaseA),squeeze(tmpPhaseB),'o');
            tmpCorr = obj.allPhaseStatistics.corrCT285H_Frascati2;
            tmpGrad = obj.allPhaseStatistics.fitCT285H_Frascati2(1);
            titString = sprintf('Mon2 vs. CT.622H: corr=%.2f, grad=%.2f',tmpCorr,tmpGrad);
            title(titString);
            xlabel('CT.285H Position [mm]');
            ylabel('Mon2 Phase [12GHz degrees]');
            grid
            
            subplot(2,2,4)
            tmpPhaseA = obj.allPhaseStatistics.goodMeanCT285H;
            tmpPhaseB = obj.allPhaseStatistics.goodMeanFrascati3;
            plot(squeeze(tmpPhaseA),squeeze(tmpPhaseB),'o');
            tmpCorr = obj.allPhaseStatistics.corrCT285H_Frascati3;
            tmpGrad = obj.allPhaseStatistics.fitCT285H_Frascati3(1);
            titString = sprintf('Mon3 vs. CT.622H: corr=%.2f, grad=%.2f',tmpCorr,tmpGrad);
            title(titString);
            xlabel('CT.285H Position [mm]');
            ylabel('Mon3 Phase [12GHz degrees]');
            grid
            
%             if (obj.usePets)
%                 subplot(2,2,4)
% 
%                 tmpPhaseA = obj.allPhaseStatistics.goodMeanFrascati3;
%                 tmpPhaseB = obj.allPhaseStatistics.goodMeanPets;
%                 plot(squeeze(tmpPhaseA),squeeze(tmpPhaseB),'o');
%                 tmpCorr = obj.allPhaseStatistics.corrPetsFrascati3_meanPulsePhase;
%                 tmpGrad = obj.allPhaseStatistics.fitPetsFrascati3_meanPulsePhase(1);
%                 titString = sprintf('Mon3 vs. PETS: corr=%.2f, grad=%.2f',tmpCorr,tmpGrad);
%                 title(titString);
%                 xlabel('Mon3 Phase [12GHz degrees]');
%                 ylabel('PETS Phase [12GHz degrees]');
%                 grid
%             end

            
            % reset back figure was in use before
            set(0,'CurrentFigure',tmpCurrentFigure);
        end
    end
    
    methods(Access=private)
        
        function updatePhaseData(obj, dataStruct)
            %tic
            % skip the first 5 acquisitions (reduces chance of NaNs etc.
            % from first few pulses)
            set(obj.lastPhaseFigureHandle,'DeleteFcn',@(h,e)obj.stopAcquisition); 
            if (obj.skipFirstFiveIndex < 5)
                obj.skipFirstFiveIndex = obj.skipFirstFiveIndex + 1;
                
                % create figures for future use whilst waiting
                if isempty(obj.lastPhaseFigureHandle)
                    obj.lastPhaseFigureHandle = figure();
                    obj.meanPhaseFigureHandle = figure();
                    obj.correlationFigureHandle = figure();
                end
                
                
                fprintf(1,'Skipping pulse %d of %d to wait for good acquisition.\n',obj.skipFirstFiveIndex,5);
                return;
            end
            
            % don't update if the beam is disabled
            if (obj.noUpdateIfBeamDisabled == true)
                
                %tmpTransmission = extractCTFSignalFromMergedData('CL.STBPM0402S.Samples.samples.value',dataStruct); % JACK - Probably needs to be changed for dots now being replaced with _ in saved data
                tmpTransmission = extractCTFSignalFromMergedData('CL_STBPM0402S.Samples.samples.value',dataStruct); % JACK - Probably needs to be changed for dots now being replaced with _ in saved data
                
                if min(tmpTransmission) > -2 % no beam in linac if true
                    if (obj.acquisitionsSinceLastMessage == 20) % print a message every 20 pulses
                        fprintf(1,'Beam disabled, plots not updating.\n');
                        obj.acquisitionsSinceLastMessage = 0;
                    else
                        obj.acquisitionsSinceLastMessage = obj.acquisitionsSinceLastMessage + 1;
                    end
                    return; % don't update anything
                else
                    if obj.acquisitionsSinceLastMessage~=20
                        obj.acquisitionsSinceLastMessage = 20;
                        fprintf(1,'Beam enabled again, plots updating.\n');
                    end
                end
                
            end
            
            %if (dataStruct.CT.SCOPE01.CH01
                        
            obj.lastData = dataStruct;

            % creates the lastTreatedPhaseData struct
            %disp('treatLastData');
            obj.treatLastData();
            
            % calculates new statistics for data history
            %disp('calculatePhaseStatistics');
            if (obj.currentPulseIndex > obj.subtractInitNPulses)
                obj.calculatePhaseStatistics();
            end
            
            % call plot function (only call every two pulses to avoid
            % backlog)
            if mod(dataStruct.seqNumber,2) == 0
                %disp('plotSignals');
                obj.plotSignals();
            end
            %toc
        end
        
        function treatLastData(obj)
            % if you want to do something on acquired data, here should go
            % your code.
            % The new data is available in obj.lastData
            
            %disp('inside treat last data');
            
            % Calculate the phases etc. from the most recent data
            [frascatiMixers,frascatiDiodes] = extractMixerDiode(obj.lastData); % JACK - Changed extractMixerDiode function to try to deal with '.' being replaced with '_'
            frascatiMixers = squeeze(frascatiMixers);
            frascatiDiodes = squeeze(frascatiDiodes);
            [~,frascatiNSamples] = size(frascatiMixers);            

            frascatiPhases = NaN*ones(obj.frascatiNMonitors,frascatiNSamples);
            for mon = 1:obj.frascatiNMonitors
                [tmpDiodes,tmpMixers] = getAlignedXCorr(squeeze(frascatiDiodes(mon,:,:)),obj.alignTo,{squeeze(frascatiMixers(mon,:,:))});
                frascatiDiodes(mon,:,:) = tmpDiodes;
                frascatiMixers(mon,:,:) = tmpMixers{1};
                
                if (obj.useMixerOverSqrtDiode == true)
                    frascatiPhases(mon,:,:) = getPhaseMixerDiode(frascatiMixers(mon,:,:),...
                                                         frascatiDiodes(mon,:,:),...
                                                         obj.frascatiCalibrationFactors(mon,1),...
                                                         obj.frascatiCalibrationFactors(mon,4)...
                                                        );
                else
                    frascatiPhases(mon,:,:) = getPhaseMixerDiode(frascatiMixers(mon,:,:),...
                                                         [],...
                                                         obj.frascatiCalibrationFactors(mon,1),...
                                                         obj.frascatiCalibrationFactors(mon,4)...
                                                        );
                end
                
                if (obj.subtractInitPhase == true);
                    frascatiPhases(mon,:,:) = frascatiPhases(mon,:,:) - obj.frascatiSubtractPhases(mon);
                end
            end
            frascatiPhases = squeeze(frascatiPhases);
            
            if (obj.usePets)
                petsI = double(extractCTFSignalFromMergedData('CE_SCOPE03_CH01.Acquisition.value.value',obj.lastData));  % JACK - Changed to try to deal with '.' being replaced with '_'
                petsQ = double(extractCTFSignalFromMergedData('CE_SCOPE03_CH02.Acquisition.value.value',obj.lastData));
                petsNSamples = length(petsI);

                [petsQ,petsI] = getAlignedXCorr(petsQ,obj.alignTo,{petsI});
                petsI = petsI{1};

                petsISensitivity = extractCTFSignalFromMergedData('CE_SCOPE03_CH01.Acquisition.sensitivity.value',obj.lastData);
                petsQSensitivity = extractCTFSignalFromMergedData('CE_SCOPE03_CH02.Acquisition.sensitivity.value',obj.lastData);

                petsI = double(petsI).*petsISensitivity;
                petsQ = double(petsQ).*petsQSensitivity;
                [petsPhase,petsPower] = getPhaseIQ(petsI, petsQ);
                petsNSamples = length(petsPhase);
                petsPhase = squeeze(petsPhase);
                if (obj.subtractInitPhase == true);
                    petsPhase = petsPhase - obj.petsSubtractPhase;
                end
            end
            
            ct285H = extractCTFSignalFromMergedData('CT_STBPI0622H.Samples.samples.value',obj.lastData);
            ct285S = extractCTFSignalFromMergedData('CT_STBPI0622S.Samples.samples.value',obj.lastData);
            ct285NSamples = length(ct285S);
            [ct285S,ct285H] = getAlignedXCorr(ct285S,obj.alignTo,{ct285H});
            ct285H = ct285H{1};
            if (obj.subtractInitPhase == true)
                % think about subtracting init position for BPM
            end
            
            if (obj.currentPulseIndex == 0) % things to initialise for the first acquired pulse
                % get time information for all the monitors and create time
                % axes (relative to mon1)
                obj.frascatiTimePerSample = extractCTFSignalFromMergedData('CT_SCOPE01_CH02.Acquisition.sampleInterval.value',obj.lastData); % JACK - Changed to try to deal with '.' being replaced with '_'
                if (obj.usePets)
                    obj.petsTimePerSample = extractCTFSignalFromMergedData('CE_SCOPE03_CH01.Acquisition.sampleInterval.value',obj.lastData);
                end
                obj.lengthArrays = obj.lengthPulseHistory + obj.ct285CrateDelay;
                
                tmpFrascatiTimeAxis = linspace(0, (obj.frascatiTimePerSample.*frascatiNSamples)-obj.frascatiTimePerSample, frascatiNSamples);
                obj.frascatiTimeAxes = {};
                for mon=1:obj.frascatiNMonitors
                    obj.frascatiTimeAxes{mon} = tmpFrascatiTimeAxis - tmpFrascatiTimeAxis(obj.frascatiStartSamples(mon));
                end
                
                if (obj.usePets)
                    tmpPetsTimeAxis = linspace(0, (obj.petsTimePerSample.*petsNSamples)-obj.petsTimePerSample, petsNSamples);
                    obj.petsTimeAxis = tmpPetsTimeAxis - tmpPetsTimeAxis(obj.petsStartSample);
                end
                
                tmpCT285Axis = linspace(0, (obj.ct285TimePerSample.*ct285NSamples)-obj.ct285TimePerSample, ct285NSamples);
                obj.ct285TimeAxis = tmpCT285Axis - tmpCT285Axis(obj.ct285StartSample);
                
                % initialise arrays
                obj.allTreatedPhaseData = struct();
                obj.allTreatedPhaseData.frascatiMixers = NaN(obj.frascatiNMonitors,obj.lengthArrays,frascatiNSamples);
                obj.allTreatedPhaseData.frascatiDiodes = NaN(obj.frascatiNMonitors,obj.lengthArrays,frascatiNSamples);
                obj.allTreatedPhaseData.frascatiPhases = NaN(obj.frascatiNMonitors,obj.lengthArrays,frascatiNSamples);
                if (obj.usePets)
                    obj.allTreatedPhaseData.petsI = NaN(obj.lengthArrays,petsNSamples);
                    obj.allTreatedPhaseData.petsQ = NaN(obj.lengthArrays,petsNSamples);
                    obj.allTreatedPhaseData.petsPower = NaN(obj.lengthArrays,petsNSamples);
                    obj.allTreatedPhaseData.petsPhase = NaN(obj.lengthArrays,petsNSamples);
                end
                
                obj.allTreatedPhaseData.strippedFrascatiPhases = [];
                obj.allTreatedPhaseData.strippedFrascatiDiodes  = [];
                obj.allTreatedPhaseData.strippedFrascatiMixers = [];
                if (obj.usePets)
                    obj.allTreatedPhaseData.strippedPetsPhase = [];
                end
                obj.allTreatedPhaseData.ct285S = NaN(obj.lengthArrays,ct285NSamples);
                obj.allTreatedPhaseData.ct285H = NaN(obj.lengthArrays,ct285NSamples);
                
                obj.goodPulses = true(1,obj.lengthArrays);
                
                obj.allPhaseStatistics = struct();
                obj.currentPulseIndex = 1;
            
            else
                obj.currentPulseIndex = obj.currentPulseIndex + 1;
            end          
            % delete the first pulse (oldest pulse) in the data history
            % when we reach the set limit for no. pulses to keep
            if (obj.currentPulseIndex > obj.lengthPulseHistory)
                obj.lastDataIndex = obj.lengthPulseHistory;

                
                % JACK!!! Need to figure out how to deal with crate delays
                % in indexing!
                
                % circshift shifts all the elements in an array by a given
                % number of indices over a given dimension
                % frascati arrays - 2nd dimensions is pulse no.
                % pets arrays - 1st dimension is pulse no.
                % AIM IS TO SHIFT THE OLDEST DATA TO THE LAST INDEX, SO THAT IT
                % IS OVERWRITTEN BELOW 
%                 obj.allTreatedPhaseData.frascatiMixers =  circshift(obj.allTreatedPhaseData.frascatiMixers,-1,2);
%                 obj.allTreatedPhaseData.frascatiDiodes = circshift(obj.allTreatedPhaseData.frascatiDiodes,-1,2);
%                 obj.allTreatedPhaseData.frascatiPhases = circshift(obj.allTreatedPhaseData.frascatiPhases,-1,2);
%                 obj.allTreatedPhaseData.petsI = circshift(obj.allTreatedPhaseData.petsI,-1,1);
%                 obj.allTreatedPhaseData.petsQ = circshift(obj.allTreatedPhaseData.petsQ,-1,1);
%                 obj.allTreatedPhaseData.petsPower = circshift(obj.allTreatedPhaseData.petsPower,-1,1);
%                 obj.allTreatedPhaseData.petsPhase = circshift(obj.allTreatedPhaseData.petsPhase,-1,1);

                obj.allTreatedPhaseData.frascatiMixers(:,1:end-1,:) =  obj.allTreatedPhaseData.frascatiMixers(:,2:end,:);
                obj.allTreatedPhaseData.frascatiDiodes(:,1:end-1,:) = obj.allTreatedPhaseData.frascatiDiodes(:,2:end,:);
                obj.allTreatedPhaseData.frascatiPhases(:,1:end-1,:) = obj.allTreatedPhaseData.frascatiPhases(:,2:end,:);
                if (obj.usePets)
                    obj.allTreatedPhaseData.petsI(1:end-1,:) = obj.allTreatedPhaseData.petsI(2:end,:);
                    obj.allTreatedPhaseData.petsQ(1:end-1,:) = obj.allTreatedPhaseData.petsQ(2:end,:);
                    obj.allTreatedPhaseData.petsPower(1:end-1,:) = obj.allTreatedPhaseData.petsPower(2:end,:);
                    obj.allTreatedPhaseData.petsPhase(1:end-1,:) = obj.allTreatedPhaseData.petsPhase(2:end,:);
                end
                obj.allTreatedPhaseData.ct285S(1:end-1,:) = obj.allTreatedPhaseData.ct285S(2:end,:);
                obj.allTreatedPhaseData.ct285H(1:end-1,:) = obj.allTreatedPhaseData.ct285H(2:end,:);


            else
                obj.lastDataIndex = obj.currentPulseIndex;
            end
            
            ct285Index = obj.lastDataIndex + obj.ct285CrateDelay;

            obj.allTreatedPhaseData.frascatiMixers(:,obj.lastDataIndex,:) = frascatiMixers;
            obj.allTreatedPhaseData.frascatiDiodes(:,obj.lastDataIndex,:) = frascatiDiodes;
            obj.allTreatedPhaseData.frascatiPhases(:,obj.lastDataIndex,:) = frascatiPhases;
            if (obj.usePets)
                obj.allTreatedPhaseData.petsI(obj.lastDataIndex,:) = petsI;
                obj.allTreatedPhaseData.petsQ(obj.lastDataIndex,:) = petsQ;
                obj.allTreatedPhaseData.petsPower(obj.lastDataIndex,:) = petsPower;
                obj.allTreatedPhaseData.petsPhase(obj.lastDataIndex,:) = petsPhase;
            end
            obj.allTreatedPhaseData.ct285S(ct285Index,:) = ct285S;
            obj.allTreatedPhaseData.ct285H(ct285Index,:) = ct285H;

            
            if (obj.currentPulseIndex==obj.subtractInitNPulses && obj.subtractInitPhase==true)
                % calculate phase to subtract 
                obj.calculateInitPhase();
                %obj.currentPulseIndex = obj.currentPulseIndex + 1;
            end
            
            
            
        end
        
        function calculateInitPhase(obj)
            %disp('oops')
            % calculate the mean initial phase to subtract (if option
            % selected)
            
            tmpMeanFrascatiPhases = zeros(1,obj.frascatiNMonitors);
            
            for mon=1:obj.frascatiNMonitors
                tmpFrascatiPhase = squeeze(obj.allTreatedPhaseData.frascatiPhases(mon,obj.goodPulses,:));
                tmpFrascatiPhase = nanmean(nanmean(tmpFrascatiPhase(:,obj.frascatiStartSamples(mon):obj.frascatiEndSamples(mon))));
                tmpMeanFrascatiPhases(mon) = tmpFrascatiPhase;
            end
            %disp(tmpMeanFrascatiPhases)
            obj.frascatiSubtractPhases = tmpMeanFrascatiPhases;
            
            if (obj.usePets)
                tmpPetsPhase = obj.allTreatedPhaseData.petsPhase(obj.goodPulses,:);
                obj.petsSubtractPhase = nanmean(nanmean(tmpPetsPhase(:,obj.petsStartSample:obj.petsEndSample)));
            end
%             % subtract the mean phase from previously acquired data
%             for mon=1:obj.frascatiNMonitors
%                 obj.allTreatedPhaseData.frascatiPhases(mon,:,:) = obj.allTreatedPhaseData.frascatiPhases(mon,:,:) -obj.frascatiSubtractPhases(mon);
%             end
%             obj.allTreatedPhaseData.petsPhase = obj.allTreatedPhaseData.petsPhase - obj.petsSubtractPhase;
        end
        
        function calculatePhaseStatistics(obj)
            %if (~isempty(obj.allTreatedPhaseData.strippedFrascatiPhases))
%             if (false)
%                 strippedFrascatiPhases = obj.allTreatedPhaseData.strippedFrascatiPhases;
%                 strippedFrascatiDiodes = obj.allTreatedPhaseData.strippedFrascatiDiodes;
%                 strippedFrascatiMixers = obj.allTreatedPhaseData.strippedFrascatiMixers;
%                 strippedPetsPhase = obj.allTreatedPhaseData.strippedPetsPhase;
%             else
            
%             end

            % Remove bad pulses
            
            strippedFrascatiPhases = obj.allTreatedPhaseData.frascatiPhases;
            strippedFrascatiDiodes  = obj.allTreatedPhaseData.frascatiDiodes;
            strippedFrascatiMixers = obj.allTreatedPhaseData.frascatiMixers;
            if (obj.usePets)
                strippedPetsPhase = obj.allTreatedPhaseData.petsPhase;  
            end
            
            for mon=1:obj.frascatiNMonitors
                tmpDiode = squeeze(strippedFrascatiDiodes(mon,:,:));
                tmpMixer = squeeze(strippedFrascatiMixers(mon,:,:));
                tmpPhase = squeeze(strippedFrascatiPhases(mon,:,:));
                tmpSampRange = obj.frascatiStartSamples(mon):obj.frascatiEndSamples(mon);
                [tmpDiode,tmpOther] = removeBadPulses(tmpDiode,tmpSampRange,{tmpMixer,tmpPhase});
                
                strippedFrascatiDiodes(mon,:,:) = tmpDiode;
                strippedFrascatiMixers(mon,:,:) = tmpOther{1};
                strippedFrascatiPhases(mon,:,:) = tmpOther{2};
            end
            
            if (obj.usePets)
                tmpPETSRange = obj.petsStartSample:obj.petsEndSample;
                strippedPetsPhase = removeBadPulses(strippedPetsPhase,tmpPETSRange);
            end
 
            strippedCT285S = obj.allTreatedPhaseData.ct285S;
            strippedCT285H = obj.allTreatedPhaseData.ct285H;
            tmp285Range = obj.ct285StartSample:obj.ct285EndSample;
            [strippedCT285S,strippedCT285H] = removeBadPulses(strippedCT285S,tmp285Range,{strippedCT285H});
            strippedCT285H = strippedCT285H{1};
            
            % Means and standard deviations for the Frascati monitors
            frascatiMeanSamplePhases = squeeze(nanmean(strippedFrascatiPhases,2));
            frascatiStdSamplePhases = squeeze(nanstd(strippedFrascatiPhases,0,2));
            frascatiMeanSampleDiodes = squeeze(nanmean(strippedFrascatiDiodes,2));
            frascatiStdSampleDiodes = squeeze(nanstd(strippedFrascatiDiodes,0,2));
            frascatiMeanSampleMixers = squeeze(nanmean(strippedFrascatiMixers,2));
            frascatiStdSampleMixers = squeeze(nanstd(strippedFrascatiMixers,0,2));

            frascatiMeanStdSamplePhases = NaN*ones(1,obj.frascatiNMonitors);
            frascatiMeanStdSampleDiodes = NaN*ones(1,obj.frascatiNMonitors);
            frascatiMeanStdSampleMixers = NaN*ones(1,obj.frascatiNMonitors);
            

            allFrascatiMeanPulsePhases = NaN*ones(obj.frascatiNMonitors,obj.lengthArrays);
            allFrascatiMeanPulseDiodes = NaN*ones(obj.frascatiNMonitors,obj.lengthArrays);
            allFrascatiMeanPulseMixers = NaN*ones(obj.frascatiNMonitors,obj.lengthArrays);
            
            for mon=1:obj.frascatiNMonitors
                frascatiMeanStdSamplePhases(mon) = squeeze(nanmean(frascatiStdSamplePhases(mon,obj.frascatiStartSamples(mon):obj.frascatiEndSamples(mon)),2));
                frascatiMeanStdSampleDiodes(mon) = squeeze(nanmean(frascatiStdSampleDiodes(mon,obj.frascatiStartSamples(mon):obj.frascatiEndSamples(mon)),2));
                frascatiMeanStdSampleMixers(mon) = squeeze(nanmean(frascatiStdSampleMixers(mon,obj.frascatiStartSamples(mon):obj.frascatiEndSamples(mon)),2));
                
                allFrascatiMeanPulsePhases(mon,:) = squeeze(nanmean(obj.allTreatedPhaseData.frascatiPhases(mon,:,obj.frascatiStartSamples(mon):obj.frascatiEndSamples(mon)),3));
                allFrascatiMeanPulseDiodes(mon,:) = squeeze(nanmean(obj.allTreatedPhaseData.frascatiDiodes(mon,:,obj.frascatiStartSamples(mon):obj.frascatiEndSamples(mon)),3));
                allFrascatiMeanPulseMixers(mon,:) = squeeze(nanmean(obj.allTreatedPhaseData.frascatiMixers(mon,:,obj.frascatiStartSamples(mon):obj.frascatiEndSamples(mon)),3));

                frascatiMeanPulsePhases(mon,:) = squeeze(nanmean(strippedFrascatiPhases(mon,:,obj.frascatiStartSamples(mon):obj.frascatiEndSamples(mon)),3));
                frascatiMeanPulseDiodes(mon,:) = squeeze(nanmean(strippedFrascatiDiodes(mon,:,obj.frascatiStartSamples(mon):obj.frascatiEndSamples(mon)),3));
                frascatiMeanPulseMixers(mon,:) = squeeze(nanmean(strippedFrascatiMixers(mon,:,obj.frascatiStartSamples(mon):obj.frascatiEndSamples(mon)),3));
                

            end
           
            frascatiStdPulsePhases = nanstd(frascatiMeanPulsePhases,0,2); 
            frascatiStdPulseDiodes = nanstd(frascatiMeanPulseDiodes,0,2); 
            frascatiStdPulseMixers = nanstd(frascatiMeanPulseMixers,0,2); 
            
             
            % Means and standard deviations for the PETS
            if (obj.usePets)
                petsMeanSamplePhase = nanmean(strippedPetsPhase); 
                allPetsMeanPulsePhase = nanmean(obj.allTreatedPhaseData.petsPhase(:,obj.petsStartSample:obj.petsEndSample),2);
                petsMeanPulsePhase = nanmean(strippedPetsPhase(:,obj.petsStartSample:obj.petsEndSample),2);
                petsStdSamplePhase = nanstd(strippedPetsPhase);
                petsMeanStdSamplePhase = nanmean(petsStdSamplePhase(obj.petsStartSample:obj.petsEndSample));
                petsStdPulsePhase = nanstd(petsMeanPulsePhase);
            end
            
            % Means and standard deviations for CT.BPM285
            ct285SMeanSample = nanmean(strippedCT285S);
            ct285SStdSample = nanstd(strippedCT285S);
            ct285SMeanStdSample = nanmean(ct285SStdSample(obj.ct285StartSample:obj.ct285EndSample));
            ct285SMeanPulse = nanmean(strippedCT285S(:,obj.ct285StartSample:obj.ct285EndSample),2);
            ct285SStdPulse = nanstd(strippedCT285S(:,obj.ct285StartSample:obj.ct285EndSample),0,2);
         
            ct285HMeanSample = nanmean(strippedCT285H);
            ct285HStdSample = nanstd(strippedCT285H);
            ct285HMeanStdSample = nanmean(ct285HStdSample(obj.ct285StartSample:obj.ct285EndSample));
            ct285HMeanPulse = nanmean(strippedCT285H(:,obj.ct285StartSample:obj.ct285EndSample),2);
            ct285HStdPulse = nanstd(ct285HMeanPulse);%strippedCT285H(:,obj.ct285StartSample:obj.ct285EndSample),0,2);

            % Correlations 
            % remove outlying means - not working?
            goodMeanCT285S = removeBadPulses(ct285SMeanPulse,1)';
            goodMeanCT285H = removeBadPulses(ct285HMeanPulse,1)';
            goodMeanFrascati1 = removeBadPulses(frascatiMeanPulsePhases(1,:),1);
            goodMeanFrascati2 = removeBadPulses(frascatiMeanPulsePhases(2,:),1);
            goodMeanFrascati3 = removeBadPulses(frascatiMeanPulsePhases(3,:),1);
            if (obj.usePets)
                goodMeanPets = removeBadPulses(petsMeanPulsePhase,1);
            end
            
            corr12_meanPulsePhase = nancorrcoef(goodMeanFrascati1,goodMeanFrascati2);           
            corr13_meanPulsePhase = nancorrcoef(goodMeanFrascati1,goodMeanFrascati3);
            corr23_meanPulsePhase = nancorrcoef(goodMeanFrascati2,goodMeanFrascati3);
            
            if (obj.usePets)
                corrPetsFrascati3_meanPulsePhase = nancorrcoef(goodMeanFrascati3,petsMeanPulsePhase');   
                corrPetsFrascati3_meanPulsePhase = corrPetsFrascati3_meanPulsePhase(1,2);
            end
            
            corrCT285S_Frascati1 = nancorrcoef(goodMeanCT285S, goodMeanFrascati1);
            corrCT285S_Frascati2 = nancorrcoef(goodMeanCT285S, goodMeanFrascati2);
            corrCT285S_Frascati3 = nancorrcoef(goodMeanCT285S, goodMeanFrascati3);

            corrCT285H_Frascati1 = nancorrcoef(goodMeanCT285H, goodMeanFrascati1);
            corrCT285H_Frascati2 = nancorrcoef(goodMeanCT285H, goodMeanFrascati2);
            corrCT285H_Frascati3 = nancorrcoef(goodMeanCT285H, goodMeanFrascati3);
            
            % Fits
            fit12_meanPulsePhase = nanpolyfit(goodMeanFrascati1,goodMeanFrascati2,1);
            fit13_meanPulsePhase = nanpolyfit(goodMeanFrascati1,goodMeanFrascati3,1);
            fit23_meanPulsePhase = nanpolyfit(goodMeanFrascati2,goodMeanFrascati3,1);
            if (obj.usePets)
                fitPetsFrascati3_meanPulsePhase = nanpolyfit(goodMeanFrascati3,petsMeanPulsePhase',1);
            end
            
            fitCT285S_Frascati1 = nanpolyfit(goodMeanCT285S, goodMeanFrascati1,1);
            fitCT285S_Frascati2 = nanpolyfit(goodMeanCT285S, goodMeanFrascati2,1);
            fitCT285S_Frascati3 = nanpolyfit(goodMeanCT285S, goodMeanFrascati3,1);
            
            fitCT285H_Frascati1 = nanpolyfit(goodMeanCT285H, goodMeanFrascati1,1);
            fitCT285H_Frascati2 = nanpolyfit(goodMeanCT285H, goodMeanFrascati2,1);
            fitCT285H_Frascati3   = nanpolyfit(goodMeanCT285H, goodMeanFrascati3,1);         

            
            % remove bad pulses and calculate new mean phase
%             if (obj.subtractInitPhase == true)
%                 calculateInitPhase(obj);
%             end
%                 
%             for p=1:obj.lengthPulseHistory
%                 isAGoodPulse = true;
% 
%                 for mon=1:obj.frascatiNMonitors
%                     %fprintf(1,'%d > %d\n', abs(frascatiMeanPulsePhases(mon,p)), 3*frascatiStdPulsePhases(mon));
%                     if (abs(allFrascatiMeanPulsePhases(mon,p)) > 3*frascatiStdPulsePhases(mon)) 
%                         isAGoodPulse = false;
%                         break;
%                     end
%                 end
% 
%                 if (isAGoodPulse)
%                     if(abs(allPetsMeanPulsePhase(p)) > 3*petsStdPulsePhase); 
%                         isAGoodPulse = false;
%                     end
%                 end
% 
%                 if (isAGoodPulse)
%                     obj.goodPulses(p) = 1;
%                 else
%                     obj.goodPulses(p) = 0;
%                 end
% 
%             end
%             %disp(sum(obj.goodPulses))
%             
%             obj.allTreatedPhaseData.strippedFrascatiPhases = obj.allTreatedPhaseData.frascatiPhases(:,obj.goodPulses,:);
%             obj.allTreatedPhaseData.strippedFrascatiDiodes  = obj.allTreatedPhaseData.frascatiDiodes(:,obj.goodPulses,:);
%             obj.allTreatedPhaseData.strippedFrascatiMixers = obj.allTreatedPhaseData.frascatiMixers(:,obj.goodPulses,:);
%             obj.allTreatedPhaseData.strippedPetsPhase = obj.allTreatedPhaseData.petsPhase(obj.goodPulses,:);

            % save everything to the allPhaseStatistics struct
            obj.allPhaseStatistics.frascatiMeanSamplePhases = frascatiMeanSamplePhases;
            obj.allPhaseStatistics.frascatiStdSamplePhases = frascatiStdSamplePhases;
            obj.allPhaseStatistics.frascatiMeanStdSamplePhases = frascatiMeanStdSamplePhases;
            obj.allPhaseStatistics.frascatiMeanPulsePhases = frascatiMeanPulsePhases;
            obj.allPhaseStatistics.frascatiStdPulsePhases = frascatiStdPulsePhases;
            
            obj.allPhaseStatistics.frascatiMeanSampleDiodes = frascatiMeanSampleDiodes;
            obj.allPhaseStatistics.frascatiStdSampleDiodes = frascatiStdSampleDiodes;
            obj.allPhaseStatistics.frascatiMeanStdSampleDiodes = frascatiMeanStdSampleDiodes;
            obj.allPhaseStatistics.frascatiMeanPulseDiodes = frascatiMeanPulseDiodes;
            obj.allPhaseStatistics.frascatiStdPulseDiodes = frascatiStdPulseDiodes;
            
            obj.allPhaseStatistics.frascatiMeanSampleMixers = frascatiMeanSampleMixers;
            obj.allPhaseStatistics.frascatiStdSampleMixers = frascatiStdSampleMixers;
            obj.allPhaseStatistics.frascatiMeanStdSampleMixers = frascatiMeanStdSampleMixers;
            obj.allPhaseStatistics.frascatiMeanPulseMixers = frascatiMeanPulseMixers;
            obj.allPhaseStatistics.frascatiStdPulseMixers = frascatiStdPulseMixers;
            
            if (obj.usePets)
                obj.allPhaseStatistics.petsMeanSamplePhase = petsMeanSamplePhase;
                obj.allPhaseStatistics.petsMeanPulsePhase = petsMeanPulsePhase;
                obj.allPhaseStatistics.petsMeanStdSamplePhase = petsMeanStdSamplePhase;
                obj.allPhaseStatistics.petsStdSamplePhase = petsStdSamplePhase;
                obj.allPhaseStatistics.petsStdPulsePhase = petsStdPulsePhase;
            end
            
            obj.allPhaseStatistics.ct285SMeanSample = ct285SMeanSample;
            obj.allPhaseStatistics.ct285SStdSample = ct285SStdSample;
            obj.allPhaseStatistics.ct285SMeanStdSample = ct285SMeanStdSample;
            obj.allPhaseStatistics.ct285SMeanPulse = ct285SMeanPulse;
            obj.allPhaseStatistics.ct285SStdPulse = ct285SStdPulse;
            
            obj.allPhaseStatistics.ct285HMeanSample = ct285HMeanSample;
            obj.allPhaseStatistics.ct285HStdSample = ct285HStdSample;
            obj.allPhaseStatistics.ct285HMeanStdSample = ct285HMeanStdSample;
            obj.allPhaseStatistics.ct285HMeanPulse = ct285HMeanPulse;
            obj.allPhaseStatistics.ct285HStdPulse = ct285HStdPulse;
            
            obj.allPhaseStatistics.goodMeanCT285S = goodMeanCT285S;
            obj.allPhaseStatistics.goodMeanCT285H = goodMeanCT285H;
            obj.allPhaseStatistics.goodMeanFrascati1 = goodMeanFrascati1;
            obj.allPhaseStatistics.goodMeanFrascati2 = goodMeanFrascati2;
            obj.allPhaseStatistics.goodMeanFrascati3 = goodMeanFrascati3;
            if (obj.usePets)
                obj.allPhaseStatistics.goodMeanPets = goodMeanPets;
            end
            
            obj.allPhaseStatistics.corr12_meanPulsePhase = corr12_meanPulsePhase;
            obj.allPhaseStatistics.corr13_meanPulsePhase = corr13_meanPulsePhase;
            obj.allPhaseStatistics.corr23_meanPulsePhase = corr23_meanPulsePhase;
            if (obj.usePets)
                obj.allPhaseStatistics.corrPetsFrascati3_meanPulsePhase = corrPetsFrascati3_meanPulsePhase;           
            end
            
            obj.allPhaseStatistics.corrCT285S_Frascati1 = corrCT285S_Frascati1;
            obj.allPhaseStatistics.corrCT285S_Frascati2 = corrCT285S_Frascati2;
            obj.allPhaseStatistics.corrCT285S_Frascati3 = corrCT285S_Frascati3;
            obj.allPhaseStatistics.corrCT285H_Frascati1 = corrCT285H_Frascati1;
            obj.allPhaseStatistics.corrCT285H_Frascati2 = corrCT285H_Frascati2;
            obj.allPhaseStatistics.corrCT285H_Frascati3 = corrCT285H_Frascati3;
            
            obj.allPhaseStatistics.fit12_meanPulsePhase = fit12_meanPulsePhase;
            obj.allPhaseStatistics.fit13_meanPulsePhase = fit13_meanPulsePhase;
            obj.allPhaseStatistics.fit23_meanPulsePhase = fit23_meanPulsePhase;
            if (obj.usePets)
                obj.allPhaseStatistics.fitPetsFrascati3_meanPulsePhase = fitPetsFrascati3_meanPulsePhase;
            end
            
            obj.allPhaseStatistics.fitCT285S_Frascati1 = fitCT285S_Frascati1;
            obj.allPhaseStatistics.fitCT285S_Frascati2 = fitCT285S_Frascati2;
            obj.allPhaseStatistics.fitCT285S_Frascati3 = fitCT285S_Frascati3;
            obj.allPhaseStatistics.fitCT285H_Frascati1 = fitCT285H_Frascati1;
            obj.allPhaseStatistics.fitCT285H_Frascati2 = fitCT285H_Frascati2;
            obj.allPhaseStatistics.fitCT285H_Frascati3 = fitCT285H_Frascati3;


            
        end
    end
    
end

