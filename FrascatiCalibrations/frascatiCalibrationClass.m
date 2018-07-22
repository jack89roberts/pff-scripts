classdef frascatiCalibrationClass < handle
    %frascatiCalibrationClass 
    % Class used to run calibrations of the Frascati phase monitors.
    
    properties

        scanType = 'auto'; % 'auto' or 'manual'. Whether scan is using remote shifters (change them automatically) or mechanical shifters (prompt wait for them to be changed)
        
        scanNPoints; % no. of shifter settings in the scan
        scanNPulsesPerPoint; % no. of pulses collected for each shifter setting

        useMixerOverSqrtDiode = false; % if true calibrate using mixer/sqrt(diode), if false calibrate using mixer only
        overwriteCalibrations = false; % if true, new zero crossings will be set and reference calibration files will be overwritten 
        saveData = true; % if true all matlab data and files will be saved in the datadir directory
        
        calibrationComment = '';
        
        phShiftNames = { ...
           'CT.PHAS-PHM0370',...
           'CC.PHAS-PHM0251',...
           'CB.PHAS-PHM0160'...
        };


        monMixerNames = {...
            'CT.SCOPE01.CH02',...
            'CT.SCOPE01.CH04',...
            'CT.SCOPE01.CH06'...
        };

        monDiodeNames = {...
            'CT.SCOPE01.CH01',...
            'CT.SCOPE01.CH03',...
            'CT.SCOPE01.CH05'...
        };
        
        beamEnabledSignal = 'CL.STBPM0402S/Samples#samples';
         
        savePath = '/user/ctf3op/PhaseFeedforward/FrascatiCalibrations/';
        calTimeStamp;
    
        scanPhShiftValues;
        scanLOPhases;

        phaseShiftFreq = 3;
        
        nMons;
        nSamples;
        
        rawMixers;
        rawDiodes;
        rawMixerSensitivities;
        rawDiodeSensitivities;
        
        alignedMixers
        alignedDiodes
        alignedMixerSqrtDiode;
        
        calSampleRange;

        zeroCrossings;
        calibrationFactors;
        calRSquare;
        calFitConf;
        newMixerSensitivities;
        newDiodeSensitivities;
        
        dataMonitor; % JAPC MatLab monitor
        signals;
    end
    properties (Access=private)
        
        cycleName = 'SCT.USER.SETUP';

        
        initPhShiftValues; % initial phase shifter settings
        scanStartPhase = 90.0; % initial value for phase shifters in the scan
        scanEndPhase = 180.0; % final value for phase shifters in the scan
        
        phShiftGetStr = 'Acquisition#phase_Val';
        phShiftSetStr = 'Setting#phase';
        phMonAcqStr = 'Acquisition#value';
        phMonSensStr = 'Acquisition#sensitivity';

        scanShifterIndex = 1;
        scanPulseIndex = 0;
        
        beamIsEnabled = 1;
        nPulsesNoBeam = 0;
                
    end
    methods
        % constructor
        function obj=frascatiCalibrationClass(scanType, scanNPulsesPerPoint, scanNPoints)
            fprintf('Creating Frascati Calibration object... ');
            
            if(strcmp(scanType,'auto') || strcmp(scanType,'manual'))
                obj.scanType = scanType;
                obj.scanNPoints = Inf;
            else
                error('scanType of %s is not valid. Must be "auto" or "manual"',scanType);
            end
            
            if nargin<2
                error('Not enough input arguments');
            elseif (strcmp(scanType,'manual') && nargin<3)
                error('Not enough input arguments')
            end
            
            obj.scanNPoints = scanNPoints;
            obj.scanNPulsesPerPoint = scanNPulsesPerPoint;
                        
            obj.nMons = length(obj.monMixerNames);
            obj.signals = cell(1,(obj.nMons.*2)+1);
            i=1;
            for mon=1:obj.nMons
                obj.signals{i} = [obj.monDiodeNames{mon} '/' obj.phMonAcqStr];
                i = i+1;
                obj.signals{i} = [obj.monMixerNames{mon} '/' obj.phMonAcqStr];
                i = i+1;
            end
            obj.signals{i} = obj.beamEnabledSignal;
            
            obj.dataMonitor = matlabJapcMonitor(obj.cycleName, obj.signals, ...
                @(e)dataReceived(obj,e),'');
            
            obj.dataMonitor.useFastStrategy(1);
            
            fprintf('done!\n');
        end
        
        function startScan(obj)
            obj.initialiseScan();
            fprintf('Starting acquisition... ');
            obj.dataMonitor.start();
            fprintf('done!\n');
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
        
        function processScanResults(obj)
            % process signals
            obj.scanLOPhases = 4*obj.scanPhShiftValues - min(4*obj.scanPhShiftValues);
            
            fprintf('Processing signals...\n');
            [obj.alignedMixers, obj.alignedDiodes, obj.calSampleRange] = ...
                frascatiCalProcessSignals(...
                    obj.rawMixers,...
                    obj.rawDiodes,...
                    obj.rawMixerSensitivities,...
                    obj.rawDiodeSensitivities);
            fprintf('Finished processing signals!\n');   
             
            % fit data
            fprintf('Fitting data...\n');
            [obj.calibrationFactors, obj.zeroCrossings, obj.calRSquare, obj.calFitConf] = ...
                frascatiCalCalculateCalibrations(...
                    obj.scanPhShiftValues,...
                    obj.alignedMixers,...
                    obj.alignedDiodes,...
                    obj.calSampleRange,...
                    obj.useMixerOverSqrtDiode,...
                    obj.phaseShiftFreq);
            fprintf('Finished fitting data!\n');
            
            % calculate optimal sensitivities - to be implemented
            
            % save data
            fprintf('Saving data... ')
            
            obj.saveCalData();
            
            zeroCrossName = sprintf('%s%s%s',obj.savePath,'frascatiZeroCrossings_',obj.calTimeStamp);
            zeroCrossFileBackup = fopen(zeroCrossName,'w');
            calibConstName = sprintf('%s%s%s',obj.savePath,'frascatiCalibrationConstants_',obj.calTimeStamp);
            calibConstFileBackup = fopen(calibConstName,'w');
            for mon=1:obj.nMons 
                fprintf(zeroCrossFileBackup,'%s, %.2f\n',...
                    obj.phShiftNames{mon},...
                    obj.zeroCrossings(mon)...
                );
                fprintf(calibConstFileBackup,'%s, %f, %f, %f, %f\n',...
                    obj.monMixerNames{mon},...
                    obj.calibrationFactors(mon,1),...
                    obj.calibrationFactors(mon,2),...
                    obj.calibrationFactors(mon,3),...
                    obj.calibrationFactors(mon,4)...
                );
            end
            fclose(zeroCrossFileBackup);
            
            fprintf(calibConstFileBackup, 'useMixerOverSqrtDiode, %d', obj.useMixerOverSqrtDiode);
            fclose(calibConstFileBackup);

            
            
            fprintf('done!\n');
            
            % Set values in machine if requested
            
            setCalInMachine = input('Save these values as reference (0 or 1)?');
            
            if (setCalInMachine)
                fprintf('Setting new shifter values and saving reference files... ')
                % files for saving reference calibrations
                zeroCrossFile = fopen('frascatiZeroCrossings','w');
                calibConstFile = fopen('frascatiCalibrationConstants','w');
                
                refNameFile = fopen('frascatiRefCalName','w');
                fprintf(refNameFile,'%s',obj.calTimeStamp);
                fclose(refNameFile);
                
                for mon=1:obj.nMons
                    % always want to set the zero crossing to the falling slope so check
                    % this here
                    dydxAtZero =...
                        obj.calibrationFactors(mon,1)*obj.calibrationFactors(mon,2)*...
                        cos(obj.calibrationFactors(mon,2)*obj.zeroCrossings(mon) + obj.calibrationFactors(mon,3)); % Ab*cos(bx+c)
                    if (dydxAtZero > 0)
                        zeroCrossingToSet = obj.zeroCrossings(mon) + (180/(12/obj.phaseShiftFreq)); % divide by 4 because of 3 GHz phase shifter vs. 12 GHz mixer
                    else
                        zeroCrossingToSet = obj.zeroCrossings(mon);
                    end
                    
                    if (strcmp(obj.scanType,'auto'))
                        obj.setPhaseShifter(mon,zeroCrossingToSet);
                    end
                    
                    fprintf(zeroCrossFile,'%s, %.2f\n',...
                        obj.phShiftNames{mon},...
                        zeroCrossingToSet...
                    );
                    fprintf(calibConstFile,'%s, %f, %f, %f, %f\n',...
                        obj.monMixerNames{mon}, obj.calibrationFactors(mon,1),...
                        obj.calibrationFactors(mon,2),...
                        obj.calibrationFactors(mon,3),...
                        obj.calibrationFactors(mon,4)...
                    );

                end
                
                fclose(zeroCrossFile);
                
                fprintf(calibConstFile, 'useMixerOverSqrtDiode, %d', obj.useMixerOverSqrtDiode);
                fclose(calibConstFile);
                
            end
            fprintf('done!\n')

           

        end
                
        function saveCalData(obj)
            varsToSave = {... 
                'calibrationComment',...  
                'useMixerOverSqrtDiode',...
                'phShiftNames',...
                'monMixerNames',...
                'monDiodeNames',...
                'calTimeStamp',...
                'scanPhShiftValues',...
                'scanLOPhases',...
                'nMons',...
                'nSamples',...
                'rawMixers',...
                'rawDiodes',...
                'rawMixerSensitivities',...
                'rawDiodeSensitivities',...
                'alignedMixers',...
                'alignedDiodes',...
                'alignedMixerSqrtDiode',...
                'calSampleRange',...
                'zeroCrossings',...
                'calibrationFactors',...
                'calRSquare',...
                'calFitConf',...
                'newMixerSensitivities',...
                'newDiodeSensitivities',...
                'phaseShiftFreq',...
                'scanType'...
            };
            saveFilePath = [obj.savePath 'frascatiCalibration_' obj.calTimeStamp];
            
            eval(sprintf('%s = obj.%s;',varsToSave{1},varsToSave{1}));
            save(saveFilePath,varsToSave{1});
            for i=2:length(varsToSave)
                eval(sprintf('%s = obj.%s;',varsToSave{i},varsToSave{i}));
                save(saveFilePath,varsToSave{i},'-append');
            end
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
            transmission = matlabDataAndSignalsHelper.simpleExtractSingleSignal(dataStruct,obj.signals{end});
            
            % DEBUG - run even if beam disabled
            %transmission = NaN;
            %warning('DEBUG: Transmission set to NaN to force scan to continue. Search for DEBUG in fracatiCalibrationClass.m and comment appropriate lines to stop this.');
            
            if (min(transmission) > -2) % NB: If min(transmission) returns NaN, beam is assumed to be present.
                obj.nPulsesNoBeam = obj.nPulsesNoBeam + 1;
                
                if (obj.beamIsEnabled || mod(obj.nPulsesNoBeam,20)==0)
                    obj.beamIsEnabled = 0;
                    fprintf('Beam Disabled. Waiting to continue scan...\n');
                end
                
                return;
            
            else
                if (~obj.beamIsEnabled)
                    obj.beamIsEnabled = 1;
                    obj.nPulsesNoBeam = 0;
                    if (obj.scanPulseIndex > 1) % repeat last measurement in case beam was disabled then
                        obj.scanPulseIndex = obj.scanPulseIndex - 1;
                    end
                    
                    fprintf('Beam Enabled. Waiting for a few pulses before continuing scan...\n');
                    obj.pauseScan();
                    pause(10);
                    obj.resumeScan();
                    
                    return;
                end
            end
            
            % check acquisition ok
            tmpDiodes = NaN(obj.nMons,obj.nSamples);
            tmpMixers = NaN(obj.nMons,obj.nSamples);
            i=1;
            for mon=1:obj.nMons
                tmpDiodes(mon,:) = matlabDataAndSignalsHelper.simpleExtractSingleSignal(dataStruct,obj.signals{i});
                if (isnan(max(tmpDiodes(mon,:))))
                    matlabDataAndSignalsHelper.simpleExtractSingleSignalError(dataStruct,obj.signals{i})
                    fprintf('Bad acquisition. Retrying...\n');
                    return;
                end
                i=i+1;
                
                tmpMixers(mon,:) = matlabDataAndSignalsHelper.simpleExtractSingleSignal(dataStruct,obj.signals{i});
                if (isnan(max(tmpMixers(mon,:))))
                    matlabDataAndSignalsHelper.simpleExtractSingleSignalError(dataStruct,obj.signals{i})
                    fprintf('Bad acquisition. Retrying...\n');
                    return;
                end
                i=i+1;
                
                
            end
            
            obj.scanPulseIndex = obj.scanPulseIndex + 1;
            
            fprintf('LO Phase %.2f degrees. Pulse %d of %d.\n',...
                obj.scanPhShiftValues(obj.scanShifterIndex),...
                obj.scanPulseIndex,...
                obj.scanNPulsesPerPoint...
            );

            % save data to arrays
            i=1;
            for mon=1:obj.nMons
                
                obj.rawDiodes(mon,obj.scanShifterIndex,obj.scanPulseIndex,:) =...
                    matlabDataAndSignalsHelper.simpleExtractSingleSignal(dataStruct,obj.signals{i});
                
                i = i+1;
                
                obj.rawMixers(mon,obj.scanShifterIndex,obj.scanPulseIndex,:) =...
                    matlabDataAndSignalsHelper.simpleExtractSingleSignal(dataStruct,obj.signals{i}); 
                
                i = i+1;
            end

            % Change phase shifters after acquiring scanNPulsesPerPoint.
            % Pause acquisition whilst changing shifters
            if (obj.scanPulseIndex == obj.scanNPulsesPerPoint)
                obj.pauseScan();
                
                if(strcmp(obj.scanType,'auto')) % Change remote shifters
                    if (obj.scanShifterIndex < obj.scanNPoints)
                        obj.scanShifterIndex = obj.scanShifterIndex + 1;
                        obj.scanPulseIndex = 0;

                        fprintf('Changing shifters... ');
                        obj.setAllPhaseShifters(obj.scanPhShiftValues(obj.scanShifterIndex));
                        pause(0.5);
                        obj.setAllPhaseShifters(obj.scanPhShiftValues(obj.scanShifterIndex));
                        pause(0.5);
                        obj.setAllPhaseShifters(obj.scanPhShiftValues(obj.scanShifterIndex));
                        pause(3);
                        fprintf('done!\n');

                        if (obj.saveData)
                            obj.saveCalData();
                        end

                        obj.resumeScan();

                    else
                        obj.stopScan();

                        if (obj.saveData)
                            obj.saveCalData();
                        end

                        fprintf('Scan complete.\n');
                        obj.processScanResults();
                    end
                    
                else % Manual scan, prompt to change shifters
                    
                    newShifterPhase = obj.changeManualShifter();

                    if(strcmp(newShifterPhase,'stop')) % end the scan
                        obj.stopScan();

                        obj.scanNPoints = obj.scanShifterIndex;
                        
                        if (obj.saveData)
                            obj.saveCalData();
                        end

                        fprintf('Scan complete.\n');
                        obj.processScanResults();
                        
                    else % continue the scan
                        obj.scanShifterIndex = obj.scanShifterIndex + 1;
                        obj.scanPulseIndex = 0;                    
                        
                        obj.scanPhShiftValues(obj.scanShifterIndex) = newShifterPhase;
                        
                        if (obj.saveData)
                            obj.saveCalData();
                        end

                        obj.resumeScan();

                    end

    
                end
                
            end
        end
        
        % initialise arrays, get and set initial phase shifter settings etc. before
        % starting scan.
        function initialiseScan(obj)
            fprintf('Initialising scan parameters...\n');

            obj.scanShifterIndex = 1;
            obj.scanPulseIndex = 0;
            obj.beamIsEnabled = 1;
            obj.nPulsesNoBeam = 0;

            if (strcmp(obj.scanType,'auto'))
                obj.scanPhShiftValues = linspace(obj.scanStartPhase,obj.scanEndPhase,obj.scanNPoints);
            else
                obj.scanPhShiftValues = [];
            end
            obj.initPhShiftValues = NaN(1,obj.nMons);            
            obj.rawMixerSensitivities = NaN(1,obj.nMons);
            obj.rawDiodeSensitivities = NaN(1,obj.nMons);
            
            for mon = 1:obj.nMons
                
                if (strcmp(obj.scanType,'auto'))
                    obj.initPhShiftValues(mon) =...
                        matlabJapc.staticGetSignal(...
                            obj.cycleName,...
                            sprintf('%s/%s',obj.phShiftNames{mon},obj.phShiftGetStr)...
                        );
                end
                
                obj.rawMixerSensitivities(mon) =...
                    matlabJapc.staticGetSignal(...
                        obj.cycleName,...
                        sprintf('%s/%s',obj.monMixerNames{mon},obj.phMonSensStr)...
                );
            
                obj.rawDiodeSensitivities(mon) =...
                    matlabJapc.staticGetSignal(...
                    obj.cycleName,...
                    sprintf('%s/%s',obj.monDiodeNames{mon},obj.phMonSensStr)...
                );
                     
            end
            
            tmpSignal =...
                matlabJapc.staticGetSignal(...
                    obj.cycleName,...
                    sprintf('%s/%s',obj.monMixerNames{mon},obj.phMonAcqStr)...
                );
            obj.nSamples = length(tmpSignal);
            
            if (strcmp(obj.scanType,'auto'))
                obj.rawMixers = NaN(obj.nMons,obj.scanNPoints,obj.scanNPulsesPerPoint, obj.nSamples);
                obj.rawDiodes = NaN(obj.nMons,obj.scanNPoints,obj.scanNPulsesPerPoint, obj.nSamples);
            else % just define arrays with 2 scan points initially, then grow array during scan
                obj.rawMixers = NaN(obj.nMons,2,obj.scanNPulsesPerPoint, obj.nSamples);
                obj.rawDiodes = NaN(obj.nMons,2,obj.scanNPulsesPerPoint, obj.nSamples);                
            end
            dV = datevec(now);
            obj.calTimeStamp = sprintf('%02d%02d%02d_%02d%02d',dV(1),dV(2),dV(3),dV(4),dV(5));
            
            if (strcmp(obj.scanType,'auto')) % set remote shifters to first value
                obj.setAllPhaseShifters(obj.scanPhShiftValues(1));
                pause(0.5);
                obj.setAllPhaseShifters(obj.scanPhShiftValues(1));
            else % prompt for first value of mechanical shifters
                newShifterPhase = obj.changeManualShifter();
                if (strcmp(newShifterPhase,'stop'))
                    error('Tried to stop scan before it started.')
                else
                    obj.scanPhShiftValues(1) = newShifterPhase;
                end
            end
            
            fprintf('Parameters initialised!\n');
        end
        
        % Set phase shifter values 
        function setPhaseShifter(obj,index,value)          
            if (isnan(value))
                warning('Phase Shifter %d (%s) not set: NaN value.',index,obj.phShiftNames{index});
                return;
            else
                matlabJapc.staticSetSignal(...
                    obj.cycleName,...
                    sprintf('%s/%s',obj.phShiftNames{index},obj.phShiftSetStr),...
                    value...
                );
            end
        end
        function setAllPhaseShifters(obj,value)
            for mon=1:obj.nMons
                obj.setPhaseShifter(mon,value);
            end
        end
        
        function newShifterPhase = changeManualShifter(obj) % Method called when shifter needs to be changed in a manual scan
            shifterIsSet = -1;
            while (shifterIsSet == -1)
                newShifterPhase = obj.promptForNewShifterVal();
                shifterIsSet = obj.isShifterSet(newShifterPhase);
            end
        end
        
        function newShifterPhase = promptForNewShifterVal(obj) % asks user for new phase shifter phase (manual scan)
            newShifterPhase = NaN;
            while(isnan(newShifterPhase))
                newShifterPhase = input('New shifter phase [or "stop" to end scan]: ','s');

                if(strcmp(newShifterPhase,'stop'))
                    break;
                else
                    newShifterPhase = str2double(newShifterPhase);
                end
            end 
        end
        
        function shifterIsSet = isShifterSet(obj, newShifterPhase) % asks user to verify that phase shifter setting is correct (manual scan)
            shifterIsSet = 0;
            while(shifterIsSet ~= 1 && shifterIsSet ~= -1)
                if (strcmp(newShifterPhase,'stop'))
                    try
                        shifterIsSet = input(sprintf('Stop the scan [1 to end, -1 to continue]? '));
                    catch
                        shifterIsSet = 0;
                    end
                else
                    try
                        shifterIsSet = input(sprintf('Shifter value %.2f set [1 to continue, -1 to change value]? ',newShifterPhase));
                    catch
                        shifterIsSet = 0;
                    end
                end
            end
        end
    end
    
end

