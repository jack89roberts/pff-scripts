classdef frascatiCalibrationClass < handle
    %frascatiCalibrationClass 
    % Class used to run calibrations of the Frascati phase monitors.
    
    properties
           
        scanNPoints; % no. of shifter settings in the scan
        scanNPulsesPerPoint; % no. of pulses collected for each shifter setting

        useMixerOverSqrtDiode = false; % if true calibrate using mixer/sqrt(diode), if false calibrate using mixer only
        overwriteCalibrations = false; % if true, new zero crossings will be set and reference calibration files will be overwritten 
        saveData = true; % if true all matlab data and files will be saved in the datadir directory
        
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
         
        savePath = 'data/frascatiCalibrations/';
        fileName;
    
        scanPhShiftValues;
        scanLOPhases;

        nMons;
        nSamples;
        
        rawMixers;
        rawDiodes;
        rawMixerSensitivities;
        rawDiodeSensitivities;
        
        alignedMixers
        alignedDiodes
        alignedMixerSqrtDiode;
        
        calStartSample;
        calEndSample;

        zeroCrossings;
        calibrationFactors;
        newMixerSensitivities;
        newDiodeSensitivities;

 
    end
    properties (Access=private)
        dataMonitor; % JAPC MatLab monitor
        cycleName = 'SCT.USER.SETUP';
        signals;
        
        scanStartPhase = 90.0; % initial value for phase shifters in the scan
        scanEndPhase = 180.0; % final value for phase shifters in the scan
        
        phShiftGetPropStr = 'Acquisition';
        phShiftGetValStr = 'phase_Val';
        phShiftSetPropStr = 'Setting';
        phShiftSetValStr = 'phase';
        phAcqPropStr = 'Acquisition';
        phAcqValStr = 'value';

    end
    methods
        % constructor
        function obj=frascatiCalibrationClass(scanNPoints, scanNPulsesPerPoint)
            obj.scanNPoints = scanNPoints;
            obj.scanNPulsesPerPoint = scanNPulsesPerPoint;
            
            
            obj.dataMonitor = matlabJapcMonitor(obj.cycleName, signals, ...
                @(e)dataReceived(obj,e),'');
        end
        
        function stopAcquisition(obj)
            obj.dataMonitor.stop();
        end
        function startAcquisition(obj)
            obj.dataMonitor.start();
        end
        function pauseAcquisition(obj)
            obj.dataMonitor.pauseOn();
        end
        function resumeAcquisition(obj)
            obj.dataMonitor.pauseOff();
        end
        function delete(obj)
            obj.stopAcquisition();
            obj.dataMonitor.delete();
        end

    end
    methods(Access=private)
        function dataReceived(obj, dataStruct)

        end
  
    end
    
end

