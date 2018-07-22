%% INPUTS
scanType = 'manual'; % 'auto' (with remote shifters) or 'manual' (mechanical)
scanNPulsesPerPoint = 3; % no. pulses acquired for each shifter setting
scanNPoints = 10; % no. points in scan, only needed for auto calibrations. Has no effect in manual calibrations.
phaseShiftFreq = 4;%0.866; % Frequency of readings on phase shifter (4GHz for new mechanical shifters, 1GHz for old Mechanical shifters, 3GHz for remote shifters, 0.866 for adjusters).
saveData = true;
calibrationComment = 'All 3 new shifters';
useMixerOverSqrtDiode = false;

%%
frascatiCalObj = frascatiCalibrationClass(scanType, scanNPulsesPerPoint, scanNPoints);
frascatiCalObj.saveData = saveData;
frascatiCalObj.calibrationComment = calibrationComment;
frascatiCalObj.useMixerOverSqrtDiode = useMixerOverSqrtDiode;
frascatiCalObj.phaseShiftFreq = phaseShiftFreq;

warning('off','backtrace') % simplifies warning messages (may help to turn back on if debugging needed)

frascatiCalObj.startScan();