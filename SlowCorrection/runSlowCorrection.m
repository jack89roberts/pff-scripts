close all; clear all;
%%
gain = -0.5; % default: 1.0. if calculated phase offset is x degrees, correction attempted is gain*x.
corrRatio = -1.2;  % default: [], ratio (765/465) of corrector strengths for orbit closure. If empty loads from ref files.
phaseShiftPerAmp = 15.21; % default: [], phase shift per amp applied to 465. If empty loads from ref files.

nAvg = 20; % default: 20, correction updates every nAvg pusles
nHistory = 1000; % default: 1000, no. of pulses to keep in history
saveEvery = 200; % default: 200, file saved with current data every saveEvery pulses

%%
corrObj = slowCorr_Class();

corrObj.gain = gain;
corrObj.corrRatio = corrRatio;
corrObj.phaseShiftPerAmp = phaseShiftPerAmp;
corrObj.nAvg = nAvg;
corrObj.nHistory = nHistory;
corrObj.saveEvery = saveEvery;

corrObj.startCorrection();