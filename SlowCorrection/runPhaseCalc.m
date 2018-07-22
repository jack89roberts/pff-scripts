close all; clear all;
%%
corrRatio = -1.2; % ratio (765/465) between corrector strengths needed for orbit closure.
nAvg = 20; % no. of pulses acquired for each point
corr465Offsets = -2:0.5:2; % offsets to apply to corrector 465

%%
phaseObj = phase_Class(corrRatio);
phaseObj.nAvg = nAvg;
phaseObj.corr465Offsets = corr465Offsets;
phaseObj.startScan();