close all; clear all;
%%
nAvg = 15; % no. of pulses acquired for each point
corr465Offsets = [-2.0 -1.0 1.0 2.0]; % offsets to apply to corrector 465
minCorrRatio = -2.0; % min765 offset = min465 offset * minCorr765Ratio
maxCorrRatio = -0.5; % max765 offset = max465 offset * maxCorr765Ratio
n765Offsets = 10; % no. of points in scan (for each setting of 465)

%%
closureObj = closure_Class();
closureObj.nAvg=nAvg; % no. of pulses acquired for each point
closureObj.corr465Offsets=corr465Offsets; % offsets to apply to corrector 465
closureObj.minCorrRatio=minCorrRatio; % min765 offset = min465 offset * minCorr765Ratio
closureObj.maxCorrRatio=maxCorrRatio; % max765 offset = max465 offset * maxCorr765Ratio
closureObj.n765Offsets=n765Offsets; % no. of points in scan (for each setting of 465)

closureObj.startScan();