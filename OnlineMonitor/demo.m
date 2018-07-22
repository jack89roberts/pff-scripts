% create object
myobj = subscribeAndPlot('SCT.USER.SETUP', ...
    {'CT.SCOPE01.CH01/Acquisition#value',...
     'CT.SCOPE01.CH02/Acquisition#value',...
     'CT.SCOPE01.CH03/Acquisition#value'});

phShiftGetPropStr = 'Acquisition';
phShiftGetValStr = 'phase_Val';
% specify some options
myobj.stopIfNaNs=false;
myobj.xLims = [110,235];
myobj.yLims = 'auto';
myobj.persistentPlot = false;
xoffsets = [0 0 0];
yoffsets = [0 1 2];
xscalingfactors = [1 1 1];
yscalingfactors = [1 1 1];

%% start
myobj.startAcquisition
%% pause
myobj.pauseAcquisition
%% resume
myobj.resumeAcquisition
%% stop
myobj.stopAcquisition
% or just close the plot windows... ;)