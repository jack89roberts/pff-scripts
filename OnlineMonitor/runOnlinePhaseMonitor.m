%% Create object
addpath('/nfs/cs-ccr-nfs6/vol29/share/ctf/data/ctfmod/MatLab/PhaseFeedforward/OnlineMonitor');
phaseMonObj = onlinePhaseMonitor();
phaseMonObj.noUpdateIfBeamDisabled = true;

% 11/12/14 Screwed up mean subtraction somehow whilst adding bad pulse
% removal etc.
phaseMonObj.subtractInitNPulses = true;

phaseMonObj.lengthPulseHistory = 300;

phaseMonObj.frascatiStartSamples = [620 620 620];
phaseMonObj.frascatiEndSamples = [680 680 680];
phaseMonObj.ct285StartSample = 575;%530;
phaseMonObj.ct285EndSample = 680;%700;
phaseMonObj.usePets = false;
phaseMonObj.petsStartSample = 125;
phaseMonObj.petsEndSample = 325;

phaseMonObj.plotSampleRange = 450:730;
%% Set start/end samples
% 
% frascatiNMonitors = phaseMonObj.frascatiNMonitors;
% frascatiStartSamples = NaN*ones(1,frascatiNMonitors);
% frascatiEndSamples = NaN*ones(1,frascatiNMonitors);
% 
% ctfSelector = 'SCT.USER.SETUP';
% 
% monMixerNames = {...
%     'CT.SCOPE01.CH02',...
%     'CT.SCOPE01.CH04',...
%     'CT.SCOPE01.CH06'...
% };

% monDiodeNames = {...
%     'CT.SCOPE01.CH01',...
%     'CT.SCOPE01.CH03',...
%     'CT.SCOPE01.CH05'...
% };
% 
% sampleSelectFigureHandle = figure();
% set(0,'CurrentFigure', sampleSelectFigureHandle);      
% 
% for mon=1:frascatiNMonitors    
%     subplot(1,2,1);
%     diodeName = [ monDiodeNames{mon} '/Acquisition#value'];
%     plot(matlabJapc.staticGetSignal(ctfSelector,diodeName));
%     title(sprintf('Mon %d Diode',mon));
%     xlabel('Sample No.');
%     ylabel('Output [au]');
%     
%     subplot(1,2,2);
%     mixerName = [ monMixerNames{mon} '/Acquisition#value'];
%     plot(matlabJapc.staticGetSignal(ctfSelector,mixerName));
%     title(sprintf('Mon %d Mixer',mon));
%     xlabel('Sample No.');
%     ylabel('Output [au]');
%     
%     fprintf(1, 'Monitor %d, Diode %s, Mixer %s:\n',mon,monDiodeNames{mon},monMixerNames{mon});
%     frascatiStartSamples(mon) = input('Start calibration at sample: ');
%     frascatiEndSamples(mon) = input('End calibration at sample: ');
%     fprintf(1,'------------------------------------\n');
% end
% 
% subplot(1,2,1)
% petsIName = 'CE.SCOPE03.CH01/Acquisition#value';
% plot(matlabJapc.staticGetSignal(ctfSelector,petsIName));
% title('PETS I');
% xlabel('Sample No.');
% ylabel('Output [au]');
% 
% subplot(1,2,2)
% petsQName = 'CE.SCOPE03.CH02/Acquisition#value';
% plot(matlabJapc.staticGetSignal(ctfSelector,petsQName));
% title('PETS Q');
% xlabel('Sample No.');
% ylabel('Output [au]');
% 
% fprintf(1, 'PETS, I = %s, Q = %s\n',petsIName,petsQName);
% petsStartSample = input('Start calibration at sample: ');
% petsEndSample = input('End calibration at sample: ');
% 
% close(sampleSelectFigureHandle);
% 
% phaseMonObj.frascatiStartSamples = frascatiStartSamples;
% phaseMonObj.frascatiEndSamples = frascatiEndSamples;
% phaseMonObj.petsStartSample = petsStartSample;
% phaseMonObj.petsEndSample = petsEndSample;


%% start
phaseMonObj.startAcquisition

%% pause
phaseMonObj.pauseAcquisition
%% resume
phaseMonObj.resumeAcquisition
%% stop
phaseMonObj.stopAcquisition
% or just close the plot windows... ;)