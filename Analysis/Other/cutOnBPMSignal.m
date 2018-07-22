% Compares results between all pulses and pulses with a +/-0.5 sigma cut
% around the mean on some BPM signal.

%dataSetName= '20150707_2037_FFGain_63_Interleaved_Odd';
%dataSetName= '20150707_2138_FFGain_43_Interleaved_Odd';
%dataSetName= '20150707_2052_FFGain_53_Interleaved_Odd';
%dataSetName= '20150707_2102_FFGain_43_Interleaved_Even';
dataSetName= '20150715_1414_FF_Gain63_Gate185_350_Int_Odd';

dataDir = '/home/jack/PhaseFeedforward/Analysis/201507';
saveDir = '/home/jack/PhaseFeedforward/Analysis/201507/CutOnSeptum';

plotXRange = [550 700];%[290 335];%[550 700];%[sampleRange(1)-20 sampleRange(end)+20]

bpmIndex = 14;%42=cc.930;%14=ct.608;%22=cr.155;
bpmSignal = 'H'; %'H','V','S'
%% load data
addpath('../');
addpath('../../');

savePlotDir = sprintf('%s/%s_BPM%s%d',saveDir,dataSetName,bpmSignal,bpmIndex);

load([dataDir '/' dataSetName '.mat']);

%%
if (strcmp(bpmSignal,'H'))
    meanBPM = meanBPMH;
elseif (strcmp(bpmSignal,'S'))
    meanBPM = meanBPMS;
elseif (strcmp(bpmSignal,'V'))
    meanBPM = meanBPMV;
end

meanCut = nanmean(meanBPM{bpmIndex});
stdCut = nanstd(meanBPM{bpmIndex});
highCut = meanCut+0.5*stdCut;
lowCut = meanCut-0.5*stdCut;

figure;
hist(meanBPM{bpmIndex},20)
hold all;
title(sprintf('%s%s',bpmNames{bpmIndex},bpmSignal))
xlabel('position [mm]')
ylabel('nPulses');
plot([lowCut lowCut],get(gca,'YLim'),'k','LineWidth',2);
plot([highCut highCut],get(gca,'YLim'),'k','LineWidth',2);
savePlot(savePlotDir,'histBPM');

cutPhases = phases;
for i=1:nPulses
    if (meanBPM{bpmIndex}(i)>highCut || meanBPM{bpmIndex}(i)<lowCut)
        cutPhases(:,i,:) = NaN;
    end
end
for mon=1:nMons
    cutPhases(mon,:,:) = cutPhases(mon,:,:)-subtractPhase(mon);
end

cutMeanPulsePhase = nanmean(cutPhases(:,:,sampleRange),3);
cutStdAlongPulse = nanstd(cutPhases,[],2);
cutMeanAlongPulse = nanmean(cutPhases,2);

figure;
plot(cutStdAlongPulse(3,:));
hold all;
plot(stdPhaseAlongPulse(3,:));
legend('Cut','All');
xlim(plotXRange);
xlabel('sample no.')
ylabel('jitter (Degrees)');
savePlot(savePlotDir,'stdAlong');

figure;
plot(cutMeanAlongPulse(3,:));
hold all;
plot(meanPhaseAlongPulse(3,:));
legend('Cut','All');
xlim(plotXRange);
xlabel('sample no.')
ylabel('phase (Degrees)');
savePlot(savePlotDir,'meanAlong');

figure;
scatter(cutMeanPulsePhase(2,:),cutMeanPulsePhase(3,:));
title(sprintf('cut, corr=%.2d',nancorrcoef(cutMeanPulsePhase(2,:),cutMeanPulsePhase(3,:))));
xlabel('phase (degrees)')
ylabel('phase (Degrees)');

savePlot(savePlotDir,'scatterAll');

figure;
scatter(meanPulsePhase(2,:),meanPulsePhase(3,:));
title(sprintf('all, corr=%.2d',corrMeanMix2Mix3));
xlabel('phase (degrees)')
ylabel('phase (Degrees)');
savePlot(savePlotDir,'scatterCut');
