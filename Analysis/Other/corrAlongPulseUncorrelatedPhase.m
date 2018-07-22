% Searching for correlation between a device and the residual uncorrelated
% phase component left over after an optimal FF correction.
close all;
%%
dataSetName= '20150715_1414_FF_Gain63_Gate185_350_Int_Odd';

dataDir = '/home/jack/PhaseFeedforward/Analysis/201507';
saveDir = '/home/jack/PhaseFeedforward/Analysis/201507/CorrAlongPulseResidualFFBPM';

bpmIndex = 22;%42=cc.930;%14=ct.608;%22=cr.155;
bpmSignal = 'H'; %'H','V','S'

samplesToUse = 550:700;%(pulseSampleRange{3}(1)+5):(pulseSampleRange{3}(end)-5);

%% load data
addpath('../');
addpath('../../');

savePlotDir = sprintf('%s/%s',saveDir,dataSetName);

load([dataDir '/' dataSetName '.mat']);

%% extract BPM signal and calculate residual phase after FF
if (strcmp(bpmSignal,'H'))
    meanBPM = meanBPMHAlongPulse;
elseif (strcmp(bpmSignal,'S'))
    meanBPM = meanBPMSAlongPulse;
elseif (strcmp(bpmSignal,'V'))
    meanBPM = meanBPMVAlongPulse;
end

gain = corrMeanMix2Mix3*(stdMeanPulsePhase(3)./stdMeanPulsePhase(2));

for mon=1:nMons
    phases(mon,:,:) = phases(mon,:,:) - subtractPhase(mon);
end

ffPhases = phases(3,:,:) - gain.*phases(2,:,:);
meanFFPhase = squeeze(nanmean(ffPhases));

%%

plotXRange = [samplesToUse(1) samplesToUse(end)];%[290 335];%[550 700];%[sampleRange(1)-20 sampleRange(end)+20]

normFF = meanFFPhase - nanmean(meanFFPhase(samplesToUse));
normFF = normFF./max(abs(normFF(samplesToUse)));

normBPM = meanBPM{bpmIndex} - nanmean(meanBPM{bpmIndex}(samplesToUse));
normBPM = normBPM./max(abs(normBPM(samplesToUse)));

corrBPMFF = nancorrcoef(normFF(samplesToUse),normBPM(samplesToUse));

if (corrBPMFF<0)
    normBPM = -normBPM;
    corrBPMFF = -corrBPMFF;
end

figure;
plot(normFF,'-o');
hold all;
plot(normBPM,'-x');
xlim(plotXRange);
title(sprintf('Residual Phase and %s%s: corr %.2f',bpmNames{bpmIndex},bpmSignal,corrBPMFF));
xlabel('Sample No.');
ylabel('Normalised Output [a.u.]');
legend('Residual Phase','BPM')
% savePlot(savePlotDir,sprintf('BPM%d%s',bpmIndex,bpmSignal));
