%% alignmentCheck.m
% Checks the sample delay that gives the best sample-by-sample
% correlation and which sample delay gives the flattest simulated corrected
% phase using a processed data file.
close all; clearvars;

ffOff = load('/home/jack/PhaseFeedforward/Analysis/201511/20151120_1538_Gain-800_R56_0.1_Interleaaved_Odd.mat');
ffOn = load('/home/jack/PhaseFeedforward/Analysis/201511/20151120_1538_Gain-800_R56_0.1_Interleaaved_Even.mat');
saveDir = '/home/jack/PhaseFeedforward/Analysis/201511/Plots/20151120_1538_Gain-800_R56_0.1/Simulation';

sampDelay = -20:20;

%% extract variables to use from data and plot initial alignment

nDelay = length(sampDelay);
colors = varycolor(nDelay);

mon2All = squeeze(ffOff.phases(2,:,:));
mon3All = squeeze(ffOff.phases(3,:,:));
mon2Mean = ffOff.meanPhaseAlongPulse(2,:);
mon3Mean = ffOff.meanPhaseAlongPulse(3,:);
mon2Std = ffOff.stdPhaseAlongPulse(2,:);
mon3Std = ffOff.stdPhaseAlongPulse(3,:);
diode2Mean = ffOff.meanDiodeAlongPulse(2,:);
diode3Mean = ffOff.meanDiodeAlongPulse(3,:);
sampleRange = ffOff.sampleRange;
nSamples = ffOff.nSamples;
xLim = [ffOff.pulseSampleRange{3}(1)-30 ffOff.pulseSampleRange{3}(end)+30];

figure;
plot(mon2Mean,'LineWidth',2);
hold all;
plot(mon3Mean,'LineWidth',2);
title({'Original Alignment' 'Mean Phase Along Pulse'});
xlim(xLim);
xlabel('Sample No.')
ylabel('Phase [degrees]');
legend('Mon2','Mon3')
format_plots;
savePlot(saveDir,'originalPhaseAlong');

figure;
plot(diode2Mean,'LineWidth',2);
hold all;
plot(diode3Mean,'LineWidth',2);
title({'Original Alignment' 'Mean Diode'});
xlim(xLim);
xlabel('Sample No.')
ylabel('Output [V]');
legend('Mon2','Mon3')
format_plots;
savePlot(saveDir,'originalDiode');

%% Find delay that gives best sample correlation

corrSamples = NaN(nDelay,nSamples);
corrSamples_err = NaN(nDelay,nSamples);

% figure;

for i=1:nDelay
    delayedMon3 = delaySignal(mon3All,sampDelay(i));
    
    for s=1:nSamples    
        [corrSamples(i,s),corrSamples_err(i,s)] = nancorrcoef(mon2All(:,s),delayedMon3(:,s));
    end
    
%     plot(corrSamples(i,:),'Color',colors(i,:))
%     hold all;
end

meanCorrSamples = nanmean(corrSamples(:,sampleRange),2);

figure;
plot(sampDelay,meanCorrSamples,'LineWidth',2);
title('Mean Sample Correlation')
xlabel('Delay Added to Mon3 [samples]')
ylabel('Correlation')
format_plots;
savePlot(saveDir,'sampCorrelationVsDelay');

[maxCorr,bestDelayInd] = max(meanCorrSamples);
bestDelay = sampDelay(bestDelayInd);

figure;
plot(mon2Mean,'LineWidth',2);
hold all;
plot(delaySignal(mon3Mean,bestDelay),'LineWidth',2);
title({'Mean Phase Along Pulse, Max Sample Correlation' sprintf('Mon3 Delay: %d samples',bestDelay)});
xlim(xLim);
xlabel('Sample No.')
ylabel('Phase [degrees]');
legend('Mon2','Mon3')
format_plots;
savePlot(saveDir,'maxCorrelPhaseAlong');

figure;
plot(diode2Mean,'LineWidth',2);
hold all;
plot(delaySignal(diode3Mean,bestDelay),'LineWidth',2);
title({'Mean Diode, Max Sample Correlation' sprintf('Mon3 Delay: %d samples',bestDelay)});
xlim(xLim);
xlabel('Sample No.')
ylabel('Output [V]');
legend('Mon2','Mon3')
format_plots;
savePlot(saveDir,'maxCorrelDiode');

%% Find delay that gives flattest PFF correction

data = load('/home/jack/PhaseFeedforward/Analysis/201511/20151120_1538_Gain-800_R56_0.1_Interleaaved_Odd.mat');

simFFResults = cell(1,nDelay);
flatness = NaN(1,nDelay);
flatness_err = NaN(1,nDelay);
stdAlong = NaN(1,nDelay);
stdMean = NaN(1,nDelay);

for i=1:nDelay
    delayedData = ffOff;
    delayedData.phases(3,:,:) = delaySignal(mon3All,sampDelay(i));
    
    simFFResults{i} = getSimulatedFF(delayedData);
    
    [flatness(i),~,flatness_err(i),~] = nanMeanStdErr(simFFResults{i}.flatnessSimFF);
    stdAlong(i) = nanmean(simFFResults{i}.stdSimFFAlongPulse(sampleRange));
    stdMean(i) = simFFResults{i}.stdSimFFPhase;
end

[realFlatness,~,realFlatness_err,~] = nanMeanStdErr(ffOn.pulsePhaseFlatness(3,:));
realStdAlong = ffOn.meanStdPhaseAlongPulse(3);
realStdMean = ffOn.stdMeanPulsePhase(3);

figure;
plot(sampDelay,flatness,'LineWidth',2)
title('Corrected Phase Flatness')
xlabel('Delay Added to Mon3 [samples]')
ylabel('Flatness [degrees]')
hold all;
plot([sampDelay(1) sampDelay(end)],[realFlatness realFlatness],'Color','k','LineWidth',2)
legend('Simulated','Real')
format_plots;
savePlot(saveDir,'flatnessVsDelay');

figure;
plot(sampDelay,stdAlong,'LineWidth',2)
title('Corrected Jitter Along Pulse')
xlabel('Delay Added to Mon3 [samples]')
ylabel('Jitter [degrees]')
hold all;
plot([sampDelay(1) sampDelay(end)],[realStdAlong realStdAlong],'Color','k','LineWidth',2)
legend('Simulated','Real')
legend('Simulated','Real')
format_plots;
savePlot(saveDir,'jitterAlongVsDelay');

figure;
plot(sampDelay,stdMean,'LineWidth',2)
title('Corrected Mean Phase Jitter')
xlabel('Delay Added to Mon3 [samples]')
ylabel('Jitter [degrees]')
hold all;
plot([sampDelay(1) sampDelay(end)],[realStdMean realStdMean],'Color','k','LineWidth',2)
format_plots;
savePlot(saveDir,'meanJitterVsDelay');

[minFlatness,bestFlatnessInd] = min(flatness);
bestFlatnessDelay = sampDelay(bestFlatnessInd);

figure;
plot(mon2Mean,'LineWidth',2);
hold all;
plot(delaySignal(mon3Mean,bestFlatnessDelay),'LineWidth',2);
plot(simFFResults{bestFlatnessInd}.meanSimFFAlongPulse,'LineWidth',2); 
plot(ffOn.meanPhaseAlongPulse(3,:),'LineWidth',2)
title({'Best Sim Flatness, Mean Phase Along Pulse'  sprintf('Mon3 Delay: %d samples',bestFlatnessDelay)});
xlim([sampleRange(1) sampleRange(end)])
xlabel('Sample No.')
ylabel('Phase [degrees]');
legend('Mon2','Delayed Mon3','SimulatedFF','RealFF')
format_plots;
savePlot(saveDir,'minFlatnessPhaseAlong');

figure;
plot(mon2Mean,'LineWidth',2);
hold all;
plot(mon3Mean,'LineWidth',2);
plot(simFFResults{sampDelay==0}.meanSimFFAlongPulse,'LineWidth',2); 
plot(ffOn.meanPhaseAlongPulse(3,:),'LineWidth',2)
title('Original Sim, Mean Phase Along Pulse');
xlim([sampleRange(1) sampleRange(end)])
xlabel('Sample No.')
ylabel('Phase [degrees]');
legend('Mon2','Mon3','SimulatedFF','RealFF')
format_plots;
savePlot(saveDir,'noDelayFlatnessPhaseAlong');

figure;
plot(mon2Mean,'LineWidth',2);
hold all;
plot(mon3Mean,'LineWidth',2);
plot(delaySignal(mon3Mean,bestFlatnessDelay),'LineWidth',2);
title('Mean Phase Along Pulse');
xlim(xLim);
xlabel('Sample No.')
ylabel('Phase [degrees]');
legend('Mon2','Mon3',sprintf('Mon3 Delayed %d samples',bestFlatnessDelay));
format_plots;
savePlot(saveDir,'minFlatnessPhaseAlong_wholePulse');

figure;
plot(diode2Mean,'LineWidth',2);
hold all;
plot(diode3Mean,'LineWidth',2);
plot(delaySignal(diode3Mean,bestFlatnessDelay),'LineWidth',2);
title({'Mean Diode, Max Sample Correlation' sprintf('Mon3 Delay: %d samples',bestDelay)});
xlim(xLim);
xlabel('Sample No.')
ylabel('Output [V]');
legend('Mon2','Mon3',sprintf('Mon3 Delayed %d samples',bestFlatnessDelay));
format_plots;
savePlot(saveDir,'minFlatnessDiode');

figure;
plot(mon2Std,'LineWidth',2);
hold all;
plot(mon3Std,'LineWidth',2);
plot(simFFResults{sampDelay==0}.stdSimFFAlongPulse,'LineWidth',2); 
plot(ffOn.stdPhaseAlongPulse(3,:),'LineWidth',2)
title('Original Sim, Jitter Along Pulse');
xlim([sampleRange(1) sampleRange(end)])
xlabel('Sample No.')
ylabel('Jitter [degrees]');
legend('Mon2','Mon3','SimulatedFF','RealFF')
format_plots;
savePlot(saveDir,'originalJitterAlong');

figure;
plot(mon2Std,'LineWidth',2);
hold all;
plot(delaySignal(mon3Std,bestFlatnessDelay),'LineWidth',2);
plot(simFFResults{bestFlatnessInd}.stdSimFFAlongPulse,'LineWidth',2); 
plot(ffOn.stdPhaseAlongPulse(3,:),'LineWidth',2)
title({'Best Sim Flatness, Jitter Along Pulse'  sprintf('Mon3 Delay: %d samples',bestFlatnessDelay)});
xlim([sampleRange(1) sampleRange(end)])
xlabel('Sample No.')
ylabel('Jitter [degrees]');
legend('Mon2','Delayed Mon3','SimulatedFF','RealFF')
format_plots;
savePlot(saveDir,'minFlatnessJitterAlong');

%%
% sampleGain = (corrSamplesMix2Mix3.*(stdPhaseAlongPulse(3,:)./stdPhaseAlongPulse(2,:)));
% meanGain = corrMeanMix2Mix3*(stdMeanPulsePhase(3)/stdMeanPulsePhase(2));
% 
% figure;
% plot(meanPhaseAlongPulse(3,:)-sampleGain.*meanPhaseAlongPulse(2,:))
% hold all
% plot(meanPhaseAlongPulse(3,:)-meanGain.*meanPhaseAlongPulse(2,:))
% 
