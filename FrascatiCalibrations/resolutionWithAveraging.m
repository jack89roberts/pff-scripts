close all; clearvars;
%%
processedDatFile = '/home/jack/PhaseFeedforward/Analysis/201511/20151116_1545_Resolution_161115.mat';
samplesToAverage = 1:2:21;
plotSaveDir = '/home/jack/PhaseFeedforward/Analysis/201511/Plots/20151116_1545_ResolutionFONT5aWithAveraging';

%%
addpath('../');
load(processedDatFile);

nPoints = length(samplesToAverage);
averagedPhases = NaN(nPoints,nMons,nPulses,nSamples);

phases(2,:,:) = delaySignal(squeeze(phases(2,:,:)),3);

for p=1:nPoints
    fprintf('Averaging %d samples...\n',samplesToAverage(p));
    newPhases = NaN(nMons,nPulses,nSamples);
    
    for m=1:nMons
        newPhases(m,:,:) = phases(m,:,:)-subtractPhase(m);
        newPhases(m,:,:) = averageSamples(squeeze(newPhases(m,:,:)),samplesToAverage(p));
    end
    
    averagedPhases(p,:,:,:) = newPhases;   
end


meanAveragedPhases = nanmean(averagedPhases,3);
diff12 = squeeze(averagedPhases(:,1,:,:)-averagedPhases(:,2,:,:));
res12AvgPhases = nanstd(diff12,1,2)./sqrt(2);

%%
figure;
for p=1:nPoints
    plot(squeeze(meanAveragedPhases(p,1,:)));
    hold all
end
title('Mon1')
%xlabel();
%savePlot(plotSaveDir,'mon1');

figure;
for p=1:nPoints
    plot(squeeze(meanAveragedPhases(p,2,:)));
    hold all
end
title('Mon2')

figure;
for p=1:nPoints
plot(res12AvgPhases(p,:))
hold all;
end
title('Resolution Along Pulse')
xlim([440 820]);
ylim([0.3 0.4]);
xlabel('Sample No.');
ylabel('Resolution [degrees]')
format_plots;
savePlot(plotSaveDir,'resolutionAlongPulse');

figure;
plot(samplesToAverage,res12AvgPhases(:,650),'o-')
title('Resolution vs. No. Samples Averaged')
xlabel('No. Samples Averaged');
ylabel('Resolution [degrees]')
xlim([1 21])
format_plots;
savePlot(plotSaveDir,'resolutionVsNAvg');

figure
plot(squeeze(averagedPhases(1,1,1,:)))
hold all
plot(squeeze(averagedPhases(1,2,1,:)))