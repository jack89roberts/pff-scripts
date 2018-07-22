load('/home/jack/PhaseFeedforward/FONTData/201610/Extracted/20161027_1545_NormalDAQCheck_271016.mat');
saveDir = '/home/jack/PhaseFeedforward/Analysis/201610/Plots/20161027_1545_NormalDAQCheck_271016';
savePlots = 1;

%%
nPulses = length(FONTData.Time);

seconds = NaN(1,nPulses);


seconds(1) = FONTData.Time(1,3);
initSecs = seconds(1);
initMins = FONTData.Time(1,2);
initHr = FONTData.Time(1,1);

for i=1:nPulses
    newSecs = FONTData.Time(i,3);
    newMins = FONTData.Time(i,2);
    newHr = FONTData.Time(i,1);
    
    seconds(i) = (newHr-initHr)*60*60 + (newMins-initMins)*60 + (newSecs-initSecs);

end

adc1 = squeeze(FONTData.ADCs(1,:,:));
[~,dupADC1] = removeDuplicatePulses(adc1);
adcDuplicates=find(dupADC1)

figure;
ax1=subplot(3,1,1:2);
plot(diff(seconds));
ylabel('Time Diff [s]')
xlabel('Pulse No.')
title('Time Stamp Diff: DAQ Saved Data')
hold all;
plot([0 nPulses],[2.400 2.400],'k','LineWidth',2);
plot([0 nPulses],[0 0],'k','LineWidth',2);
plot([0 nPulses],[1.800 1.800],'k--','LineWidth',1.5);
plot([0 nPulses],[0.600 0.600],'k--','LineWidth',1.5);
ax2=subplot(3,1,3);
plot(dupADC1);
title('Duplicate Pulses')
xlabel('Pulse No.')
linkaxes([ax1,ax2],'x')
savePlot(saveDir,'fontTimeStamps');