nDelay = -20:20;
data = load('/home/jack/PhaseFeedforward/Analysis/201511/20151120_1538_Gain-800_R56_0.1_Interleaaved_Odd.mat');
colors = varycolor(length(nDelay));
sampleRange = data.sampleRange;
flatness = NaN(1,length(nDelay));
figure;
for i=1:length(nDelay)
    down = squeeze(data.phases(3,:,:));
    down = delaySignal(down,nDelay(i));
    delayedData = data;
    delayedData.phases(3,:,:) = down;
    
    simFFResults = getSimulatedFF(delayedData);
    plot(simFFResults.meanSimFFAlongPulse,'Color',colors(i,:));
    hold all;
    
    flatness(i) = nanstd(simFFResults.meanSimFFAlongPulse(sampleRange));
end

figure;
plot(nDelay,flatness)

%%
figure;
subplot(1,2,1)
plot(data.meanPhaseAlongPulse(2,:));
hold all;
plot(data.meanPhaseAlongPulse(3,:));
simFFResults = getSimulatedFF(data,[],6);
plot(simFFResults.meanSimFFAlongPulse);
title('orig')

subplot(1,2,2)
plot(data.meanPhaseAlongPulse(2,:))
hold all;
plot(delaySignal(data.meanPhaseAlongPulse(3,:),6));
down = squeeze(data.phases(3,:,:));
down = delaySignal(down,6);
delayedData = data;
delayedData.phases(3,:,:) = down;
simFFResults = getSimulatedFF(delayedData,[],6);
plot(simFFResults.meanSimFFAlongPulse);
title('delayed')
