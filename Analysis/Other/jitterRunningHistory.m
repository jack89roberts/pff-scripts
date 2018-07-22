% load PFF off dataset to data{1} and PFF on to data{2} first

repRate=2.4; % pulses/s (2.4 for interleaved)
lastGoodPulse = 831;%863;%831;

%% Min jitter across different no. pulses in these datasets
avgArr = 25:1:830;

nPulses = min([data{1}.nPulses data{2}.nPulses]);
minJitterFFOn = NaN(1,length(avgArr));
minJitterFFOff = NaN(1,length(avgArr));
minJitterFFOnInd = NaN(1,length(avgArr));
minJitterFFOffInd = NaN(1,length(avgArr));
maxJitterRatio = NaN(1,length(avgArr));
maxJitterRatioInd = NaN(1,length(avgArr));

for i=1:length(avgArr)
    nAvg = avgArr(i);
    fprintf('%d\n',nAvg);

    stdNextNPulses_FFOn = NaN(1,nPulses);
    stdNextNPulses_FFOff = NaN(1,nPulses);

    stopPulse = min([data{1}.nPulses data{2}.nPulses])-nAvg+1;
    
    p=1;
    while ( p <= stopPulse )
        stdNextNPulses_FFOn(p) = nanstd(data{2}.meanPulsePhase(3,p:(p+nAvg-1)));
        stdNextNPulses_FFOff(p) = nanstd(data{1}.meanPulsePhase(3,p:(p+nAvg-1)));
        p=p+1;
    end

    [minJitterFFOn(i),minJitterFFOnInd(i)] = min(stdNextNPulses_FFOn(1:(lastGoodPulse-nAvg)));
    [minJitterFFOff(i),minJitterFFOffInd(i)] = min(stdNextNPulses_FFOff(1:(lastGoodPulse-nAvg)));

    ratioNextNPulses = stdNextNPulses_FFOff./stdNextNPulses_FFOn;
    [maxJitterRatio(i),maxJitterRatioInd(i)] = max(ratioNextNPulses(1:(lastGoodPulse-nAvg)));
end

figure;
subplot(3,1,1:2)

plot(repRate.*avgArr,minJitterFFOn,'r','LineWidth',2);
hold all;
plot(repRate.*avgArr,minJitterFFOff,'b','LineWidth',2);
plot([avgArr(1) avgArr(end)].*repRate, [0.2 0.2],'k','LineWidth',2)
plot([avgArr(1) avgArr(end)].*repRate, [0.3 0.3],'k--','LineWidth',1)
legend('PFF On','PFF Off')
xlim([avgArr(1) avgArr(end)].*repRate)
xlabel('Time Period [s]')
ylabel('Phase Jitter [degrees]')
title('Lowest Phase Jitter vs. Time Period')

subplot(3,1,3)
plot(repRate.*avgArr,maxJitterRatio,'r','LineWidth',2);
title('Max Jitter Reduction Factor vs. Time Period')
xlabel('Time Period [s]')
ylabel('PFF Jitter: Off/On')
xlim([avgArr(1) avgArr(end)].*repRate)
ylim([1 6])
% figure;
% plot(repRate.*avgArr,maxJitterRatioInd,'r','LineWidth',2);

figure;
plot(repRate.*avgArr,minJitterFFOnInd,'r','LineWidth',2);
xlim([avgArr(1) avgArr(end)].*repRate)


%% Jitter across n pulses across this dataset
nAvg = 50;

stdNextNPulses_FFOn = NaN(1,nPulses);
stdNextNPulses_FFOff = NaN(1,nPulses);
corrNextNPulses_FFOn = NaN(1,nPulses);
corrNextNPulses_FFOff = NaN(1,nPulses);
stdRatFFOff = NaN(1,nPulses);
 
stopPulse = min([data{1}.nPulses data{2}.nPulses])-nAvg+1;

p=1;
while ( p <= stopPulse )
    stdNextNPulses_FFOn(p) = nanstd(data{2}.meanPulsePhase(3,p:(p+nAvg-1)));
    stdNextNPulses_FFOff(p) = nanstd(data{1}.meanPulsePhase(3,p:(p+nAvg-1)));
    corrNextNPulses_FFOn(p) = nancorrcoef(data{2}.meanPulsePhase(1,p:(p+nAvg-1)),data{2}.meanPulsePhase(3,p:(p+nAvg-1)));
    corrNextNPulses_FFOff(p) = nancorrcoef(data{1}.meanPulsePhase(1,p:(p+nAvg-1)),data{1}.meanPulsePhase(3,p:(p+nAvg-1)));
    stdRatFFOff(p) = stdNextNPulses_FFOff(p)./nanstd(data{1}.meanPulsePhase(1,p:(p+nAvg-1)));
    p=p+1;
end

timeAx = ((1:nPulses)-1)*repRate;

plotStart = timeAx(1);
plotEnd = timeAx(lastGoodPulse-nAvg);

figure;
subplot(3,1,1:2)
plot(timeAx,stdNextNPulses_FFOn,'r','LineWidth',2)
hold all
plot(timeAx,stdNextNPulses_FFOff,'b','LineWidth',2)
plot([timeAx(1) timeAx(end)], [0.2 0.2],'k','LineWidth',2)
plot([timeAx(1) timeAx(end)], [0.3 0.3],'k--','LineWidth',1)
xlim([plotStart plotEnd])
legend('PFF On','PFF Off')
title('Phase Jitter of Next 2 Minutes')
ylabel('Phase Jitter [degrees]')
xlabel('Time [s]');

subplot(3,1,3)
plot(timeAx,stdNextNPulses_FFOff./stdNextNPulses_FFOn,'k','LineWidth',2)
xlim([plotStart plotEnd])
title('Jitter Reduction Factor of Next 2 Minutes')
ylabel('PFF Jitter: Off/')
xlabel('Time [s]');

figure;
subplot(3,1,1:2)
plot(timeAx,corrNextNPulses_FFOn,'r','LineWidth',2)
hold all
plot(timeAx,corrNextNPulses_FFOff,'b','LineWidth',2)
xlim([plotStart plotEnd])
legend('PFF On','PFF Off')
title('U/S - D/S Correlation of Next 2 Minutes')
ylabel('Correlation Coefficient')
xlabel('Time [s]');

subplot(3,1,3)
plot(timeAx,stdRatFFOff,'r','LineWidth',2)
xlim([plotStart plotEnd])
title('U/S - D/S Jitter Ratio of Next 2 Minutes (PFF Off)')
ylabel('DS / US Jitter')
xlabel('Time [s]');


%% plot a sub-set of the data
pulseRange = 385:625;%294:343; %385:625; %1:899;%385:625; 

subPFFOff = data{1}.meanPulsePhase(3,pulseRange);
subPFFOn = data{2}.meanPulsePhase(3,pulseRange);
% subPFFOff = subPFFOff - nanmean(subPFFOff);
% subPFFOn = subPFFOn - nanmean(subPFFOn);

% PFF Off
figure;
upOff = data{1}.meanPulsePhase(1,pulseRange);
downOff = data{1}.meanPulsePhase(3,pulseRange);
plot(timeAx(pulseRange),downOff-nanmean(downOff),'b','LineWidth',2)
hold all
plot(timeAx(pulseRange),upOff-nanmean(upOff),'g','LineWidth',2)
legend('Downstream','Upstream')
xlim([plotStart plotEnd])
xlabel('Time [s]')
ylabel('Phase [degrees]')
title('Mean Upstream and Downstream Phase (PFF Off)')

figure;
stdUpOff = nanstd(upOff);
stdDownOff = nanstd(downOff);
histogram(downOff-nanmean(downOff),'FaceColor','b')
hold all
histogram(upOff-nanmean(upOff),'FaceColor','g')
legend(sprintf('Downstream (std=%.2f^o)',stdDownOff),sprintf('Upstream (std=%.2f^o)',stdUpOff))
title('Up and Downstream Phase Distribution (PFF Off)')
xlabel('Phase [degrees]')
ylabel('No. Pulses')

figure;
scatterColourProgression(upOff,downOff);%;,'o',8,1)
xlabel('Upstream Phase [degrees]')
ylabel('Downstream Phase [degrees]')
xlim([min(upOff) max(upOff)])
ylim([min(downOff) max(downOff)])
title(sprintf('Upstrem-Downstream Phase Correlation, PFF Off (%.2f)',nancorrcoef(upOff,downOff)))

% PFF ON
figure;
plot(timeAx(pulseRange),subPFFOff,'b','LineWidth',2);
hold all;
plot(timeAx(pulseRange),subPFFOn,'r','LineWidth',2);
legend('PFF Off','PFF On')
xlabel('Time [s]')
ylabel('Phase [degrees]')
title('Mean Downstream Phase')

figure;
stdDownOff = nanstd(subPFFOff);
stdDownOn = nanstd(subPFFOn);
histogram(subPFFOff,'FaceColor','b')
hold all
histogram(subPFFOn,'FaceColor','r')
legend(sprintf('PFF Off (std=%.2f^o)',stdDownOff),sprintf('PFF On (std=%.2f^o)',stdDownOn))
title('Downstream Phase Distribution')
xlabel('Phase [degrees]')
ylabel('No. Pulses')

figure;
upOn = data{2}.meanPulsePhase(1,pulseRange);
scatterColourProgression(upOn,subPFFOn);%;,'o',8,1)
xlabel('Upstream Phase [degrees]')
ylabel('Downstream Phase [degrees]')
xlim([min(upOn) max(upOn)])
ylim([min(subPFFOn) max(subPFFOn)])
title(sprintf('Upstrem-Downstream Phase Correlation, PFF On (%.2f)',nancorrcoef(upOn,subPFFOn)))


figure;
plot(upOff,subPFFOff,'o','MarkerFaceColor','b');
hold all;
plot(upOn,subPFFOn,'o','MarkerFaceColor','r')
xlabel('Upstrema Phase [degrees]')
ylabel('Downstream Phase [degrees]')
title('Upstrem-Downstream Phase Correlation')
legend(sprintf('PFF Off (corr=%.2f)',nancorrcoef(upOff,subPFFOff)), sprintf('PFF On (corr=%.2f)',nancorrcoef(upOn,subPFFOn)));
