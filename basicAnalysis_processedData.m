close all;
%%
%processedDataFile = '/home/jack/PhaseFeedforward/CTFData/201412/20141216_1810_Gate325to450_Gain30_WiggleAmplitudeAlong3/processed/20141216_1810_Gate325to450_Gain30_WiggleAmplitudeAlong3.mat';
%processedDataFile = '/home/jack/PhaseFeedforward/CTFData/201411/20141128_1252_FF_Gate340to460_GainK1_0_GainK2_0/processed/20141128_1252_FF_Gate340to460_GainK1_0_GainK2_0.mat';
%load(processedDataFile);
%dataPath = fileparts(processedDataFile);

datasetName = '20141128_1252_FF_Gate340to460_GainK1_0_GainK2_0';%'20141128_1140_ConstKick_Gate_390_420_DAC1_800_DAC2_-800';
processedData = loadProcessedData({datasetName});
processedData = processedData{1};
dataPath = processedData.dataDir;

saveDir = [dataPath '/plots/'];
if (~exist(saveDir,'dir'))
    mkdir(saveDir);
end

removeThreshold = 0.01; % remove pulses that are removeThreshold % away from max on diode
skipPercent = 20;
legendString = {'Mon1','Mon2','Mon3','PETS'};
plotColours = {'b','r','g','m'};
mon12StartSample = 150;
mon12EndSample = 240;
%mon3Delay = 660; %delay between mon1/2 and mon3 in ns
mon3Delay = 648;
%petsDelay = -748;
petsDelay = -742; % delay between PETS and Mon3 (will change depending on oasis setup) in ns.

plotData = true;

%%

diodes = processedData.phaseData.frascatiDiodes;
mixers = processedData.phaseData.frascatiMixers;
phases =  processedData.phaseData.frascatiPhases;
pets_Phase = processedData.phaseData.petsPhase;
pets_I = processedData.phaseData.petsI;
pets_Q = processedData.phaseData.petsQ;
frascatiTimePerSample = processedData.phaseData.frascatiTimePerSample;
petsTimePerSample = processedData.phaseData.petsTimePerSample;

[nMonitors,nPulses,nSamples] = size(mixers);
[~,petsNSamples] = size(pets_I);

% % to help find mon3 delay
% crosscorr(squeeze(meanPhases(1,:)),squeeze(meanPhases(3,:)),400)
% 
% % to help find pets Delay
% tmpFrascati = meanPhases(3,:);%-mean(meanPhases(3,startSamples(3):endSamples(3)));
% tmpPETS = pets_meanPhases; %- mean(pets_meanPhases(petsStartSample:petsEndSample));
modSamples = linspace(1,nSamples,nSamples.*(frascatiTimePerSample./petsTimePerSample));
% yy=spline(1:nSamples,tmpFrascati,linspace(1,nSamples,nSamples.*(frascatiTimePerSample./petsTimePerSample)));
% crosscorr(yy,tmpPETS,900)
% [~,maxInd] = max(crosscorr(yy,tmpPETS,900));
% delaySamples = -900:900;
% delaySamples = delaySamples(maxInd);
% 
% figure;
% plot(4*(1:nSamples),meanPhases(3,:));%-mean(meanPhases(3,startSamples(3):endSamples(3))))
% hold all;
% plot((1:petsNSamples)-delaySamples,pets_meanPhases);%-mean(pets_meanPhases(petsStartSample:petsEndSample)))
% plot((1:petsNSamples),pets_meanPhases);%-mean(pets_meanPhases(petsStartSample:petsEndSample)))

transmission = NaN*ones(nMonitors,nPulses);
isGoodPulseMon = false(nMonitors,nPulses);

for mon = 1:nMonitors
    [startSamp,endSamp] = getSignalRange(squeeze(nanmean(diodes(mon,:,:))),20);
    transmission(mon,:) = squeeze(nanmean(diodes(mon,:,startSamp:endSamp),3));
    isGoodPulseMon(mon,:) = (transmission(mon,:)-min(transmission(mon,:))) < removeThreshold;
end

isGoodPulseAll = ones(1,nPulses);
for mon = 1:nMonitors
    isGoodPulseAll = isGoodPulseAll & isGoodPulseMon(mon,:);
end

mixers = mixers(:,isGoodPulseAll,:);
diodes = diodes(:,isGoodPulseAll,:);
phases = phases(:,isGoodPulseAll,:);
pets_I = pets_I(isGoodPulseAll,:);
pets_Q = pets_Q(isGoodPulseAll,:);
pets_Phase = pets_Phase(isGoodPulseAll,:);
[~,nPulses,~] = size(mixers);


pulseWidth = mon12EndSample-mon12StartSample;
mon3StartSample = mon12StartSample + round(mon3Delay./frascatiTimePerSample);
mon3EndSample = mon3StartSample + pulseWidth;
startSamples = [mon12StartSample mon12StartSample mon3StartSample];
endSamples = [mon12EndSample mon12EndSample mon3EndSample];
frascatiTimeAxes = NaN*ones(nMonitors,nSamples);
for mon=1:nMonitors
    frascatiTimeAxes(mon,:) = ((1:nSamples)-startSamples(mon)).*frascatiTimePerSample;
end

[~,mod3StartInd] = min(abs(modSamples-startSamples(3)));
petsStartSample = mod3StartInd+petsDelay;%mon3StartSample - round(petsDelay./petsTimePerSample);
petsEndSample = petsStartSample + round(pulseWidth.*(frascatiTimePerSample./petsTimePerSample));
petsTimeAxis = ((1:petsNSamples)-petsStartSample).*petsTimePerSample;



% pulseStartSamples = NaN*ones(nMonitors,nPulses);
% pulseEndSamples = NaN*ones(nMonitors,nPulses);
% for mon=1:nMonitors
%     for p=1:nPulses
%         [pulseStartSamples(mon,p),pulseEndSamples(mon,p)] = getSignalRange(squeeze(diodes(mon,p,:)),skipPercent);
%     end
% end
% startSamples = round(nanmean(pulseStartSamples,2));
% endSamples = round(nanmean(pulseEndSamples,2));
% frascatiTimeAxes = NaN*ones(nMonitors,nSamples);
% for mon=1:nMonitors
%     frascatiTimeAxes(mon,:) = ((1:nSamples)-startSamples(mon)).*frascatiTimePerSample;
% end

% 
% petsPulseStartSamples = NaN*ones(nPulses,1);
% petsPulseEndSamples = NaN*ones(nPulses,1);
% for p=1:nPulses
%     [iStart,iEnd] = getSignalRange(squeeze(pets_I(p,:)),skipPercent);
%     %[qStart,qEnd] = getSignalRange(squeeze(pets_Q(p,:)),skipPercent);
%     %petsPulseStartSamples(p) = nanmean([iStart qStart]);
%     %petsPulseEndSamples(p) = nanmean([iEnd qEnd]);
%     petsPulseStartSamples(p) = iStart;
%     petsPulseEndSamples(p) = iEnd;
% end
% petsStartSample = round(nanmean(petsPulseStartSamples));
% petsEndSample = round(nanmean(petsPulseEndSamples));
% 
% petsTimeAxis = ((1:petsNSamples)-petsStartSample).*petsTimePerSample;

% 
% figure;
% plotHandles = NaN*ones(1,4);
% for mon=1:nMonitors
%     plotHandles(mon) = plot(pulseStartSamples(mon,:),plotColours{mon});
%     hold all;
%     plot(pulseEndSamples(mon,:),plotColours{mon});
% end
% plotHandles(4) = plot(petsPulseStartSamples,plotColours{4});
% plot(petsPulseEndSamples,plotColours{4});
% title('Calculated Pulse Start and End Points');
% xlabel('Pulse No.');
% ylabel('Start/End Sample');
% legend(plotHandles,legendString);

%% Frascati
meanDiodes = squeeze(nanmean(diodes,2));
meanMixers = squeeze(nanmean(mixers,2));
meanPhases = squeeze(nanmean(phases,2)); % mean per sample
stdPhases = squeeze(nanstd(phases,0,2)); % std per sample

for mon=1:nMonitors
    tmpSubtractMean = nanmean(meanPhases(mon,startSamples(mon):endSamples(mon)));
    meanPhases(mon,:) = meanPhases(mon,:) - tmpSubtractMean;
    
    for p=1:nPulses
        phases(mon,p,:) = phases(mon,p,:)-tmpSubtractMean;
    end
end



meanStdPhase = NaN*ones(1,nMonitors);
for mon=1:nMonitors
    meanStdPhase(mon) = squeeze(nanmean(stdPhases(mon,startSamples(mon):endSamples(mon)),2)); % mean sample std
end
% 
diffPhases12 = phases(1,:,startSamples(1):endSamples(1)) - phases(2,:,startSamples(2):endSamples(2));
diffPhases13 = phases(1,:,startSamples(1):endSamples(1)) - phases(3,:,startSamples(3):endSamples(3));
diffPhases23 = phases(2,:,startSamples(2):endSamples(2)) - phases(3,:,startSamples(3):endSamples(3));
% 
resolution12 = squeeze(nanstd(diffPhases12,0,3))./sqrt(2);
resolution13 = squeeze(nanstd(diffPhases13,0,3))./sqrt(2);
resolution23 = squeeze(nanstd(diffPhases23,0,3))./sqrt(2);
meanResolution12 = mean(resolution12);
meanResolution13 = mean(resolution13);
meanResolution23 = mean(resolution23);
% 
meanPulseDiodes = NaN*ones(nMonitors,nPulses);
meanPulseMixers = NaN*ones(nMonitors,nPulses);
meanPulsePhases = NaN*ones(nMonitors,nPulses); % mean per pulse
for mon = 1:nMonitors
    meanPulseDiodes(mon,:) = squeeze(nanmean(diodes(mon,:,startSamples(mon):endSamples(mon)),3));
    meanPulseMixers(mon,:) = squeeze(nanmean(mixers(mon,:,startSamples(mon):endSamples(mon)),3));
    meanPulsePhases(mon,:) = squeeze(nanmean(phases(mon,:,startSamples(mon):endSamples(mon)),3));
end
stdPulseDiodes =std(meanPulseDiodes,0,2); 
stdPulseMixers =std(meanPulseMixers,0,2);
stdPulsePhases =std(meanPulsePhases,0,2); % error bars on mean pulse phase
% 
corr12_meanPulsePhase = corrcoef(meanPulsePhases(1,:),meanPulsePhases(2,:));
corr13_meanPulsePhase = corrcoef(meanPulsePhases(1,:),meanPulsePhases(3,:));
corr23_meanPulsePhase = corrcoef(meanPulsePhases(2,:),meanPulsePhases(3,:));
corr12_meanPulsePhase = corr12_meanPulsePhase(1,2);
corr13_meanPulsePhase = corr13_meanPulsePhase(1,2);
corr23_meanPulsePhase = corr23_meanPulsePhase(1,2);

fit12_meanPulsePhase = polyfit(meanPulsePhases(1,:),meanPulsePhases(2,:),1);
fit13_meanPulsePhase = polyfit(meanPulsePhases(1,:),meanPulsePhases(3,:),1);
fit23_meanPulsePhase = polyfit(meanPulsePhases(2,:),meanPulsePhases(3,:),1);
fit12_meanPulsePhase = fit12_meanPulsePhase(1);
fit13_meanPulsePhase = fit13_meanPulsePhase(1);
fit23_meanPulsePhase = fit23_meanPulsePhase(1);
% 
meanDiffPhase12 = meanPhases(1,startSamples(1):endSamples(1)) - meanPhases(2,startSamples(2):endSamples(2));
meanDiffPhase13 = meanPhases(1,startSamples(1):endSamples(1)) - meanPhases(3,startSamples(3):endSamples(3));
meanDiffPhase23 = meanPhases(2,startSamples(2):endSamples(2)) - meanPhases(3,startSamples(3):endSamples(3));
% 

phases_diffToMean = NaN*ones(nMonitors,nPulses,nSamples);
for mon=1:nMonitors
    for p=1:nPulses
        tmpMeanPhase = nanmean(squeeze(phases(mon,p,startSamples(mon):endSamples(mon))));
        phases_diffToMean(mon,p,:) = phases(mon,p,:) - tmpMeanPhase;
    end
end

corr12_samples = NaN*ones(1,pulseWidth+1);
corr13_samples = NaN*ones(1,pulseWidth+1);
corr23_samples = NaN*ones(1,pulseWidth+1);
corr12_diffToMean = NaN*ones(1,pulseWidth+1);
corr13_diffToMean = NaN*ones(1,pulseWidth+1);
corr23_diffToMean = NaN*ones(1,pulseWidth+1);
fit12_samples = NaN*ones(1,pulseWidth+1);
fit13_samples = NaN*ones(1,pulseWidth+1);
fit23_samples = NaN*ones(1,pulseWidth+1);
fit12_diffToMean = NaN*ones(1,pulseWidth+1);
fit13_diffToMean = NaN*ones(1,pulseWidth+1);
fit23_diffToMean = NaN*ones(1,pulseWidth+1);
for s=1:pulseWidth+1
    sMon1 = startSamples(1)+s-1;
    sMon2 = startSamples(2)+s-1;
    sMon3 = startSamples(3)+s-1;

    tmp12 = corrcoef(squeeze(phases(1,:,sMon1)),squeeze(phases(2,:,sMon2)));
    tmp13 = corrcoef(squeeze(phases(1,:,sMon1)),squeeze(phases(3,:,sMon3)));
    tmp23 = corrcoef(squeeze(phases(2,:,sMon2)),squeeze(phases(3,:,sMon3)));
    corr12_samples(s) = tmp12(1,2);
    corr13_samples(s) = tmp13(1,2);
    corr23_samples(s) = tmp23(1,2);
    
    tmp12 = polyfit(squeeze(phases(1,:,sMon1)),squeeze(phases(2,:,sMon2)),1);
    tmp13 = polyfit(squeeze(phases(1,:,sMon1)),squeeze(phases(3,:,sMon3)),1);
    tmp23 = polyfit(squeeze(phases(2,:,sMon2)),squeeze(phases(3,:,sMon3)),1);
    fit12_samples(s) = tmp12(1);
    fit13_samples(s) = tmp13(1);
    fit23_samples(s) = tmp23(1);

    tmp12 = corrcoef(squeeze(phases_diffToMean(1,:,sMon1)),squeeze(phases_diffToMean(2,:,sMon2)));
    tmp13 = corrcoef(squeeze(phases_diffToMean(1,:,sMon1)),squeeze(phases_diffToMean(3,:,sMon3)));
    tmp23 = corrcoef(squeeze(phases_diffToMean(2,:,sMon2)),squeeze(phases_diffToMean(3,:,sMon3)));
    corr12_diffToMean(s) = tmp12(1,2);
    corr13_diffToMean(s) = tmp13(1,2);
    corr23_diffToMean(s) = tmp23(1,2);
    
    tmp12 = polyfit(squeeze(phases_diffToMean(1,:,sMon1)),squeeze(phases_diffToMean(2,:,sMon2)),1);
    tmp13 = polyfit(squeeze(phases_diffToMean(1,:,sMon1)),squeeze(phases_diffToMean(3,:,sMon3)),1);
    tmp23 = polyfit(squeeze(phases_diffToMean(2,:,sMon2)),squeeze(phases_diffToMean(3,:,sMon3)),1);
    fit12_diffToMean(s) = tmp12(1);
    fit13_diffToMean(s) = tmp13(1);
    fit23_diffToMean(s) = tmp23(1);

end

%% PETS
pets_meanPhases = nanmean(pets_Phase); % mean per sample
tmpSubtractMean = nanmean(pets_meanPhases(petsStartSample:petsEndSample));
pets_meanPhases = pets_meanPhases - tmpSubtractMean;
for p=1:nPulses
    pets_Phase(p,:) = pets_Phase(p,:) - tmpSubtractMean;
end
pets_meanPulsePhases = nanmean(pets_Phase(:,petsStartSample:petsEndSample),2); % mean per pulse
pets_stdPhases = nanstd(pets_Phase); % std per sample
pets_meanStdPhase = nanmean(pets_stdPhases(petsStartSample:petsEndSample)); % mean sample std
pets_stdPulsePhase = std(pets_meanPulsePhases); 

corr_PetsFrascati3_meanPulsePhase = corrcoef(meanPulsePhases(3,:),pets_meanPulsePhases);
corr_PetsFrascati3_meanPulsePhase = corr_PetsFrascati3_meanPulsePhase(1,2);
fit_PetsFrascati3_meanPulsePhase = polyfit(meanPulsePhases(3,:)',pets_meanPulsePhases,1);
fit_PetsFrascati3_meanPulsePhase = fit_PetsFrascati3_meanPulsePhase(1);

corr_PetsFrascati1_meanPulsePhase = corrcoef(meanPulsePhases(1,:),pets_meanPulsePhases);
corr_PetsFrascati1_meanPulsePhase = corr_PetsFrascati1_meanPulsePhase(1,2);
fit_PetsFrascati1_meanPulsePhase = polyfit(meanPulsePhases(1,:)',pets_meanPulsePhases,1);
fit_PetsFrascati1_meanPulsePhase = fit_PetsFrascati1_meanPulsePhase(1);

corr_PetsFrascati2_meanPulsePhase = corrcoef(meanPulsePhases(2,:),pets_meanPulsePhases);
corr_PetsFrascati2_meanPulsePhase = corr_PetsFrascati2_meanPulsePhase(1,2);
fit_PetsFrascati2_meanPulsePhase = polyfit(meanPulsePhases(2,:)',pets_meanPulsePhases,1);
fit_PetsFrascati2_meanPulsePhase = fit_PetsFrascati2_meanPulsePhase(1);


%% Plots
if (~plotData)
    return;
end

figure;
plotHandles = NaN*ones(1,4);
for mon=1:nMonitors
    plotHandles(mon) = plot(frascatiTimeAxes(mon,startSamples(mon):endSamples(mon)),squeeze(meanPhases(mon,startSamples(mon):endSamples(mon))),plotColours{mon});
    hold all;
end
plotHandles(4) = plot(petsTimeAxis(petsStartSample:petsEndSample),pets_meanPhases(petsStartSample:petsEndSample),plotColours{4});
legend(plotHandles,legendString);
xlabel('Time [ns]');
ylabel('Phase [degrees]');
title('Mean Pulse Phase');
set(gca,'FontSize',18)
set(findall(gcf,'type','text'),'FontSize',18)
print([saveDir 'MeanPhaseVsSampleNo.png'],'-dpng');

figure;
plotHandles = NaN*ones(1,4);
legStrs = cell(1,4);
for mon=1:nMonitors
    plotHandles(mon) = plot(frascatiTimeAxes(mon,startSamples(mon):endSamples(mon)),squeeze(stdPhases(mon,startSamples(mon):endSamples(mon))),plotColours{mon});
    hold all;
    legStrs{mon} = sprintf('Mon%d (mean: %.2f^o)',mon,mean(stdPhases(mon,startSamples(mon):endSamples(mon))));
end
plotHandles(4) = plot(petsTimeAxis(petsStartSample:petsEndSample),pets_stdPhases(petsStartSample:petsEndSample),plotColours{4});
legStrs{4} = sprintf('PETS (mean: %.2f^o)',mean(pets_stdPhases(petsStartSample:petsEndSample)));
legend(plotHandles,legStrs);
xlabel('Time [ns]');
ylabel('Std Phase [degrees]');
title('Phase Jitter');
set(gca,'FontSize',18)
set(findall(gcf,'type','text'),'FontSize',18)
print([saveDir 'StdPhaseVsSampleNo.png'],'-dpng');

figure;
tmpTime = frascatiTimePerSample.*((1:length(resolution12))-frascatiTimePerSample);
plot(tmpTime,resolution12,'b');
hold all;
plot(tmpTime,resolution13,'r');
plot(tmpTime,resolution23,'g');
title('Resolution vs. Sample No.');
leg1Str = sprintf('Mon1-Mon2 (%.2f^o)',meanResolution12);
leg2Str = sprintf('Mon1-Mon3 (%.2f^o)',meanResolution13);
leg3Str = sprintf('Mon2-Mon3 (%.2f^o)',meanResolution23);
legend(leg1Str,leg2Str,leg3Str);
xlabel('Sample No.');
ylabel('Resolution [12GHz Degrees]');
print([saveDir 'Resolution.png'],'-dpng');

figure;
plot(meanPulsePhases(1,:),meanPulsePhases(2,:),'o');
xlabel('Mon1 Phase [12GHz Degrees]');
ylabel('Mon2 Phase [12GHz Degrees]');
title(sprintf('CT Mon1 - CT Mon2: Corr: %.2f, Grad: %.2f',corr12_meanPulsePhase,fit12_meanPulsePhase));
set(gca,'FontSize',18)
set(findall(gcf,'type','text'),'FontSize',18)
grid on;
print([saveDir 'CorrelationMeanMon1Mon2.png'],'-dpng');

figure;
plot(meanPulsePhases(1,:),meanPulsePhases(3,:),'o');
hold all
xlabel('Mon1 Phase [12GHz Degrees]');
ylabel('Mon3 Phase [12GHz Degrees]');
title(sprintf('CT Mon1 - CLEX Mon3: Corr: %.2f, Grad: %.2f',corr13_meanPulsePhase,fit13_meanPulsePhase));
set(gca,'FontSize',18)
set(findall(gcf,'type','text'),'FontSize',18)
grid on;
print([saveDir 'CorrelationMeanMon1Mon3.png'],'-dpng');

% 
figure;
plot(meanPulsePhases(2,:),meanPulsePhases(3,:),'o');
hold all
xlabel('Mon2 Phase [12GHz Degrees]');
ylabel('Mon3 Phase [12GHz Degrees]');
title(sprintf('Mon2 - Mon3: Corr: %.2f, Grad: %.2f',corr23_meanPulsePhase,fit23_meanPulsePhase));
print([saveDir 'CorrelationMeanMon2Mon3.png'],'-dpng');

% 
% 
figure;
plot(meanPulsePhases(1,:),pets_meanPulsePhases,'o');
xlabel('Mon3 Phase [degrees]');
ylabel('PETS Phase [degrees]');
title(sprintf('CT Mon1 - CLEX PETS: Corr: %.2f, Grad: %.2f',corr_PetsFrascati1_meanPulsePhase,fit_PetsFrascati1_meanPulsePhase));
set(gca,'FontSize',18)
set(findall(gcf,'type','text'),'FontSize',18)
grid on;
print([saveDir 'CorrelationMeanMon1PETS.png'],'-dpng');

figure;
plot(meanPulsePhases(2,:),pets_meanPulsePhases,'o');
xlabel('Mon3 Phase [12GHz Degrees]');
ylabel('PETS Phase [12GHz Degrees]');
title(sprintf('Mon2 - PETS: Corr: %.2f, Grad: %.2f',corr_PetsFrascati2_meanPulsePhase,fit_PetsFrascati2_meanPulsePhase));
print([saveDir 'CorrelationMeanMon2PETS.png'],'-dpng');

figure;
plot(meanPulsePhases(3,:),pets_meanPulsePhases,'o');
xlabel('Mon3 Phase [12GHz Degrees]');
ylabel('PETS Phase [12GHz Degrees]');
set(gca,'FontSize',18)
set(findall(gcf,'type','text'),'FontSize',18)
title(sprintf('CLEX Mon3 - CLEX PETS: Corr: %.2f, Grad: %.2f',corr_PetsFrascati3_meanPulsePhase,fit_PetsFrascati3_meanPulsePhase));
print([saveDir 'CorrelationMeanMon3PETS.png'],'-dpng');

figure;
plot(corr12_samples,'b');
hold all;
plot(corr13_samples,'r');
plot(corr23_samples,'g');
legend('Mon1-Mon2','Mon1-Mon3','Mon2-Mon3','Location','best');
xlabel('Sample No.');
ylabel('Corr Coeff');
title('Correlation vs. Sample No.');
print([saveDir 'CorrelationVsSampleNo.png'],'-dpng');


figure;
plot(corr12_diffToMean,'b');
hold all;
plot(corr13_diffToMean,'r');
plot(corr23_diffToMean,'g');
legend('Mon1-Mon2','Mon1-Mon3','Mon2-Mon3','Location','best');
xlabel('Sample No.');
ylabel('Corr Coeff');
title('Correlation in Difference to Mean Phase vs. Sample No.');
print([saveDir 'CorrelationDiffToMeanVsSampleNo.png'],'-dpng');

% 
figure;
plot(meanPulsePhases(1,:),'b');
hold all;
plot(meanPulsePhases(2,:)+5,'r')
plot(meanPulsePhases(3,:)+10,'g')
plot(pets_meanPulsePhases+15,'m');
title('Mean Phase vs. Time');
leg1 = sprintf('Mon1 (std = %.2f^o)',stdPulsePhases(1)); 
leg2 = sprintf('Mon2 (std = %.2f^o)',stdPulsePhases(2)); 
leg3 = sprintf('Mon3 (std = %.2f^o)',stdPulsePhases(3)); 
leg4 = sprintf('PETS (std = %.2f^o)',pets_stdPulsePhase); 
legend(leg1,leg2,leg3,leg4);
xlabel('Time [Pulse No.]');
ylabel('Phase [12GHz Degrees]');
print([saveDir 'MeanPhaseVsTime.png'],'-dpng');


figure;
plot(meanPulsePhases(2,:)-meanPulsePhases(1,:),'b');
title(sprintf('Mon2-Mon1 Mean Phase vs. Time (std: %.2f^o)',std(meanPulsePhases(2,:)-meanPulsePhases(1,:))));
xlabel('Time [Pulse No.]');
ylabel('Phase [degrees]');
ylim([-5 5]);
print([saveDir 'DiffVsTime_Mon1Mon2.png'],'-dpng');


figure;
plot(meanPulsePhases(3,:)-meanPulsePhases(1,:),'b');
title(sprintf('Mon3-Mon1 Mean Phase vs. Time (std: %.2f^o)',std(meanPulsePhases(3,:)-meanPulsePhases(1,:))));
xlabel('Time [Pulse No.]');
ylabel('Phase [degrees]');
ylim([-5 5]);
print([saveDir 'DiffVsTime_Mon1Mon3.png'],'-dpng');

figure;
plot(meanPulsePhases(3,:)-meanPulsePhases(2,:),'b');
title(sprintf('Mon3-Mon2 Mean Phase vs. Time (std: %.2f^o)',std(meanPulsePhases(3,:)-meanPulsePhases(2,:))));
xlabel('Time [Pulse No.]');
ylabel('Phase [degrees]');
ylim([-5 5]);
print([saveDir 'DiffVsTime_Mon2Mon3.png'],'-dpng');

figure;
plot(pets_meanPulsePhases'-meanPulsePhases(3,:),'b');
title(sprintf('PETS-Mon3 Mean Phase vs. Time (std: %.2f^o)',std(pets_meanPulsePhases'-meanPulsePhases(3,:))));
xlabel('Time [Pulse No.]');
ylabel('Phase [degrees]');
ylim([-5 5]);
print([saveDir 'DiffVsTime_PETSMon3.png'],'-dpng');

figure;
plot(pets_meanPulsePhases'-meanPulsePhases(1,:),'b');
title(sprintf('PETS-Mon1 Mean Phase vs. Time (std: %.2f^o)',std(pets_meanPulsePhases'-meanPulsePhases(1,:))));
xlabel('Time [Pulse No.]');
ylabel('Phase [degrees]');
ylim([-5 5]);
print([saveDir 'DiffVsTime_PETSMon1.png'],'-dpng');

figure;
plot(pets_meanPulsePhases'-meanPulsePhases(2,:),'b');
title(sprintf('PETS-Mon2 Mean Phase vs. Time (std: %.2f^o)',std(pets_meanPulsePhases'-meanPulsePhases(2,:))));
xlabel('Time [Pulse No.]');
ylabel('Phase [degrees]');
ylim([-5 5]);
print([saveDir 'DiffVsTime_PETSMon2.png'],'-dpng');

figure;
diffP


% % 
% figure;
% plot(meanDiffPhase12);
% hold all;
% plot(meanDiffPhase13);
% plot(meanDiffPhase23);
% title('Mean Diff Phase');
% legend('1-2','1-3','2-3');
% xlabel('Time [Pulse No.]');
% ylabel('Phase [12GHz Degrees]');
% % 

% figure;
% shadedErrorBar(frascatiTimeAxes(1,startSamples(1):endSamples(1)),meanPhases(1,startSamples(1):endSamples(1)),stdPhases(1,startSamples(1):endSamples(1)),'b',1)
% hold all
% shadedErrorBar(frascatiTimeAxes(2,startSamples(2):endSamples(2)),meanPhases(2,startSamples(2):endSamples(2)),stdPhases(2,startSamples(2):endSamples(2)),'r',1)
% shadedErrorBar(frascatiTimeAxes(3,startSamples(3):endSamples(3)),meanPhases(3,startSamples(3):endSamples(3)),stdPhases(3,startSamples(3):endSamples(3)),'g',1)
% 
% figure;
% shadedErrorBar(petsTimeAxis(petsStartSample:petsEndSample),pets_meanPhases(petsStartSample:petsEndSample),pets_stdPhases(petsStartSample:petsEndSample),'m',1)
