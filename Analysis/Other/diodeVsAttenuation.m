close all; clearvars;

addpath('../../');
addpath('../');
format_plots;
%% input parameters
dataSetNames = {...
    '20150716_1459_Mon2DiodeCheck_0dB_MinMixer',...
    '20150716_1501_Mon2DiodeCheck_0dB_MaxMixer',...
    '20150716_1509_Mon2DiodeCheck_1dB_MaxMixer',...
    '20150716_1509_Mon2DiodeCheck_1dB_MinMixer',...
    '20150716_1512_Mon2DiodeCheck_2dB_MinMixer',...
    '20150716_1514_Mon2DiodeCheck_2dB_MaxMixer',...
    '20150716_1515_Mon2DiodeCheck_3dB_MaxMixer',...
    '20150716_1516_Mon2DiodeCheck_3dB_MinMixer',...
    '20150716_1518_Mon2DiodeCheck_4dB_MinMixer',...
    '20150716_1519_Mon2DiodeCheck_4dB_MaxMixer',...
    '20150716_1521_Mon2DiodeCheck_5dB_MaxMixer',...
    '20150716_1522_Mon2DiodeCheck_5dB_MinMixer',...
    '20150716_1524_Mon2DiodeCheck_6dB_MinMixer',...
    '20150716_1525_Mon2DiodeCheck_6dB_MaxMixer',...
    '20150716_1528_Mon2DiodeCheck_7dB_MaxMixer',...
    '20150716_1529_Mon2DiodeCheck_7dB_MinMixer',...
    '20150716_1530_Mon2DiodeCheck_8dB_MinMixer',...
    '20150716_1531_Mon2DiodeCheck_8dB_MaxMixer',...
    '20150716_1536_Mon2DiodeCheck_9dB_MaxMixer',...
    '20150716_1537_Mon2DiodeCheck_9dB_MinMixer',...
    '20150716_1542_Mon2DiodeCheck_10dB_MinMixer',...
    '20150716_1543_Mon2DiodeCheck_10dB_MaxMixer',...
    '20150716_1545_Mon2DiodeCheck_11dB_MaxMixer',...
    '20150716_1546_Mon2DiodeCheck_11dB_MinMixer',...
    '20150716_1547_Mon2DiodeCheck_12dB_MinMixer',...
    '20150716_1548_Mon2DiodeCheck_12dB_MaxMixer',...
    '20150716_1550_Mon2DiodeCheck_13dB_MaxMixer',...
    '20150716_1550_Mon2DiodeCheck_13dB_MinMixer',...
    '20150716_1552_Mon2DiodeCheck_14dB_MinMixer',...
    '20150716_1553_Mon2DiodeCheck_14dB_MaxMixer',...
    '20150716_1555_Mon2DiodeCheck_15dB_MaxMixer',...
    '20150716_1555_Mon2DiodeCheck_15dB_MinMixer',...
    '20150716_1558_Mon2DiodeCheck_16dB_MinMixer',...
    '20150716_1559_Mon2DiodeCheck_16dB_MaxMixer',...
    '20150716_1601_Mon2DiodeCheck_17dB_MaxMixer',...
    '20150716_1602_Mon2DiodeCheck_17dB_MinMixer',...
    ...%'20150716_1604_Mon2DiodeCheck_18dB_MinMixer',...
    ...%'20150716_1605_Mon2DiodeCheck_18dB_MaxMixer',...
    ...%'20150716_1607_Mon2DiodeCheck_19dB_MaxMixer',...
    '20150716_1608_Mon2DiodeCheck_19dB_MaxMixer',...
    '20150716_1609_Mon2DiodeCheck_19dB_MinMixer',...
    '20150716_1612_Mon2DiodeCheck_20dB_MinMixer',...
    '20150716_1613_Mon2DiodeCheck_20dB_MaxMixer'...
    ...%'20150716_1616_Mon2DiodeCheck_25dB_MaxMixer',...
    ...%'20150716_1617_Mon2DiodeCheck_25dB_MinMixer'...
};

sampleRange = 300:330;
nSamples = [];
nPulses = [];
monIndex = 2;
nominalInputPowerdBm = 25;

mixerPointsToFit = 5:19;
diodePointsToFit = 15:19;
diodeZoomPlotPoints = [19 12];

saveDir = '/home/jack/PhaseFeedforward/Analysis/201507/DiodeVsAttenuation';
savePlots = 0;

%% load data + initial processing
nDataSets = length(dataSetNames);
attenuation = NaN(1,nDataSets);
mixerMinOrMax = NaN(1,nDataSets);

for d=1:length(dataSetNames)

    dataSetName = dataSetNames{d};
    splitName = strsplit(dataSetName,'_');
    attenuation(d) = str2double(splitName{4}(1:end-2));
    mixStatus = splitName{5}(1:3);
    if (strcmp(mixStatus,'Min'))
        mixerMinOrMax(d) = -1;
    elseif (strcmp(mixStatus,'Max'))
        mixerMinOrMax(d) = 1;
    end
          
    [CTFData, ~, dataDir] = loadMergedData(dataSetName);
    [tmpMixers,tmpDiodes] = extractMixerDiode(CTFData);
    tmpMixers = squeeze(tmpMixers(monIndex,:,:));
    tmpDiodes = squeeze(tmpDiodes(monIndex,:,:));
    
    if (isempty(nSamples))
        [nPulses,nSamples] = size(tmpDiodes);
        diodes = NaN(nDataSets,nPulses,nSamples); 
        mixers = NaN(nDataSets,nPulses,nSamples);
    end

    mixers(d,:,:) = tmpMixers;
    diodes(d,:,:) = tmpDiodes;
    
end
clear CTFData;

inputPowerdBm = nominalInputPowerdBm-attenuation;
inputPowerWatts = dBmToWatts(inputPowerdBm);
inputVoltage = dBmToVolts(inputPowerdBm);

meanDiodes = squeeze(nanmean(diodes,2));
meanMixers = squeeze(nanmean(mixers,2));

maxMixInd = mixerMinOrMax == 1;
minMixInd = mixerMinOrMax == -1;
diodesMaxMix = meanDiodes(maxMixInd,:);
diodesMinMix = meanDiodes(minMixInd,:);
mixersMax = meanMixers(maxMixInd,:);
mixersMin = meanMixers(minMixInd,:);
attenuationMaxMix = attenuation(maxMixInd);
attenuationMinMix = attenuation(minMixInd);
inputPowerMaxMixWatts = inputPowerWatts(maxMixInd);
inputPowerMinMixWatts = inputPowerWatts(minMixInd);
inputPowerMaxMixdBm = inputPowerdBm(maxMixInd);
inputPowerMinMixdBm = inputPowerdBm(minMixInd);
inputVoltageMaxMix = inputVoltage(maxMixInd);
inputVoltageMinMix = inputVoltage(minMixInd);

meanDiodesMaxMix = abs(nanmean(diodesMaxMix(:,sampleRange),2));
meanDiodesMinMix = abs(nanmean(diodesMinMix(:,sampleRange),2));
meanAmpDiodes = mean([meanDiodesMaxMix, meanDiodesMinMix],2);

meanSqrtDiodeMaxMix = sqrt(abs(nanmean(diodesMaxMix(:,sampleRange),2)));
meanSqrtDiodeMinMix = sqrt(abs(nanmean(diodesMinMix(:,sampleRange),2)));
meanAmpSqrtDiode = mean([meanSqrtDiodeMaxMix, meanSqrtDiodeMinMix],2);

meanMixerMax = abs(nanmean(mixersMax(:,sampleRange),2));
meanMixerMin = abs(nanmean(mixersMin(:,sampleRange),2));
meanAmpMixer = mean([meanMixerMax, meanMixerMin],2);

%% plots
figure;
for a=1:length(attenuationMaxMix)
    atten = attenuationMaxMix(a);
    
    plot(diodesMaxMix(a,:),'LineWidth',lineWidthBig);
    hold all;
    plot(diodesMinMix(a,:),'LineWidth',lineWidthBig);
    xlabel('Sample No.');
    ylabel('Output [V]');
    title(sprintf('DIODE: %.0f dB',atten));
    legend('At Max Mixer','At Min Mixer')
    xlim([200 400])
    ylim([-0.17 0.01])
    format_plots;
    if savePlots; savePlot(saveDir,sprintf('diode_%.0fdB',atten)); end;
    hold off;
    
    plot(-mixersMax(a,:),'LineWidth',lineWidthBig);
    hold all;
    plot(mixersMin(a,:),'LineWidth',lineWidthBig);
    xlabel('Sample No.');
    ylabel('Output [V]');
    title(sprintf('MIXER: %.0f dB',atten));
    legend('Max Mixer','Min Mixer')
    xlim([200 400])
    ylim([-1.25 0.1])
    format_plots;
    if savePlots; savePlot(saveDir,sprintf('mixer_%.0fdB',atten)); end;
    hold off;
end

figure;
plot(inputVoltageMaxMix,meanDiodesMaxMix,'o-','LineWidth',lineWidthBig);
hold all
plot(inputVoltageMinMix,meanDiodesMinMix,'o-','LineWidth',lineWidthBig);
drawnow();
legend('At Max Mixer','At Min Mixer','Location','SouthEast')
xlabel('Input [V]');
ylabel('Output [V]');
title('DIODE');
format_plots;
if savePlots; savePlot(saveDir,'diodeVsVoltage'); end;

figure;
plot(inputVoltageMaxMix,meanSqrtDiodeMaxMix,'o-','LineWidth',lineWidthBig);
hold all
plot(inputVoltageMinMix,meanSqrtDiodeMinMix,'o-','LineWidth',lineWidthBig);
drawnow();
legend('At Max Mixer','At Min Mixer','Location','SouthEast')
xlabel('Input [V]');
ylabel('Output [V]');
title('SQRT DIODE');
format_plots;
if savePlots; savePlot(saveDir,'sqrtDiodeVsVoltage'); end;

figure;
plot(inputVoltageMaxMix,meanMixerMax,'o-','LineWidth',lineWidthBig);
hold all;
plot(inputVoltageMinMix,meanMixerMin,'o-','LineWidth',lineWidthBig);
drawnow();
legend('Max Mixer','Min Mixer','Location','NorthWest')
xlabel('Input [V]');
ylabel('Output [V]');
title('MIXER');
format_plots;
if savePlots; savePlot(saveDir,'mixerVsVoltage'); end;

figure;
plot(meanDiodesMaxMix,meanMixerMax);
hold all;
plot(meanDiodesMinMix,meanMixerMin);
plot(meanDiodes,meanMixers)

% figure;
% plot(attenuationMaxMix,meanDiodesMaxMix,'o-','LineWidth',lineWidthBig);
% hold all
% plot(attenuationMinMix,meanDiodesMinMix,'o-','LineWidth',lineWidthBig);
% legend('At Max Mixer','At Min Mixer')
% xlabel('Attenuation [dB]');
% ylabel('Output [V]');
% title('DIODE');
% format_plots;
% if savePlots; savePlot(saveDir,'diodeVsAttenuation'); end;
% 
% figure;
% plot(attenuationMaxMix,meanSqrtDiodeMaxMix,'o-','LineWidth',lineWidthBig);
% hold all
% plot(attenuationMinMix,meanSqrtDiodeMinMix,'o-','LineWidth',lineWidthBig);
% legend('At Max Mixer','At Min Mixer')
% xlabel('Attenuation [dB]');
% ylabel('Output [V]');
% title('SQRT DIODE');
% format_plots;
% if savePlots; savePlot(saveDir,'sqrtDiodeVsAttenuation'); end;
% 
% figure;
% plot(attenuationMaxMix,meanMixersMax,'o-','LineWidth',lineWidthBig);
% hold all;
% plot(attenuationMaxMix,meanMixersMin,'o-','LineWidth',lineWidthBig);
% legend('Max Mixer','Min Mixer')
% xlabel('Attenuation [dB]');
% ylabel('Output [V]');
% title('MIXER (log axis)');
% format_plots;
% if savePlots; savePlot(saveDir,'mixerVsAttenuation'); end;
% 
% figure;
% semilogy(attenuationMaxMix,meanMixersMax,'o-','LineWidth',lineWidthBig);
% hold all;
% semilogy(attenuationMaxMix,meanMixersMin,'o-','LineWidth',lineWidthBig);
% legend('Max Mixer','Min Mixer')
% xlabel('Attenuation [dB]');
% ylabel('Output [V]');
% title('MIXER');
% format_plots;
% if savePlots; savePlot(saveDir,'mixerVsAttenuationLog'); end;
% 
% 
% figure;
% plot(inputPowerMaxMixWatts,meanDiodesMaxMix,'o-','LineWidth',lineWidthBig);
% hold all
% plot(inputPowerMinMixWatts,meanDiodesMinMix,'o-','LineWidth',lineWidthBig);
% legend('At Max Mixer','At Min Mixer')
% xlabel('Estimated Input Power [W]');
% ylabel('Output [V]');
% title('DIODE');
% format_plots;
% if savePlots; savePlot(saveDir,'diodeVsPower'); end;
% 
% figure;
% plot(inputPowerMaxMixWatts,meanSqrtDiodeMaxMix,'o-','LineWidth',lineWidthBig);
% hold all
% plot(inputPowerMinMixWatts,meanSqrtDiodeMinMix,'o-','LineWidth',lineWidthBig);
% legend('At Max Mixer','At Min Mixer')
% xlabel('Estimated Input Power [W]');
% ylabel('Output [V]');
% title('SQRT DIODE');
% format_plots;
% if savePlots; savePlot(saveDir,'sqrtDiodeVsPower'); end;
% 
% 
% figure;
% plot(inputPowerMaxMixWatts,meanMixersMax,'o-','LineWidth',lineWidthBig);
% hold all;
% plot(inputPowerMinMixWatts,meanMixersMin,'o-','LineWidth',lineWidthBig);
% legend('Max Mixer','Min Mixer')
% xlabel('Estimated Input Power [W]');
% ylabel('Output [V]');
% title('MIXER');
% format_plots;
% if savePlots; savePlot(saveDir,'mixerVsPower'); end;

%% fits

linFitMixer = polyfit(inputVoltageMaxMix(mixerPointsToFit),meanAmpMixer(mixerPointsToFit)',1);
linFitDiode = polyfit(inputVoltageMaxMix(diodePointsToFit),meanAmpDiodes(diodePointsToFit)',1);
linFitSqrtDiode = polyfit(inputVoltageMaxMix(diodePointsToFit),meanAmpSqrtDiode(diodePointsToFit)',1);
linFitSqrtDiodedBm = polyfit(inputPowerMaxMixdBm(diodePointsToFit),meanAmpSqrtDiode(diodePointsToFit)',1);

plotFitMixer = linFitMixer(1).*inputVoltageMaxMix + linFitMixer(2);
plotFitDiode = linFitDiode(1).*inputVoltageMaxMix + linFitDiode(2);
plotFitSqrtDiode = linFitSqrtDiode(1).*inputVoltageMaxMix + linFitSqrtDiode(2);
plotFitSqrtDiodedBm = linFitSqrtDiodedBm(1).*inputPowerMaxMixdBm + linFitSqrtDiodedBm(2);

figure;
plot(inputVoltageMaxMix,meanAmpDiodes,'o-','LineWidth',1);
hold all
plot(inputVoltageMaxMix,plotFitDiode,'--','LineWidth',2);
drawnow();
legend('Data','Fit','Location','SouthEast')
xlabel('Input [V]');
ylabel('Output [V]');
title('DIODE');
format_plots;
if savePlots; savePlot(saveDir,'fitDiodeVsVoltage'); end;

xlim(inputVoltageMaxMix(diodeZoomPlotPoints));
if savePlots; savePlot(saveDir,'fitDiodeVsVoltage_zoom'); end;

figure;
plot(inputVoltageMaxMix,meanAmpSqrtDiode,'o-','LineWidth',1);
hold all
plot(inputVoltageMaxMix,plotFitSqrtDiode,'--','LineWidth',2);
drawnow();
legend('Data','Fit','Location','SouthEast')
xlabel('Input [V]');
ylabel('Output [V]');
title('SQRT DIODE');
format_plots;
if savePlots; savePlot(saveDir,'fitSqrtDiodeVsVoltage'); end;

xlim(inputVoltageMaxMix(diodeZoomPlotPoints));
if savePlots; savePlot(saveDir,'fitSqrtDiodeVsVoltage_zoom'); end;

figure;
plot(inputPowerMaxMixdBm,meanAmpSqrtDiode,'o-','LineWidth',1);
hold all
plot(inputPowerMaxMixdBm,plotFitSqrtDiodedBm,'--','LineWidth',2);
drawnow();
legend('Data','Fit','Location','SouthEast')
xlabel('Input Power [dBm]');
ylabel('Output [V]');
title('SQRT DIODE vs. Input Power');
format_plots;
if savePlots; savePlot(saveDir,'fitSqrtDiodeVsdBm'); end;

xlim(inputPowerMaxMixdBm(diodeZoomPlotPoints));
if savePlots; savePlot(saveDir,'fitSqrtDiodeVsdBm_zoom'); end;

figure;
plot(inputVoltageMaxMix,meanAmpMixer,'o-','LineWidth',1);
hold all;
plot(inputVoltageMaxMix,plotFitMixer,'--','LineWidth',2);
drawnow();
legend('Max Mixer','Min Mixer','Location','NorthWest')
xlabel('Input [V]');
ylabel('Output [V]');
title('MIXER');
format_plots;
if savePlots; savePlot(saveDir,'fitMixerVsVoltage'); end;

%% estimate phase resolution

%% save results to a structure;
if savePlots; save([saveDir '/processedData.mat']); end;
    