% amplifierMonitoring_scan
%
% Takes amplifier monitoring signals and makes some plots/fits, e.g.
% linearity, flatness along pulse etc.
%
% Then does the same for downstream phase signals and compares results.

close all; clearvars;
%% inputs

dataSetNames = {...
'20151123_1407_ConstKick_DAC_-4000_Interleaved',...
'20151123_1410_ConstKick_DAC_-3500_Interleaved',...
'20151123_1413_ConstKick_DAC_-3000_Interleaved',...
'20151123_1415_ConstKick_DAC_-2500_Interleaved',...
'20151123_1417_ConstKick_DAC_-2000_Interleaved',...
'20151123_1420_ConstKick_DAC_-1500_Interleaved',...
'20151123_1422_ConstKick_DAC_-1000_Interleaved',...
'20151123_1425_ConstKick_DAC_-500_Interleaved',...
'20151123_1428_ConstKick_DAC_0_Interleaved',...
'20151123_1430_ConstKick_DAC_500_Interleaved',...
'20151123_1433_ConstKick_DAC_1000_Interleaved',...
'20151123_1436_ConstKick_DAC_1500_Interleaved',...
'20151123_1438_ConstKick_DAC_2000_Interleaved',...
'20151123_1440_ConstKick_DAC_2500_Interleaved',...
'20151123_1443_ConstKick_DAC_3000_Interleaved',...
'20151123_1445_ConstKick_DAC_3500_Interleaved',...
'20151123_1452_ConstKick_DAC_4000_Interleaved'...
};
dataDir = '/home/jack/PhaseFeedforward/CTFData/201511';
saveDir = '/home/jack/PhaseFeedforward/Analysis/201511/Plots/20151123_1407_ConstKick';

ampSampRange = 270:330; % used to calculate means
ampFullPulseRange = [207:239 254:380 388:393]; % used to calculate deviations along pulse (range to be used after alignment)
delayL = 9;

phaseSampRange = 480:570;

ampInCounts = -4000:500:4000;

pointsToFit = 4:14;
%% load data
addpath('../../');

nDataSets = length(dataSetNames);

ampInVolts = ampInCounts*(2/4096);

% figure;
for ds=1:nDataSets
    load([dataDir '/' dataSetNames{ds} '/Merged/' dataSetNames{ds} '.mat']);
    
    % Initialise arrays
    if (ds==1)
        nPulses = length(CTFData)/2;
        
        ampNSamp = length(CTFData(1).CT_SCOPE02_CH05.Acquisition.value.value);
        kickOnLA = NaN(nDataSets,nPulses,ampNSamp);
        kickOnLB = NaN(nDataSets,nPulses,ampNSamp);
        kickOnRA = NaN(nDataSets,nPulses,ampNSamp);
        kickOnRB = NaN(nDataSets,nPulses,ampNSamp);
        kickOffLA = NaN(nDataSets,nPulses,ampNSamp);
        kickOffLB = NaN(nDataSets,nPulses,ampNSamp);
        kickOffRA = NaN(nDataSets,nPulses,ampNSamp);
        kickOffRB = NaN(nDataSets,nPulses,ampNSamp);
        ampSampInt = CTFData(1).CT_SCOPE02_CH05.Acquisition.sampleInterval.value;
        ampTimeAxis = (0:(ampNSamp-1)).*ampSampInt;
        
        phasNSamp = length(CTFData(1).CT_SCOPE01_CH06.Acquisition.value.value);
        phasSampInt = CTFData(1).CT_SCOPE01_CH06.Acquisition.sampleInterval.value;
        kickOnPhase = NaN(nDataSets,nPulses,phasNSamp);
        kickOffPhase = NaN(nDataSets,nPulses,phasNSamp);
    end
    
    % =============== Amplifier monitoring signals =======================
    
%     % ASSUMES SENSITIVITY THE SAME FOR ALL PULSES AND SIGNALS
%     sensitivity = CTFData(1).CT_SCOPE02_CH05.Acquisition.sensitivity.value;
%     ampLA =sensitivity.*double(extractCTFSignalFromMergedData('CT_SCOPE02_CH05.Acquisition.value.value',CTFData));
%     ampLB =sensitivity.*double(extractCTFSignalFromMergedData('CT_SCOPE02_CH06.Acquisition.value.value',CTFData));
%     ampRA =sensitivity.*double(extractCTFSignalFromMergedData('CT_SCOPE02_CH07.Acquisition.value.value',CTFData));
%     ampRB =sensitivity.*double(extractCTFSignalFromMergedData('CT_SCOPE02_CH08.Acquisition.value.value',CTFData));
    ampMon = extractAmpMon(CTFData);
    ampLA = squeeze(ampMon(1,:,:));
    ampLB = squeeze(ampMon(2,:,:));
    ampRA = squeeze(ampMon(3,:,:));
    ampRB = squeeze(ampMon(4,:,:));

    % delay amp L signals to align with amp R signals
    ampLA = delaySignal(ampLA,delayL);
    ampLB = delaySignal(ampLB,delayL);
    
%     % CONVERT MONITORING OUTPUT TO REAL AMPLIFIER OUTPUT (approx factor
%     % 115 with 12 dB added)
%     Now implemented in extractAmpMon
%     ampLA = ampLA*115*4;
%     ampLB = ampLB*115*4;
%     ampRA = ampRA*115*4;
%     ampRB = ampRB*115*4;
    
    % ===================== Phase ========================================    
    CTFData = removePulsesNoBeam(CTFData);
    
    [frascatiCalTimeStamp,~,~] = loadInitSettings([dataDir '/' dataSetNames{ds} '/initSettings.dat']);
    calDir = ['/home/jack/PhaseFeedforward/CTFData/' frascatiCalTimeStamp(1:6) '/FrascatiCalibrations'];
    [calibrationConstants, useDiode] = loadFrascatiCalibrationConstants(sprintf([calDir '/frascatiCalibrationConstants_' frascatiCalTimeStamp]));
    calAmp = calibrationConstants(3,1);
    calOff = calibrationConstants(3,4);

    [mix,dio]=extractMixerDiode(CTFData);
    mix = squeeze(mix(3,:,:));
    dio = squeeze(dio(3,:,:));
    phase = getPhaseMixerDiode(mix,dio,calAmp,calOff,useDiode);
    phase(1:2:end,:) = removeBadPulses(phase(1:2:end,:),phaseSampRange);
    phase(2:2:end,:) = removeBadPulses(phase(2:2:end,:),phaseSampRange);
    
    % =========== Split in to kick on and kick off =======================
    % define kick on as being set of even/odd pulses that has the max
    % absolute output on the monitoring signals.

    if (max(max(abs(ampLA(1:2:end,ampSampRange)))) > max(max(abs(ampLA(2:2:end,ampSampRange)))))        
        % odd is kick on
        kickOnLA(ds,:,:) = ampLA(1:2:end,:);
        kickOnLB(ds,:,:) = ampLB(1:2:end,:);
        kickOnRA(ds,:,:) = ampRA(1:2:end,:);
        kickOnRB(ds,:,:) = ampRB(1:2:end,:);

        kickOffLA(ds,:,:) = ampLA(2:2:end,:);
        kickOffLB(ds,:,:) = ampLB(2:2:end,:);
        kickOffRA(ds,:,:) = ampRA(2:2:end,:);
        kickOffRB(ds,:,:) = ampRB(2:2:end,:);
        
        kickOnPhase(ds,:,:) = phase(1:2:end,:);
        kickOffPhase(ds,:,:) = phase(2:2:end,:);
    else
        % even is kick on
        kickOnLA(ds,:,:) = ampLA(2:2:end,:);
        kickOnLB(ds,:,:) = ampLB(2:2:end,:);
        kickOnRA(ds,:,:) = ampRA(2:2:end,:);
        kickOnRB(ds,:,:) = ampRB(2:2:end,:);

        kickOffLA(ds,:,:) = ampLA(1:2:end,:);
        kickOffLB(ds,:,:) = ampLB(1:2:end,:);
        kickOffRA(ds,:,:) = ampRA(1:2:end,:);
        kickOffRB(ds,:,:) = ampRB(1:2:end,:);  
        
        kickOnPhase(ds,:,:) = phase(2:2:end,:);
        kickOffPhase(ds,:,:) = phase(1:2:end,:);
    end
    
%     subplot(2,2,1)
%     plot(ampLA(1,:));
%     hold all;
%     plot(ampLA(100,:));
%     hold off;
%     legend('ODD','EVEN')
%     title('LA')
%     
%     subplot(2,2,2)
%     plot(ampLB(1,:));
%     hold all;
%     plot(ampLB(100,:));
%     hold off;
%     legend('ODD','EVEN')
%     title('LB')    
%     
%     subplot(2,2,3)
%     plot(ampRA(1,:));
%     hold all;
%     plot(ampRA(100,:));
%     hold off;
%     legend('ODD','EVEN')
%     title('RA')   
%     
%     subplot(2,2,4)
%     plot(ampRB(1,:));
%     hold all;
%     plot(ampRB(100,:));
%     hold off;
%     legend('ODD','EVEN')
%     title('RB')
%     
%     plot(ampLA(1,:));
%     hold all
%     plot(ampLA(2,:));
%     legend('ODD','EVEN')
%     hold off;
%     input(sprintf('I think: %s',kickOnIs))

%     plot(abs(diff(ampLA(:,270)))>0.3)
%     input('next...')
    
    clear CTFData ampMon ampLA ampLB ampRA ampRB phase;
end

%% Amp monitoring: take differences, means
diffL = kickOnLA-kickOnLB;
diffR = kickOnRA-kickOnRB;

[sampleMeanLA,~,sampleMeanLA_err] = nanMeanStdErr(kickOnLA,2);
[sampleMeanLB,~,sampleMeanLB_err] = nanMeanStdErr(kickOnLB,2);
[sampleMeanRA,~,sampleMeanRA_err] = nanMeanStdErr(kickOnRA,2);
[sampleMeanRB,~,sampleMeanRB_err] = nanMeanStdErr(kickOnRB,2);

sampleMeanDiffL = sampleMeanLA-sampleMeanLB;
sampleMeanDiffR = sampleMeanRA-sampleMeanRB;
sampleMeanDiffL_err = sqrt(sampleMeanLA_err.^2 + sampleMeanLB_err.^2);
sampleMeanDiffR_err = sqrt(sampleMeanRA_err.^2 + sampleMeanRB_err.^2);

pulseMeanDiffL = squeeze(nanmean(diffL(:,:,ampSampRange),3));
pulseMeanDiffR = squeeze(nanmean(diffR(:,:,ampSampRange),3));
pulseMeanLA = squeeze(nanmean(kickOnLA(:,:,ampSampRange),3));
pulseMeanLB = squeeze(nanmean(kickOnLB(:,:,ampSampRange),3));
pulseMeanRA = squeeze(nanmean(kickOnRA(:,:,ampSampRange),3));
pulseMeanRB = squeeze(nanmean(kickOnRB(:,:,ampSampRange),3));

[meanKickL,~,meanKickL_err] = nanMeanStdErr(pulseMeanDiffL,2);
[meanKickR,~,meanKickR_err] = nanMeanStdErr(pulseMeanDiffR,2);
[meanKickLA,~,meanKickLA_err] = nanMeanStdErr(pulseMeanLA,2);
[meanKickLB,~,meanKickLB_err] = nanMeanStdErr(pulseMeanLB,2);
[meanKickRA,~,meanKickRA_err] = nanMeanStdErr(pulseMeanRA,2);
[meanKickRB,~,meanKickRB_err] = nanMeanStdErr(pulseMeanRB,2);

fprintf('OUTPUT @ MIN INPUT\n')
fprintf('AmpLA: %.2f %c %.2f V out per V in\n',meanKickLA(1),char(177),meanKickLA_err(1));
fprintf('AmpLB: %.2f %c %.2f V out per V in\n',meanKickLB(1),char(177),meanKickLB_err(1));
fprintf('AmpRA: %.2f %c %.2f V out per V in\n',meanKickRA(1),char(177),meanKickRA_err(1));
fprintf('AmpRB: %.2f %c %.2f V out per V in\n',meanKickRB(1),char(177),meanKickRB_err(1));

fprintf('\nOUTPUT @ MAX INPUT\n')
fprintf('AmpLA: %.2f %c %.2f V out per V in\n',meanKickLA(end),char(177),meanKickLA_err(end));
fprintf('AmpLB: %.2f %c %.2f V out per V in\n',meanKickLB(end),char(177),meanKickLB_err(end));
fprintf('AmpRA: %.2f %c %.2f V out per V in\n',meanKickRA(end),char(177),meanKickRA_err(end));
fprintf('AmpRB: %.2f %c %.2f V out per V in\n',meanKickRB(end),char(177),meanKickRB_err(end));

[fitMeanKickL,fitMeanKickL_rsq,fitMeanKickL_conf] = nanpolyfit(ampInVolts(pointsToFit),meanKickL(pointsToFit),1,1./meanKickL_err(pointsToFit).^2);
[fitMeanKickR,fitMeanKickR_rsq,fitMeanKickR_conf] = nanpolyfit(ampInVolts(pointsToFit),meanKickR(pointsToFit),1,1./meanKickR_err(pointsToFit).^2);
[fitMeanKickLA,fitMeanKickLA_rsq,fitMeanKickLA_conf] = nanpolyfit(ampInVolts(pointsToFit),meanKickLA(pointsToFit),1,1./meanKickLA_err(pointsToFit).^2);
[fitMeanKickLB,fitMeanKickLB_rsq,fitMeanKickLB_conf] = nanpolyfit(ampInVolts(pointsToFit),meanKickLB(pointsToFit),1,1./meanKickLB_err(pointsToFit).^2);
[fitMeanKickRA,fitMeanKickRA_rsq,fitMeanKickRA_conf] = nanpolyfit(ampInVolts(pointsToFit),meanKickRA(pointsToFit),1,1./meanKickRA_err(pointsToFit).^2);
[fitMeanKickRB,fitMeanKickRB_rsq,fitMeanKickRB_conf] = nanpolyfit(ampInVolts(pointsToFit),meanKickRB(pointsToFit),1,1./meanKickRB_err(pointsToFit).^2);

fitMeanKickL_err = (fitMeanKickL-fitMeanKickL_conf(1,:))/2;
fitMeanKickR_err = (fitMeanKickR-fitMeanKickR_conf(1,:))/2;
fitMeanKickLA_err = (fitMeanKickLA-fitMeanKickLA_conf(1,:))/2;
fitMeanKickLB_err = (fitMeanKickLB-fitMeanKickLB_conf(1,:))/2;
fitMeanKickRA_err = (fitMeanKickRA-fitMeanKickRA_conf(1,:))/2;
fitMeanKickRB_err = (fitMeanKickRB-fitMeanKickRB_conf(1,:))/2;


% [fitMeanKickL,fitMeanKickL_err] = linFit(ampInVolts(pointsToFit),meanKickL(pointsToFit),meanKickL_err(pointsToFit));
% [fitMeanKickR,fitMeanKickR_err] = linFit(ampInVolts(pointsToFit),meanKickR(pointsToFit),meanKickR_err(pointsToFit));
% [fitMeanKickLA,fitMeanKickLA_err] = linFit(ampInVolts(pointsToFit),meanKickLA(pointsToFit),meanKickLA_err(pointsToFit));
% [fitMeanKickLB,fitMeanKickLB_err] = linFit(ampInVolts(pointsToFit),meanKickLB(pointsToFit),meanKickLB_err(pointsToFit));
% [fitMeanKickRA,fitMeanKickRA_err] = linFit(ampInVolts(pointsToFit),meanKickRA(pointsToFit),meanKickRA_err(pointsToFit));
% [fitMeanKickRB,fitMeanKickRB_err] = linFit(ampInVolts(pointsToFit),meanKickRB(pointsToFit),meanKickRB_err(pointsToFit));

fprintf('\nFITTED GRADIENT\n')
fprintf('AmpLA: %.2f %c %.2f V out per V in\n',fitMeanKickLA(1),char(177),fitMeanKickLA_err(1));
fprintf('AmpLB: %.2f %c %.2f V out per V in\n',fitMeanKickLB(1),char(177),fitMeanKickLB_err(1));
fprintf('AmpRA: %.2f %c %.2f V out per V in\n',fitMeanKickRA(1),char(177),fitMeanKickRA_err(1));
fprintf('AmpRB: %.2f %c %.2f V out per V in\n',fitMeanKickRB(1),char(177),fitMeanKickRB_err(1));

fprintf('\nFITTED OFFSET\n')
fprintf('AmpLA: %.2f %c %.2f V out per V in\n',fitMeanKickLA(2),char(177),fitMeanKickLA_err(2));
fprintf('AmpLB: %.2f %c %.2f V out per V in\n',fitMeanKickLB(2),char(177),fitMeanKickLB_err(2));
fprintf('AmpRA: %.2f %c %.2f V out per V in\n',fitMeanKickRA(2),char(177),fitMeanKickRA_err(2));
fprintf('AmpRB: %.2f %c %.2f V out per V in\n',fitMeanKickRB(2),char(177),fitMeanKickRB_err(2));

residL = meanKickL - polyval(fitMeanKickL,ampInVolts)';
residR = meanKickR - polyval(fitMeanKickR,ampInVolts)';
residLA = meanKickLA - polyval(fitMeanKickLA,ampInVolts)';
residLB = meanKickLB - polyval(fitMeanKickLB,ampInVolts)';
residRA = meanKickRA - polyval(fitMeanKickRA,ampInVolts)';
residRB = meanKickRB - polyval(fitMeanKickRB,ampInVolts)';

% average deviation from mean along pulse
[wholeMeanL,flatnessDiffL] = nanMeanStdErr(sampleMeanDiffL(:,ampFullPulseRange),2);
[wholeMeanR,flatnessDiffR] = nanMeanStdErr(sampleMeanDiffR(:,ampFullPulseRange),2);
avgDevL = NaN(1,nDataSets);
avgDevR = NaN(1,nDataSets);
for i=1:nDataSets
    avgDevL(i) = nanmean(abs(sampleMeanDiffL(i,ampFullPulseRange)-wholeMeanL(i)));
    avgDevR(i) = nanmean(abs(sampleMeanDiffR(i,ampFullPulseRange)-wholeMeanR(i)));
end
peakDevL = abs(max(sampleMeanDiffL(:,ampFullPulseRange),[],2)-min(sampleMeanDiffL(:,ampFullPulseRange),[],2));
peakDevR = abs(max(sampleMeanDiffR(:,ampFullPulseRange),[],2)-min(sampleMeanDiffR(:,ampFullPulseRange),[],2));

% average deviation from mean of ampL+ampR = residual orbit closure error
sampleAmpSum = sampleMeanDiffL+sampleMeanDiffR;
ampSumMean = nanmean(sampleAmpSum(:,ampFullPulseRange),2);
avgDevAmpSum = NaN(1,nDataSets);
for i=1:nDataSets
    avgDevAmpSum(i) = nanmean(abs(sampleAmpSum(i,ampFullPulseRange)-ampSumMean(i)));
end
peakDevAmpSum = abs(max(sampleAmpSum(:,ampFullPulseRange),[],2)-min(sampleAmpSum(:,ampFullPulseRange),[],2));


%% Amp monitoring: plots
figure;
plot(ampTimeAxis,sampleMeanLA(13,:),'b','LineWidth',2);
hold all;
plot(ampTimeAxis,sampleMeanLB(13,:),'Color',[0 0.8 0],'LineWidth',2);
plot(ampTimeAxis,sampleMeanDiffL(13,:),'r','LineWidth',2);
title('AmpL')
xlabel('Time [ns]')
ylabel('Output [V]')
legend('A','B','Diff','Location','NorthWest')
format_plots;
% savePlot(saveDir,'AmpL_Traces');

figure;
plot(ampTimeAxis,sampleMeanRA(13,:),'b','LineWidth',2);
hold all;
plot(ampTimeAxis,sampleMeanRB(13,:),'Color',[0 0.8 0],'LineWidth',2);
plot(ampTimeAxis,sampleMeanDiffR(13,:),'r','LineWidth',2);
title('AmpR')
xlabel('Time [ns]')
ylabel('Residual [V]')
legend('A','B','Diff','Location','NorthWest')
format_plots;
% savePlot(saveDir,'AmpR_Traces');

figure;
plot(ampInVolts,avgDevL,'b','LineWidth',2);
hold all;
plot(ampInVolts,avgDevR,'r','LineWidth',2);
plot(ampInVolts,peakDevL,'b--','LineWidth',2);
plot(ampInVolts,peakDevR,'r--','LineWidth',2);
xlabel('Input [V]')
ylabel('Deviation from Flat [V]')
legend('L Avg','R Avg','L Peak','R Peak')
title('Voltage Variations Along Output Pulse')
format_plots;
% savePlot(saveDir,'AmpFlatness')

figure;
cols = varycolor(nDataSets);
for i=1:nDataSets
    plot(ampTimeAxis,sampleAmpSum(i,:),'LineWidth',2,'Color',cols(i,:));
    hold all;
end
title({'Reisdual Kick' '(AmpLA-AmpLB)+(AmpRA-AmpRB)'})
xlabel('Time [ns]')
ylabel('Output [V]')
format_plots;
set(gcf, 'Colormap', cols);
figColBar = colorbar;
allTicks = linspace(0,1,nDataSets+1);
set(figColBar,'YTick',allTicks(1:2:nDataSets));
set(figColBar,'YTickLabel',[-2.0 -1.5 -1.0 -0.5 0.0 0.5 1.0 1.5 2.0]);
ylabel(figColBar,'Input [V]');
% savePlot(saveDir,'residualKick_Traces');

figure;
plot(ampInVolts,avgDevAmpSum,'b','LineWidth',2);
hold all;
plot(ampInVolts,peakDevAmpSum,'b--','LineWidth',2);
xlabel('Input [V]')
ylabel('Deviation from Flat [V]')
legend('Avg','Peak')
title('Reisdual Kick: Deviation from Flat')
format_plots;
% savePlot(saveDir,'residualKick_Flatness')

figure; 
for i=1:nDataSets
    plot(ampInVolts(i),abs(avgDevAmpSum(i)./ampSumMean(i)),'bo','LineWidth',2);
    hold all;
    plot(ampInVolts(i),abs(peakDevAmpSum(i)./ampSumMean(i)),'bx','LineWidth',2);
end
xlabel('Input [V]')
ylabel('Relative Deviation from Flat')
title('Reisdual Kick: Relative Deviation from Flat')
format_plots;

figure;
plot(ampInVolts([1:8 10:17]),avgDevL([1:8 10:17])./abs(wholeMeanL([1:8 10:17])'),'b','LineWidth',2);
hold all;
plot(ampInVolts([1:8 10:17]),avgDevR([1:8 10:17])./abs(wholeMeanR([1:8 10:17])'),'r','LineWidth',2);
plot(ampInVolts([1:8 10:17]),peakDevL([1:8 10:17])./abs(wholeMeanL([1:8 10:17])),'b--','LineWidth',2);
plot(ampInVolts([1:8 10:17]),peakDevR([1:8 10:17])./abs(wholeMeanR([1:8 10:17])),'r--','LineWidth',2);
xlabel('Input [V]')
ylabel('Relative Deviation from Flat')
legend('L Avg','R Avg','L Peak','R Peak')
title('Voltage Variations Along Output Pulse')
format_plots;

% figure;
% plot(ampInVolts,flatnessDiffL);
% hold all;
% plot(ampInVolts,flatnessDiffR);
% 
% figure;
% plot(ampInVolts,flatnessDiffL./abs(wholeMeanL));
% hold all;
% plot(ampInVolts,flatnessDiffR./abs(wholeMeanR));


figure;
plot(ampInVolts,meanKickL,'bo','MarkerFaceColor','b');
hold all;
plot(ampInVolts,meanKickR,'ro','MarkerFaceColor','r');
legend(sprintf('AmpL (%.2f V per 1000cnts)',fitMeanKickL(1).*1000),sprintf('AmpR (%.2fV per 1000cnts)',fitMeanKickR(1).*1000),'Location','North');
plot(ampInVolts,polyval(fitMeanKickL,ampInVolts),'b');
plot(ampInVolts,polyval(fitMeanKickR,ampInVolts),'r');
xlabel('Input [V]')
ylabel('Output [V]')
format_plots
%savePlot(saveDir,'AmpOutvsDAC')

figure;
errorbar(ampInVolts,meanKickLA,meanKickLA_err,'bo','MarkerFaceColor','b');
hold all;
errorbar(ampInVolts,meanKickLB,meanKickLB_err,'ro','MarkerFaceColor','r');
errorbar(ampInVolts,meanKickRA,meanKickRA_err,'go','MarkerFaceColor','g');
errorbar(ampInVolts,meanKickRB,meanKickRB_err,'mo','MarkerFaceColor','m');
% legend(sprintf('AmpLA (%.0fV per 1000cnts)',fitMeanKickLA(1).*1000),...
%        sprintf('AmpLB (%.0fV per 1000cnts)',fitMeanKickLB(1).*1000),...
%        sprintf('AmpRA (%.0fV per 1000cnts)',fitMeanKickRA(1).*1000),...
%        sprintf('AmpRB (%.0fV per 1000cnts)',fitMeanKickRB(1).*1000),...       
%        'Location','North');
legend('AmpLA','AmpLB','AmpRA','AmpRB','Location','Best')
plot(ampInVolts,polyval(fitMeanKickLA,ampInVolts),'b');
plot(ampInVolts,polyval(fitMeanKickLB,ampInVolts),'r');
plot(ampInVolts,polyval(fitMeanKickRA,ampInVolts),'g');
plot(ampInVolts,polyval(fitMeanKickRB,ampInVolts),'m');
xlabel('Input [V]')
ylabel('Output [V]')
format_plots
xlim([-2 2])
ylim([-800 800])
title('Amplifier Output vs. Input')
% savePlot(saveDir,'AmpOutvsDAC_allStrips')


figure;
plot(ampInVolts,residLA,'bo-','MarkerFaceColor','b');
hold all;
plot(ampInVolts,residLB,'ro-','MarkerFaceColor','r');
plot(ampInVolts,residRA,'go-','MarkerFaceColor','g');
plot(ampInVolts,residRB,'mo-','MarkerFaceColor','m');
% legend(sprintf('AmpLA (%.0fV per 1000cnts)',fitMeanKickLA(1).*1000),...
%        sprintf('AmpLB (%.0fV per 1000cnts)',fitMeanKickLB(1).*1000),...
%        sprintf('AmpRA (%.0fV per 1000cnts)',fitMeanKickRA(1).*1000),...
%        sprintf('AmpRB (%.0fV per 1000cnts)',fitMeanKickRB(1).*1000),...       
%        'Location','North');
legend('AmpLA','AmpLB','AmpRA','AmpRB','Location','Best')
xlabel('Input [V]')
ylabel('Difference to Fit [V]')
format_plots
xlim([-2 2])
title('Difference to Linear Fit')
% savePlot(saveDir,'AmpOutvsDAC_allStrips_residual')

figure;
plot(ampInVolts,residLA./meanKickLA,'bo-','MarkerFaceColor','b');
hold all;
plot(ampInVolts,residLB./meanKickLB,'ro-','MarkerFaceColor','r');
plot(ampInVolts,residRA./meanKickRA,'go-','MarkerFaceColor','g');
plot(ampInVolts,residRB./meanKickRB,'mo-','MarkerFaceColor','m');
% legend(sprintf('AmpLA (%.0fV per 1000cnts)',fitMeanKickLA(1).*1000),...
%        sprintf('AmpLB (%.0fV per 1000cnts)',fitMeanKickLB(1).*1000),...
%        sprintf('AmpRA (%.0fV per 1000cnts)',fitMeanKickRA(1).*1000),...
%        sprintf('AmpRB (%.0fV per 1000cnts)',fitMeanKickRB(1).*1000),...       
%        'Location','North');
legend('AmpLA','AmpLB','AmpRA','AmpRB','Location','Best')
xlabel('Input [V]')
ylabel('Relative Difference to Fit')
format_plots
xlim([-2 2])
%savePlot(saveDir,'AmpOutvsDAC')

% 
% figure;
% plot(kick,6*meanKickL/1300,'bo-','MarkerFaceColor','b');
% hold all;
% plot(kick,6*meanKickR/1300,'ro-','MarkerFaceColor','r');
% title('very rough expected phase shift')
% 
% figure;
% plot(kick,6*residL/1300,'bo-','MarkerFaceColor','b');
% hold all;
% plot(kick,6*residR/1300,'ro-','MarkerFaceColor','r');
% title('very rough aprox phase residual')

% % AmpOutvsDAC in scale equivalent to one strip
% figure;
% plot(kick,meanKickL*230,'bo','MarkerFaceColor','b');
% hold all;
% plot(kick,meanKickR*230,'ro','MarkerFaceColor','r');
% legend(sprintf('AmpL (%.0f V per 1000cnts)',fitMeanKickL(1).*1000*230),sprintf('AmpR (%.0fV per 1000cnts)',fitMeanKickR(1).*1000*230),'Location','North');
% plot(kick,polyval(fitMeanKickL,kick)*230,'b');
% plot(kick,polyval(fitMeanKickR,kick)*230,'r');
% xlabel('DAC1 Output [counts]')
% ylabel('Amplifier Monitoring Output [V]')
% format_plots
% %savePlot(saveDir,'AmpOutvsDAC')

%% Phase: calculate means etc.

% diffPhase = kickOnPhase-kickOffPhase;
% phaseDiffSampMean = squeeze(nanmean(diffPhase,2));
% phaseDiffPulseMean = squeeze(nanmean(diffPhase(:,:,phaseSampRange),3));
% [phaseDiffMean,~,phaseDiffMean_err] = nanMeanStdErr(phaseDiffPulseMean,2);

[phaseOnSampMean,~,phaseOnSampMean_err] = nanMeanStdErr(kickOnPhase,2);
phaseOnPulseMean = squeeze(nanmean(kickOnPhase(:,:,phaseSampRange),3));
[phaseOnMean,~,phaseOnMean_err] = nanMeanStdErr(phaseOnPulseMean,2);

[phaseOffSampMean,~,phaseOffSampMean_err] = nanMeanStdErr(kickOffPhase,2);
phaseOffPulseMean = squeeze(nanmean(kickOffPhase(:,:,phaseSampRange),3));
[phaseOffMean,~,phaseOffMean_err] = nanMeanStdErr(phaseOffPulseMean,2);

phaseDiffSampMean = phaseOnSampMean-phaseOffSampMean;
phaseDiffSampMean_err = sqrt(phaseOnSampMean_err.^2 + phaseOffSampMean_err.^2);
phaseDiffPulseMean = phaseOnPulseMean-phaseOffPulseMean;
phaseDiffMean = phaseOnMean-phaseOffMean;
phaseDiffMean_err = sqrt(phaseOnMean_err.^2 + phaseOffMean_err.^2);


[phaseDiffFit,phaseDiffFit_rsq,phaseDiffFit_conf] = nanpolyfit(ampInVolts(pointsToFit),phaseDiffMean(pointsToFit),1,1./phaseDiffMean_err(pointsToFit).^2);
phaseDiffFit_err = (phaseDiffFit-phaseDiffFit_conf(1,:))/2;

%% Call MADX model to get expected phase shift from amplifier output voltages.

myPath ='/home/jack/Documents/MATLAB/PhaseFeedforward/ctfmod/PhaseFeedforward/R56/MADXv2';
addpath('/home/jack/Documents/MATLAB/madx2Matlab');
origPath = pwd;
cd(myPath);

machineModel = madx2matlab(fullfile(myPath,'tl2From400.madx'),...
    '',...
    fullfile('/usr/local/bin/madx'));

machineModel.twissCommandLine=[...
    'select, flag=twiss, clear;',char(10),...
    'select, flag=twiss, column=name, s, BETX, ALFX, L, ANGLE, x, px, dx, BETY, ALFY, y, py, dy,mux,muy, t, pt,re56, t566, re52;',char(10),...
    'twiss, RMATRIX,',char(10),...
    'BETX=initialBETX,ALFX=initialALFX,',char(10),...
    'BETY=initialBETY,ALFY=initialALFY,',char(10),...
    'X=0,Y=0,',char(10),...
    'DELTAP=initialDeltaP;',char(10),...
    ];

% init parameters found in old file for exit of vertical chicane
initialBETX =    4.060520895;
initialBETY =    3.9900588;
initialALFX =    -0.7650527317;
initialALFY =    -5.974166392;
initialX = 0;
initialY = 0;
initialDeltaP = 0.0;
machineModel.setModelParameters('initialBETX', initialBETX);
machineModel.setModelParameters('initialBETY', initialBETY);
machineModel.setModelParameters('initialALFX', initialALFX);
machineModel.setModelParameters('initialALFY', initialALFY);
machineModel.setModelParameters('initialX', initialX);
machineModel.setModelParameters('initialY', initialY);
machineModel.setModelParameters('initialDELTAP', initialDeltaP);

machineModel.setModelParameters('cc.kcktot', 0.001); % angle deflection from both kickers

nominalOptic = machineModel.computeOptics();

figure;
subplot(2,2,1)
plot([nominalOptic.DATA.S], [nominalOptic.DATA.DX])
title('DX')
subplot(2,2,2)
plot([nominalOptic.DATA.S], [nominalOptic.DATA.RE56]);
title('R56')
subplot(2,2,3)
plot([nominalOptic.DATA.S], (360/0.025)*[nominalOptic.DATA.T]);
title('phase')
subplot(2,2,4)
plot([nominalOptic.DATA.S], [nominalOptic.DATA.X]);
title('X')

meanKickAng = (meanKickL/(2*1260))*0.001; %1.26 kV to each strip = 1 mrad kick at 135 MeV, meanKickL is potential difference -> meanKick/2 is mean absolute each strip first kicker
meanModelPhase = NaN(1,nDataSets);
for i=1:nDataSets

    machineModel.setModelParameters('initialBETX', initialBETX);
    machineModel.setModelParameters('initialBETY', initialBETY);
    machineModel.setModelParameters('initialALFX', initialALFX);
    machineModel.setModelParameters('initialALFY', initialALFY);
    machineModel.setModelParameters('initialX', initialX);
    machineModel.setModelParameters('initialY', initialY);
    machineModel.setModelParameters('initialDELTAP', initialDeltaP);

    machineModel.setModelParameters('cc.kcktot', meanKickAng(i));

    auxOptic = machineModel.computeOptics();
    auxPhas = [auxOptic.DATA.T];
    meanModelPhase(i) = (360/0.025)*auxPhas(end);
end

meanModelPhase = -meanModelPhase;

[fitModel,fitModel_rsq,fitModel_conf] = nanpolyfit(ampInVolts(pointsToFit),meanModelPhase(pointsToFit),1);

% figure;
% plot(meanKickL/2,meanModelPhase);
% 
cd(origPath);
%% Phase: plots

figure;
errorbar(ampInVolts,phaseOffMean,phaseOffMean_err,'bo','MarkerFaceColor','b');
hold all;
errorbar(ampInVolts,phaseOnMean,phaseOnMean_err,'ro','MarkerFaceColor','r');
legend('KickOff','KickOn')
xlabel('Input [V]')
ylabel('Phase [degrees]')
title('Phase Shift vs. Amplifier Input Voltage')
xlim([-2 2])
format_plots;


figure;
errorbar(ampInVolts,phaseDiffMean,phaseDiffMean_err,'bo','MarkerFaceColor','b');
hold all;
plot(ampInVolts,polyval(phaseDiffFit,ampInVolts),'b','LineWidth',1.5)
plot(ampInVolts,meanModelPhase,'r','LineWidth',2)
legend('DATA','FIT','MODEL')
xlabel('Input [V]')
ylabel('Phase [degrees]')
title('Phase Shift vs. Amplifier Input Voltage')
xlim([-2 2])
format_plots;
% savePlot(saveDir,'phaseVsAmpVoltage')

figure;
plot([-2 2],[0 0],'k','LineWidth',2)
hold all;
residErr = sqrt((ampInVolts'.*phaseDiffFit_err(1)).^2 + phaseDiffMean_err.^2);
a=errorbar(ampInVolts',phaseDiffMean-polyval(phaseDiffFit,ampInVolts)',residErr,'bo-','MarkerFaceColor','b');
b=errorbar(ampInVolts,phaseDiffMean'-meanModelPhase,phaseDiffMean_err','ro-','MarkerFaceColor','r');
legend([a,b],'FIT','MODEL')
ylabel('Residual [degrees]')
xlabel('Input [V]')
title('Residuals Between Phase Shift and Fit/Model')
format_plots;
% savePlot(saveDir,'phaseVsAmpVoltage_residuals')

% %% phase - differences in kick along pulse - SLOW (thousands of MADX calls)!!!
% 
% origPath = pwd;
% cd(myPath);
% sampKickAng = (sampleMeanDiffL/(2*1260))*0.001; %1.26 kV to each strip = 1 mrad kick at 135 MeV, meanKickL is potential difference -> meanKick/2 is mean absolute each strip first kicker
% sampModelPhase = NaN(nDataSets,ampNSamp);
% fprintf('CALLING MADX...\n')
% for ds=1:nDataSets
%     fprintf('Dataset %d of %d...\n',ds,nDataSets);
%     for i=1:ampNSamp
% 
%         machineModel.setModelParameters('initialBETX', initialBETX);
%         machineModel.setModelParameters('initialBETY', initialBETY);
%         machineModel.setModelParameters('initialALFX', initialALFX);
%         machineModel.setModelParameters('initialALFY', initialALFY);
%         machineModel.setModelParameters('initialX', initialX);
%         machineModel.setModelParameters('initialY', initialY);
%         machineModel.setModelParameters('initialDELTAP', initialDeltaP);
% 
%         machineModel.setModelParameters('cc.kcktot', sampKickAng(ds,i));
% 
%         auxOptic = machineModel.computeOptics();
%         auxPhas = [auxOptic.DATA.T];
%         sampModelPhase(ds,i) = (360/0.025)*auxPhas(end);
%     end
% end
% sampModelPhase = -sampModelPhase;
% sampModelPhase_err = (sampleMeanDiffL_err./abs(sampleMeanDiffL)).*abs(sampModelPhase);
% 
% cd(origPath);
% 
% %% phase - differences in kick along pulse - plots
% % figure;
% % for i=1:nDataSets
% %     errorbar(phaseOnSampMean(i,:),phaseOnSampMean_err(i,:)); 
% %     hold all;
% %     errorbar(phaseOffSampMean(i,:),phaseOffSampMean_err(i,:)); 
% %     errorbar(phaseDiffSampMean(i,:),phaseDiffSampMean_err(i,:));
% %     hold off;
% %     title(sprintf('%.2f',ampInVolts(i)))
% %     input('next...')
% % end
% 
% % figure;
% % errorbar(phaseDiffSampMean',phaseDiffSampMean_err');
% 
% figure;
% for i=1:nDataSets
%     errorbar(phasSampInt*(1:phasNSamp),phaseDiffSampMean(i,:),phaseDiffSampMean_err(i,:),'b');
%     hold all;
%     errorbar(ampSampInt*(1:ampNSamp),sampModelPhase(i,:),sampModelPhase_err(i,:),'r'); 
%     hold off;
%     title({'Phase Shift Along Pulse (NOT ALIGNED)', sprintf('Amp Input: %.2f V',ampInVolts(i))})
%     xlabel('Time [ns]')
%     ylabel('Phase [degrees]')
%     legend('DATA','MODEL')
%     xlim([1000 4000])
%     format_plots;
%     savePlot([saveDir '/alongPulse/'],sprintf('ampV_%.2f',ampInVolts(i)));
% end
