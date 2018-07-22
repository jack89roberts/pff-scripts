% orbitClosure
%
% Takes constant kick data and makes some orbit closure plots

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
savePlots = 0;

ampSampRange = 270:330; % used to calculate means
ampFullPulseRange = [207:239 254:380 388:393]; % used to calculate deviations along pulse (range to be used after alignment)
delayL = 9;

bpmSampRange = 555:698;

ampInCounts = -4000:500:4000;

pointsToFit = 4:14;

bpmNames = {...
    'CC_SVBPM0365',...
    'CC_SVBPM0435',...
    'CC_SVBPI0535',...
    'CC_SVBPI0645',...
    'CC_SVBPI0685',...
    'CC_SVBPI0735',...
    'CC_SVBPM0845',...
    'CC_SVBPM0930'...
};
bpmPos = [3.1831 6.2205 10.1413 13.4795 15.5532 17.0922 20.8066 25.3145];

%% load data
addpath('../../');

nDataSets = length(dataSetNames);
nBPMs = length(bpmNames);

ampInVolts = ampInCounts*(2/4096);

% figure;
for ds=1:nDataSets
    load([dataDir '/' dataSetNames{ds} '/Merged/' dataSetNames{ds} '.mat']);
    
    % Initialise arrays
    if (ds==1)
        nPulses = length(CTFData)/2; % divide by 2 - interleaved mode, splitting pulses in to different arrays later
        
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
        
        bpmNSamp = length(CTFData(1).CC_SVBPI0645V.Samples.samples.value);
        bpmSampInt = 5.2083;
        kickOnBPMH = NaN(nDataSets,nBPMs,nPulses,bpmNSamp);
        kickOffBPMH = NaN(nDataSets,nBPMs,nPulses,bpmNSamp);
        kickOnBPMV = NaN(nDataSets,nBPMs,nPulses,bpmNSamp);
        kickOffBPMV = NaN(nDataSets,nBPMs,nPulses,bpmNSamp);
        kickOnBPMS = NaN(nDataSets,nBPMs,nPulses,bpmNSamp);
        kickOffBPMS = NaN(nDataSets,nBPMs,nPulses,bpmNSamp);

    end
    
    % =============== Amplifier monitoring signals =======================
    
    % ASSUMES SENSITIVITY THE SAME FOR ALL PULSES AND SIGNALS
    sensitivity = CTFData(1).CT_SCOPE02_CH05.Acquisition.sensitivity.value;
    ampLA =sensitivity.*double(extractCTFSignalFromMergedData('CT_SCOPE02_CH05.Acquisition.value.value',CTFData));
    ampLB =sensitivity.*double(extractCTFSignalFromMergedData('CT_SCOPE02_CH06.Acquisition.value.value',CTFData));
    ampRA =sensitivity.*double(extractCTFSignalFromMergedData('CT_SCOPE02_CH07.Acquisition.value.value',CTFData));
    ampRB =sensitivity.*double(extractCTFSignalFromMergedData('CT_SCOPE02_CH08.Acquisition.value.value',CTFData));
    
    % delay amp L signals to align with amp R signals
    ampLA = delaySignal(ampLA,delayL);
    ampLB = delaySignal(ampLB,delayL);
    
    % CONVERT MONITORING OUTPUT TO REAL AMPLIFIER OUTPUT (approx factor
    % 115 with 12 dB added)
    ampLA = ampLA*115*4;
    ampLB = ampLB*115*4;
    ampRA = ampRA*115*4;
    ampRB = ampRB*115*4;
    
    % ===================== BPM ========================================    
    CTFData = removePulsesNoBeam(CTFData);
    
    bpmH = NaN(nBPMs,nPulses*2,bpmNSamp);
    bpmS = NaN(nBPMs,nPulses*2,bpmNSamp);
    bpmV = NaN(nBPMs,nPulses*2,bpmNSamp);
    for b=1:nBPMs
        [bpmH(b,:,:), bpmS(b,:,:), bpmV(b,:,:)] = extractBPMFromCTFData( bpmNames{b}, CTFData);
        [bpmS(b,:,:),others] = getAlignedXCorr(squeeze(bpmS(b,:,:)),'end',{squeeze(bpmH(b,:,:)) squeeze(bpmV(b,:,:))});
        bpmH(b,:,:) = others{1};
        bpmV(b,:,:) = others{2};
        
        [bpmS(b,1:2:end,:),others] = removeBadPulses(squeeze(bpmS(b,1:2:end,:)),bpmSampRange,{squeeze(bpmH(b,1:2:end,:)) squeeze(bpmV(b,1:2:end,:))});
        bpmH(b,1:2:end,:) = others{1};
        bpmV(b,1:2:end,:) = others{2};
        
        [bpmS(b,2:2:end,:),others] = removeBadPulses(squeeze(bpmS(b,2:2:end,:)),bpmSampRange,{squeeze(bpmH(b,2:2:end,:)) squeeze(bpmV(b,2:2:end,:))});
        bpmH(b,2:2:end,:) = others{1};
        bpmV(b,2:2:end,:) = others{2};

    end
    
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
        
        kickOnBPMH(ds,:,:,:) = bpmH(:,1:2:end,:);
        kickOnBPMV(ds,:,:,:) = bpmV(:,1:2:end,:);
        kickOnBPMS(ds,:,:,:) = bpmS(:,1:2:end,:);
        
        kickOffBPMH(ds,:,:,:) = bpmH(:,2:2:end,:);
        kickOffBPMV(ds,:,:,:) = bpmV(:,2:2:end,:);
        kickOffBPMS(ds,:,:,:) = bpmS(:,2:2:end,:);

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
        
        kickOnBPMH(ds,:,:,:) = bpmH(:,2:2:end,:);
        kickOnBPMV(ds,:,:,:) = bpmV(:,2:2:end,:);
        kickOnBPMS(ds,:,:,:) = bpmS(:,2:2:end,:);
        
        kickOffBPMH(ds,:,:,:) = bpmH(:,1:2:end,:);
        kickOffBPMV(ds,:,:,:) = bpmV(:,1:2:end,:);
        kickOffBPMS(ds,:,:,:) = bpmS(:,1:2:end,:);

    end
        
    clear CTFData ampLA ampLB ampRA ampRB bpmH bpmV bpmS;
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



%% BPMs: calculate means etc.

[bpmHOnSampMean,~,bpmHOnSampMean_err] = nanMeanStdErr(kickOnBPMH,3);
[bpmSOnSampMean,~,bpmSOnSampMean_err] = nanMeanStdErr(kickOnBPMS,3);
[bpmVOnSampMean,~,bpmVOnSampMean_err] = nanMeanStdErr(kickOnBPMV,3);
[bpmHOffSampMean,~,bpmHOffSampMean_err] = nanMeanStdErr(kickOffBPMH,3);
[bpmSOffSampMean,~,bpmSOffSampMean_err] = nanMeanStdErr(kickOffBPMS,3);
[bpmVOffSampMean,~,bpmVOffSampMean_err] = nanMeanStdErr(kickOffBPMV,3);

bpmHOnPulseMean = squeeze(nanmean(kickOnBPMH(:,:,:,bpmSampRange),4));
bpmVOnPulseMean = squeeze(nanmean(kickOnBPMV(:,:,:,bpmSampRange),4));
bpmSOnPulseMean = squeeze(nanmean(kickOnBPMS(:,:,:,bpmSampRange),4));
bpmHOffPulseMean = squeeze(nanmean(kickOffBPMH(:,:,:,bpmSampRange),4));
bpmVOffPulseMean = squeeze(nanmean(kickOffBPMV(:,:,:,bpmSampRange),4));
bpmSOffPulseMean = squeeze(nanmean(kickOffBPMS(:,:,:,bpmSampRange),4));

[bpmHOnMean,~,bpmHOnMean_err] = nanMeanStdErr(bpmHOnPulseMean,3);
[bpmVOnMean,~,bpmVOnMean_err] = nanMeanStdErr(bpmVOnPulseMean,3);
[bpmSOnMean,~,bpmSOnMean_err] = nanMeanStdErr(bpmSOnPulseMean,3);
[bpmHOffMean,~,bpmHOffMean_err] = nanMeanStdErr(bpmHOffPulseMean,3);
[bpmVOffMean,~,bpmVOffMean_err] = nanMeanStdErr(bpmVOffPulseMean,3);
[bpmSOffMean,~,bpmSOffMean_err] = nanMeanStdErr(bpmSOffPulseMean,3);

bpmHDiffSampMean = bpmHOnSampMean-bpmHOffSampMean;
bpmHDiffSampMean_err = sqrt(bpmHOnSampMean_err.^2 + bpmHOffSampMean_err.^2);
bpmHDiffPulseMean = bpmHOnPulseMean-bpmHOffPulseMean;
bpmHDiffMean = bpmHOnMean-bpmHOffMean;
bpmHDiffMean_err = sqrt(bpmHOnMean_err.^2 + bpmHOffMean_err.^2);

bpmVDiffSampMean = bpmVOnSampMean-bpmVOffSampMean;
bpmVDiffSampMean_err = sqrt(bpmVOnSampMean_err.^2 + bpmVOffSampMean_err.^2);
bpmVDiffPulseMean = bpmVOnPulseMean-bpmVOffPulseMean;
bpmVDiffMean = bpmVOnMean-bpmVOffMean;
bpmVDiffMean_err = sqrt(bpmVOnMean_err.^2 + bpmVOffMean_err.^2);

bpmSDiffSampMean = bpmSOnSampMean-bpmSOffSampMean;
bpmSDiffSampMean_err = sqrt(bpmSOnSampMean_err.^2 + bpmSOffSampMean_err.^2);
bpmSDiffPulseMean = bpmSOnPulseMean-bpmSOffPulseMean;
bpmSDiffMean = bpmSOnMean-bpmSOffMean;
bpmSDiffMean_err = sqrt(bpmSOnMean_err.^2 + bpmSOffMean_err.^2);


bpmHDiffFit = NaN(nBPMs,2);
bpmHDiffFit_rsq = NaN(1,nBPMs);
bpmHDiffFit_err = NaN(nBPMs,2);

bpmVDiffFit = NaN(nBPMs,2);
bpmVDiffFit_rsq = NaN(1,nBPMs);
bpmVDiffFit_err = NaN(nBPMs,2);

bpmSDiffFit = NaN(nBPMs,2);
bpmSDiffFit_rsq = NaN(1,nBPMs);
bpmSDiffFit_err = NaN(nBPMs,2);
for b=1:nBPMs
    [bpmHDiffFit(b,:),bpmHDiffFit_rsq(b),tmpConf] = nanpolyfit(ampInVolts(pointsToFit),bpmHDiffMean(pointsToFit,b),1,1./bpmHDiffMean_err(pointsToFit,b).^2);
    bpmHDiffFit_err(b,:) = (bpmHDiffFit(b,:)-tmpConf(1,:))/2;
    
    [bpmVDiffFit(b,:),bpmVDiffFit_rsq(b),tmpConf] = nanpolyfit(ampInVolts(pointsToFit),bpmVDiffMean(pointsToFit,b),1,1./bpmVDiffMean_err(pointsToFit,b).^2);
    bpmVDiffFit_err(b,:) = (bpmVDiffFit(b,:)-tmpConf(1,:))/2;
    
    [bpmSDiffFit(b,:),bpmSDiffFit_rsq(b),tmpConf] = nanpolyfit(ampInVolts(pointsToFit),bpmSDiffMean(pointsToFit,b),1,1./bpmSDiffMean_err(pointsToFit,b).^2);
    bpmSDiffFit_err(b,:) = (bpmSDiffFit(b,:)-tmpConf(1,:))/2;
end

%%
xTickLab = strrep(bpmNames,'_','.');
xTickLab = strrep(xTickLab,'SV','');
xTickLab = strrep(xTickLab,'BPM','');
xTickLab = strrep(xTickLab,'BPI','');
xTickLab = strrep(xTickLab,'.0','.');

figure
myCols = varycolor(nDataSets);
for i=1:nDataSets
    %shadedErrorBar(1:nBPMs,bpmHDiffMean(i,:),bpmHDiffMean_err(i,:),{'Color',myCols(i,:)},0.9);
    plot(1:nBPMs,bpmHDiffMean(i,:),'LineWidth',2,'Color',myCols(i,:),'LineWidth',2);
    hold all;
end
set(gca,'XTickLabel',xTickLab)
set(gca,'XTick',1:8)
xlim([1 8])
ylim([-1.5 1.5])
plot([2.5 2.5],[-1.5 1.5],'k');
plot([6.5 6.5],[-1.5 1.5],'k');
rotateXLabels(gca,45)
xlabel('BPM Index')
ylabel('Position [mm]')
title('Horizontal Position')
set(gcf, 'Colormap', myCols);
figColBar = colorbar;
allTicks = linspace(0,1,nDataSets+1);
set(figColBar,'YTick',allTicks(1:2:nDataSets));
set(figColBar,'YTickLabel',[-2.0 -1.5 -1.0 -0.5 0.0 0.5 1.0 1.5 2.0]);
ylabel(figColBar,'Input [V]');
format_plots;
if savePlots
    savePlot(saveDir,'HOrbitVsInput');
end;


figure
myCols = varycolor(nDataSets);
for i=1:nDataSets
    shadedErrorBar(1:nBPMs,bpmVDiffMean(i,:),bpmVDiffMean_err(i,:),{'Color',myCols(i,:)},0.9);
%     plot(1:nBPMs,bpmHDiffMean(i,:),'LineWidth',2);%,'Color',myCols(i,:),'LineWidth',2);
    hold all;
end
set(gca,'XTickLabel',xTickLab)
set(gca,'XTick',1:8)
xlim([1 8])
ylim([-1.5 1.5])
plot([2.5 2.5],[-1.5 1.5],'k');
plot([6.5 6.5],[-1.5 1.5],'k');
rotateXLabels(gca,45)
xlabel('BPM Index')
ylabel('Position [mm]')
title('Vertical Position')
set(gcf, 'Colormap', myCols);
figColBar = colorbar;
allTicks = linspace(0,2,nDataSets+1);
set(figColBar,'YTick',allTicks(1:2:nDataSets));
set(figColBar,'YTickLabel',[-2.0 -1.5 -1.0 -0.5 0.0 0.5 1.0 1.5 2.0]);
ylabel(figColBar,'Input [V]');
format_plots;
if savePlots; savePlot(saveDir,'VOrbitVsInput'); end;

% figure
% myCols = varycolor(nDataSets);
% for i=1:nDataSets
%     shadedErrorBar(1:nBPMs,bpmSDiffMean(i,:),bpmSDiffMean_err(i,:),{'Color',myCols(i,:)});
% %     plot(1:nBPMs,bpmHDiffMean(i,:),'Color',myCols(i,:),'LineWidth',2);
%     hold all;
% end
% 
% figure;
% for b=1:nBPMs
%     plot(ampInVolts,bpmHDiffMean(:,b));
%     hold all;
% end
% figure;
% for b=1:nBPMs
%     errorbar(ampInVolts,bpmVDiffMean(:,b),bpmVDiffMean_err(:,b));
%     hold all;
% end
% figure;
% for b=1:nBPMs
%     errorbar(ampInVolts,bpmSDiffMean(:,b),bpmSDiffMean_err(:,b));
%     hold all;
% end

%%
%% Call MADX model to get expected orbit from amplifier output voltages.

myPath ='/home/jack/Documents/MATLAB/PhaseFeedforward/ctfmod/PhaseFeedforward/R56/MADXv2';
addpath('/home/jack/Documents/MATLAB/madx2Matlab');
origPath = pwd;
cd(myPath);

machineModel = madx2matlab(fullfile(myPath,'tl2From300.madx'),...
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

% init parameters found in old file for entrance to vertical chicane
initialBETX =    4.730253952;
initialBETY =    0.1000845332;
initialALFX =   1.652390883;
initialALFY =   -2.83608518;
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

 %1.26 kV to each strip = 1 mrad kick at 135 MeV, meanKickL is potential difference -> meanKick/2 is mean absolute each strip first kicker
meanKickAng = (fitMeanKickL(1)/(2*1260))*0.001;
machineModel.setModelParameters('cc.kcktot', meanKickAng); % angle deflection from both kickers

nominalOptic = machineModel.computeOptics();
nomModX = [nominalOptic.DATA.X].*1000;
nomModX_err = nomModX.*(fitMeanKickL_err(1)./fitMeanKickL(1));

% real amplifier voltage ratio
machineModel.setModelParameters('cc.kckratio', abs(fitMeanKickR(1)./fitMeanKickL(1)));

ampRatOptic = machineModel.computeOptics();
ampRatModX = [ampRatOptic.DATA.X].*1000;
ampRatModX_err = ampRatModX.*(fitMeanKickL_err(1)./fitMeanKickL(1));


% real quad currents taken from initconditions file (scaled to 140 MeV)
machineModel.setModelParameters('cc.kckratio', 1);
machineModel.setModelParameters('CC.IQDL0330', 24.173000*(140/135));
machineModel.setModelParameters('CC.IQFH0350', 44.312000*(140/135));
machineModel.setModelParameters('CC.IQDL0430', 0.153000*(140/135));
machineModel.setModelParameters('CC.IQFL0450', 0.000000*(140/135));
machineModel.setModelParameters('CC.IQDH0490', -5.704000*(140/135));
machineModel.setModelParameters('CC.IQDL0550', 28.125000*(140/135));
machineModel.setModelParameters('CC.IQFL0530', 21.628000*(140/135));
machineModel.setModelParameters('CC.IQFL0570', 30.343000*(140/135));
machineModel.setModelParameters('CC.IQFL0620', -1.151000*(140/135));
machineModel.setModelParameters('CC.IQDL0650', -10.407000*(140/135));
machineModel.setModelParameters('CC.IQFL0680', -5.919000*(140/135));
machineModel.setModelParameters('CC.IQFL0730', -6.941000*(140/135));
machineModel.setModelParameters('CC.IQDL0750', -14.540000*(140/135));
machineModel.setModelParameters('CC.IQDH0790', 6.293000*(140/135));
machineModel.setModelParameters('CC.IQDD0820', 60.266998*(140/135));
machineModel.setModelParameters('CC.IQFD0840', 41.581001*(140/135));
machineModel.setModelParameters('CC.IQFL0910', -12.694000*(140/135));
machineModel.setModelParameters('CC.IQDL0920', -10.872000*(140/135));

realQuadOptic = machineModel.computeOptics();
rqModX = [realQuadOptic.DATA.X].*1000;
rqModX_err = nomModX.*(fitMeanKickL_err(1)./fitMeanKickL(1));

cd(origPath);

figure;
subplot(2,2,1)
plot([nominalOptic.DATA.S], [nominalOptic.DATA.DX])
hold all;
plot([realQuadOptic.DATA.S], [realQuadOptic.DATA.DX])
title('DX')
legend('NOMNIAL','ACTUALQUADS')
subplot(2,2,2)
plot([nominalOptic.DATA.S], [nominalOptic.DATA.RE56]);
hold all;
plot([realQuadOptic.DATA.S], [realQuadOptic.DATA.RE56]);
title('R56')
legend('NOMNIAL','ACTUALQUADS')
subplot(2,2,3)
plot([nominalOptic.DATA.S], (360/0.025)*[nominalOptic.DATA.T]);
hold all;
plot([realQuadOptic.DATA.S], (360/0.025)*[realQuadOptic.DATA.T]);
title('phase')
legend('NOMNIAL','ACTUALQUADS')
subplot(2,2,4)
plot([nominalOptic.DATA.S], [nominalOptic.DATA.X]);
hold all;
plot([realQuadOptic.DATA.S], [realQuadOptic.DATA.X]);
title('X')
legend('NOMNIAL','ACTUALQUADS')


figure;
plot([nominalOptic.DATA.S], nomModX,'k','LineWidth',2);
hold all;
plot([realQuadOptic.DATA.S], rqModX,'b--','LineWidth',2);
errorbar(bpmPos,bpmHDiffFit(:,1),bpmHDiffFit_err(:,1),'ro','MarkerFaceColor','r');
legend('Nominal Optics','Real Quad Currents','BPM Data')
xlabel('S [m]')
ylabel('X [mm]')
title('TL2 Horizontal Position, 1V Amp Input')
format_plots;
plot([8.3646 8.3646],[-1.5 1],'k');
plot([19.7419 19.7419],[-1.5 1],'k');
% if savePlots; savePlot(saveDir,'orbClosureVsQuadModel'); end;

figure;
plot([nominalOptic.DATA.S], nomModX,'k','LineWidth',2);
hold all;
plot([ampRatOptic.DATA.S], ampRatModX,'g--','LineWidth',2);
errorbar(bpmPos,bpmHDiffFit(:,1),bpmHDiffFit_err(:,1),'ro','MarkerFaceColor','r');
legend('Nominal Optics','Real Amplifier Ratio','BPM Data')
xlabel('S [m]')
ylabel('X [mm]')
title('TL2 Horizontal Position, 1V Amp Input')
format_plots;
plot([8.3646 8.3646],[-1.5 1],'k');
plot([19.7419 19.7419],[-1.5 1],'k');
if savePlots; savePlot(saveDir,'orbClosureVsAmpModel'); end;
