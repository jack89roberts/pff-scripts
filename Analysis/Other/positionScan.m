clearvars; close all;
%% input parameters
%dataDir = '/home/jack/PhaseFeedforward/CTFData/201512/20151207_1258_PhMonPosScan_Horizontal';

dataDir = '/home/jack/PhaseFeedforward/CTFData/201603/20160316_1144_PosScanCT_H';
% dataDir = '/home/jack/PhaseFeedforward/CTFData/201603/20160316_1333_PosScanCT_H';

phaseSampleRange = 650:700;%[];%630:680;
frascatiCalTimeStamp = [];

% phaseDelayMon2 = 0;
% monPosRange = [];
% monPosDelayMon2 = [];
% bpmRange = [];
% bpmDelays = [];

monPositions = [0 0.4]; % positions of monitors relative to Mon1

bpmNames = {'CT_SVBPM0335','CT_SVBPM0430','CT_SVBPI0487','CT_SVBPI0495','CT_SVBPM0515'};
bpmPositions = [-2.3728 2.5241 7.6377 11.7109 16.6984]; % positions of BPMs relative to Mon1
bpmIndsPosCalc = [2 3]; % 2 bpms to use to back propagate position to monitors

corrNames = {'CT_DD0360','CT_DF0365'}; % without H/V
refCorrIndex = 1; % which corrector to use for splitting values etc.

%% load data (each file individiually to allow large datasets to be processed)
addpath('../../');

fileList = dir([dataDir '/*.mat']);
nPulses = length(fileList); 

nBPMs = length(bpmNames);
nSamplesBPM = 384;
sampFreqBPM = 10.417;

nCorrs = length(corrNames);

if (isempty(frascatiCalTimeStamp))
    [frascatiCalTimeStamp,~,~] = loadInitSettings([dataDir '/initSettings.dat']);
    calDir = ['/home/jack/PhaseFeedforward/CTFData/' frascatiCalTimeStamp(1:6) '/FrascatiCalibrations'];
    [calibrationConstants, useMixerOverSqrtDiode] = loadFrascatiCalibrationConstants(sprintf([calDir '/frascatiCalibrationConstants_' frascatiCalTimeStamp]));
end

fprintf('Loading %d files...\n',nPulses);
tic
arraysInit = 0;
for p=1:nPulses
    if (mod(p,round(nPulses/20))==0)
        fprintf('%d percent (%.2f seconds for %d files)...\n',round(100*p/nPulses),toc,p);
    end
    
    try
    load([dataDir '/' fileList(p).name]);
    
    % transmission check: skip if sum in 402 worse than -2A
    tmpTransmission = min(myDataStruct.CL_SVBPM0402S.Samples.samples.value);
%     if (transmission < -2)
    
        % extract signals from data struct
        [tmpMixers,tmpDiodes] = extractMixerDiode(myDataStruct);
        tmpPhases{1} = getPhaseMixerDiode(tmpMixers(1,:),[],calibrationConstants(1,1),calibrationConstants(1,4));
        tmpPhases{2} = getPhaseMixerDiode(tmpMixers(2,:),[],calibrationConstants(2,1),calibrationConstants(2,4));
        [tmpMonPosH,tmpMonPosV] = extractMonPos(myDataStruct);
        tmpBPMH = NaN(nBPMs,nSamplesBPM);
        tmpBPMV = NaN(nBPMs,nSamplesBPM);
        tmpBPMS = NaN(nBPMs,nSamplesBPM);
        for b=1:nBPMs
            [tmpBPMH(b,:), tmpBPMS(b,:), tmpBPMV(b,:)] = extractBPMFromCTFData(bpmNames{b},myDataStruct);
        end
        tmpCorrH = NaN(1,nCorrs);
        tmpCorrV = NaN(1,nCorrs);
        for c=1:nCorrs    
            [tmpCorrH(c),tmpCorrV(c)] = extractCorrFromCTFData(corrNames{c},myDataStruct);
        end
        
        % initialise arrays
        if (~arraysInit)
            [~,nSamplesPhase] = size(tmpMixers);
            [~,nSamplesMonPos] = size(tmpMonPosH);
%             mixers = NaN(2,nPulses,nSamplesPhase);
            phases = NaN(2,nPulses,nSamplesPhase);
            diodes = NaN(2,nPulses,nSamplesPhase);
            monPosH = NaN(2,nPulses,nSamplesMonPos);
            monPosV = NaN(2,nPulses,nSamplesMonPos);
            corrsH = NaN(nCorrs,nPulses);
            corrsV = NaN(nCorrs,nPulses);
            bpmH = NaN(nBPMs,nPulses,nSamplesBPM);
            bpmV = NaN(nBPMs,nPulses,nSamplesBPM);
            bpmS = NaN(nBPMs,nPulses,nSamplesBPM);
            transmission = NaN(1,nPulses);
            
            sampFreqPhase = myDataStruct.CT_SCOPE01_CH01.Acquisition.sampleInterval.value;
            sampFreqMonPos = myDataStruct.CT_SCOPE02_CH01.Acquisition.sampleInterval.value;
            
            arraysInit = 1;
        end
    
        % save data to arrays
%         mixers(:,p,:) = tmpMixers(1:2,:);
        diodes(:,p,:) = tmpDiodes(1:2,:);
        phases(1,p,:) = tmpPhases{1};
        phases(2,p,:) = tmpPhases{2};
        monPosH(:,p,:) = tmpMonPosH;
        monPosV(:,p,:) = tmpMonPosV;
        corrsH(:,p) = tmpCorrH;
        corrsV(:,p) = tmpCorrV;
        bpmH(:,p,:) = tmpBPMH;
        bpmV(:,p,:) = tmpBPMV;
        bpmS(:,p,:) = tmpBPMS;
        transmission(p) = tmpTransmission;
        
%     end
    catch
        fprintf('Error with file %d: %s\n',p,fileList(p).name);
    end
    
end
clear myDataStruct tmpMixers tmpDiodes tmpPhases tmpMonPosH tmpMonPosV tmpBPMH tmpBPMV tmpBPMS tmpCorrH tmpCorrV;
fprintf('%d percent (%.2f seconds for %d files).\n',round(100*p/nPulses),toc,p);

%% align signals
fprintf('Aligning signals...\n');

[diodes(1,:,:), others, phasePulseRange] = getAlignedXCorr(squeeze(diodes(1,:,:)),'end',{squeeze(phases(1,:,:))});
phases(1,:,:) = others{1};
if (isempty(phasePulseRange))
    [diodes(2,:,:),others, phasePulseRange] = getAlignedXCorr(squeeze(diodes(2,:,:)),'end',{squeeze(phases(2,:,:))});
else
    [diodes(2,:,:),others] = getAlignedXCorr(squeeze(diodes(2,:,:)),'end',{squeeze(phases(2,:,:))});    
end
phases(2,:,:) = others{1};

[monPosH(1,:,:), others,monPosPulseRange] = getAlignedXCorr(squeeze(monPosH(1,:,:)),'end',{squeeze(monPosV(1,:,:))});
monPosV(1,:,:) = others{1};
[monPosH(2,:,:), others] = getAlignedXCorr(squeeze(monPosH(2,:,:)),'end',{squeeze(monPosV(2,:,:))});
monPosV(2,:,:) = others{1};
   
for b=1:nBPMs
    [bpmS(b,:,:),others,tmpBPMPulseRange] = getAlignedXCorr(squeeze(bpmS(b,:,:)),'end',{squeeze(bpmH(b,:,:)),squeeze(bpmV(b,:,:))});
    bpmH(b,:,:) = others{1};
    bpmV(b,:,:) = others{2};
    
    if b==1
        bpmPulseRange = tmpBPMPulseRange;
    end  
end
clear tmpBPMPulseRange;

%% plot/ask for sample ranges/delays
fprintf('Selecting sample range...\n');

figure;
subplot(1,2,1)
plot(squeeze(diodes(1,1,:)));
hold all
plot(squeeze(diodes(2,1,:)));
legend('Mon1','Mon2')
xlabel('Sample No.')
ylabel('Output')
title('DIODE')
subplot(1,2,2)
plot(squeeze(phases(1,1,:)));
hold all
plot(squeeze(phases(2,1,:)));
legend('Mon1','Mon2')
xlabel('Sample No.')
ylabel('Output')
title('PHASE')
if (isempty(phaseSampleRange))
    phaseSampleRange = input('Sample Range: ');
end
subplot(1,2,1)
plot([phaseSampleRange(1) phaseSampleRange(1)],get(gca,'YLim'),'k');
plot([phaseSampleRange(end) phaseSampleRange(end)],get(gca,'YLim'),'k');
subplot(1,2,2)
plot([phaseSampleRange(1) phaseSampleRange(1)],get(gca,'YLim'),'k');
plot([phaseSampleRange(end) phaseSampleRange(end)],get(gca,'YLim'),'k');

startDelay = (phaseSampleRange(1)-phasePulseRange(1)).*sampFreqPhase;
endDelay = (phasePulseRange(end)-phaseSampleRange(end)).*sampFreqPhase;

monPosStart = monPosPulseRange(1)+round(startDelay./sampFreqMonPos);
monPosEnd = monPosPulseRange(end)-round(endDelay./sampFreqMonPos);
monPosSampleRange = monPosStart:monPosEnd;

bpmStart = bpmPulseRange(1)+round(startDelay./sampFreqBPM);
bpmEnd = bpmPulseRange(end)-round(endDelay./sampFreqBPM);
bpmSampleRange = bpmStart:bpmEnd;

figure;
subplot(1,2,1)
plot(squeeze(monPosH(1,1,:)))
hold all;
plot(squeeze(monPosH(2,1,:)))
plot([monPosSampleRange(1) monPosSampleRange(1)],get(gca,'YLim'),'k');
plot([monPosSampleRange(end) monPosSampleRange(end)],get(gca,'YLim'),'k');
legend('Mon1','Mon2')
xlabel('Sample No.')
ylabel('Output')
title('MON POS H')
subplot(1,2,2)
plot(squeeze(monPosV(1,1,:)))
hold all;
plot(squeeze(monPosV(2,1,:)))
plot([monPosSampleRange(1) monPosSampleRange(1)],get(gca,'YLim'),'k');
plot([monPosSampleRange(end) monPosSampleRange(end)],get(gca,'YLim'),'k');
legend('Mon1','Mon2')
xlabel('Sample No.')
ylabel('Output')
title('MON POS V')

figure;
for b=1:nBPMs
    plot(squeeze(bpmS(b,1,:)))
    hold all;
end
legend(bpmNames)
xlabel('Sample No.')
ylabel('Output')
title('BPM S')
plot([bpmSampleRange(1) bpmSampleRange(1)],get(gca,'YLim'),'k');
plot([bpmSampleRange(end) bpmSampleRange(end)],get(gca,'YLim'),'k');

%% split arrays in to different corrector settings
fprintf('Splitting in to corrector settings...\n');

% [sortedCorr,sortedOrder] = sort(corrsH(1,:));
% sortedCorr = sortedCorr(~isnan(sortedCorr));
% sortedOrder = sortedOrder(~isnan(sortedCorr));
% corrChangedInd = find(diff(sortedCorr)>0.05);

% find points where corrector value changes by more than 0.05 A
corrChangedIndH = find(abs(diff(corrsH(refCorrIndex,:)))>0.05);
corrChangedIndV = find(abs(diff(corrsV(refCorrIndex,:)))>0.05);

% determine whether the corrector is being changed horizontally or
% vertically
isHorizontal = length(corrChangedIndH)>length(corrChangedIndV);
if (isHorizontal)
    corrChangedInd = corrChangedIndH;
else
    corrChangedInd = corrChangedIndV;
end

nCorrSettings = length(corrChangedInd)+1;

splitPhases = cell(2,nCorrSettings,nSamplesPhase);
splitMonPosH = cell(1,nCorrSettings);
splitMonPosV = cell(1,nCorrSettings);
splitCorrsH = cell(1,nCorrSettings);
splitCorrsV = cell(1,nCorrSettings);
splitBPMH = cell(1,nCorrSettings);
splitBPMV = cell(1,nCorrSettings);
splitBPMS = cell(1,nCorrSettings);
splitTransmission = cell(1,nCorrSettings);
splitNPulses = NaN(1,nCorrSettings);

for s=1:nCorrSettings
    if s==1
        startInd = 1;
    else
        startInd = corrChangedInd(s-1)+1;
    end
    
    if s==nCorrSettings
        endInd = nPulses;
    else
        endInd = corrChangedInd(s);
    end
    
    pulseRange = startInd:endInd;
    
    splitPhases{s} = phases(:,pulseRange,:);
    splitMonPosH{s} = monPosH(:,pulseRange,:);
    splitMonPosV{s} = monPosV(:,pulseRange,:);
    splitCorrsH{s} = corrsH(:,pulseRange);
    splitCorrsV{s} = corrsV(:,pulseRange);
    splitBPMH{s} = bpmH(:,pulseRange,:);
    splitBPMV{s} = bpmV(:,pulseRange,:);
    splitBPMS{s} = bpmS(:,pulseRange,:);
    splitTransmission{s} = transmission(pulseRange);
    splitNPulses(s) = length(pulseRange);
end

goodCorrSettings = splitNPulses > 10; % For plots, fits etc. only include split sections with more than 10 pulses

%% remove pulses no transmission and bad pulses at each corrector setting
fprintf('Removing bad pulses...\n');

for s=1:nCorrSettings
    
    [~,N,~] = size(splitPhases{s});
    
    % remove pulses no beam
    for p=1:N        
        if splitTransmission{s}(p)>-2
            splitPhases{s}(:,p,:) = NaN;
            splitMonPosH{s}(:,p,:) = NaN;    
            splitMonPosV{s}(:,p,:) = NaN;    
            splitCorrsH{s}(p) = NaN;    
            splitCorrsV{s}(p) = NaN;    
            splitBPMH{s}(:,p,:) = NaN;    
            splitBPMV{s}(:,p,:) = NaN;    
            splitBPMS{s}(:,p,:) = NaN;    
            splitTransmission{s}(p) = NaN;    
        end
    end
    
    % remove outliers
    splitPhases{s}(1,:,:) = removeBadPulses(squeeze(splitPhases{s}(1,:,:)),phaseSampleRange);
    splitPhases{s}(2,:,:) = removeBadPulses(squeeze(splitPhases{s}(2,:,:)),phaseSampleRange);
    splitMonPosH{s}(1,:,:) = removeBadPulses(squeeze(splitMonPosH{s}(1,:,:)),monPosSampleRange);
    splitMonPosH{s}(2,:,:) = removeBadPulses(squeeze(splitMonPosH{s}(2,:,:)),monPosSampleRange);
    splitMonPosV{s}(1,:,:) = removeBadPulses(squeeze(splitMonPosV{s}(1,:,:)),monPosSampleRange);
    splitMonPosV{s}(2,:,:) = removeBadPulses(squeeze(splitMonPosV{s}(2,:,:)),monPosSampleRange);    
    for bpm=1:nBPMs
        % for bad transmission in BPM, remove H and V as well
        [splitBPMS{s}(b,:,:),others] = removeBadPulses(squeeze(splitBPMS{s}(b,:,:)),bpmSampleRange,{squeeze(splitBPMH{s}(b,:,:)),squeeze(splitBPMV{s}(b,:,:))});
        splitBPMH{s}(b,:,:) = others{1};
        splitBPMV{s}(b,:,:) = others{2};
        splitBPMH{s}(b,:,:) = removeBadPulses(squeeze(splitBPMH{s}(b,:,:)),bpmSampleRange);
        splitBPMV{s}(b,:,:) = removeBadPulses(squeeze(splitBPMV{s}(b,:,:)),bpmSampleRange);
    end
    for c=1:nCorrs
        splitCorrsH{s}(c,:) = removeBadPulses(splitCorrsH{s}(c,:));    
        splitCorrsV{s}(c,:) = removeBadPulses(splitCorrsV{s}(c,:));   
    end

end

%% calculate means and errors at each corr setting
fprintf('Calculating means...\n');

meanPhasesAlong = NaN(2,nCorrSettings,nSamplesPhase);
meanMonPosHAlong = NaN(2,nCorrSettings,nSamplesMonPos);
meanMonPosVAlong = NaN(2,nCorrSettings,nSamplesMonPos);
meanBPMHAlong = NaN(nBPMs,nCorrSettings,nSamplesBPM);
meanBPMVAlong = NaN(nBPMs,nCorrSettings,nSamplesBPM);
meanBPMSAlong = NaN(nBPMs,nCorrSettings,nSamplesBPM);
meanCorrsH = NaN(nCorrs,nCorrSettings);
meanCorrsV = NaN(nCorrs,nCorrSettings);
meanPhasesAlong_err = NaN(2,nCorrSettings,nSamplesPhase);
meanMonPosHAlong_err = NaN(2,nCorrSettings,nSamplesMonPos);
meanMonPosVAlong_err = NaN(2,nCorrSettings,nSamplesMonPos);
meanBPMHAlong_err = NaN(nBPMs,nCorrSettings,nSamplesBPM);
meanBPMVAlong_err = NaN(nBPMs,nCorrSettings,nSamplesBPM);
meanBPMSAlong_err = NaN(nBPMs,nCorrSettings,nSamplesBPM);
meanCorrsH_err = NaN(nCorrs,nCorrSettings);
meanCorrsV_err = NaN(nCorrs,nCorrSettings);

for s=1:nCorrSettings
    [meanPhasesAlong(:,s,:),~,meanPhasesAlong_err(:,s,:),~] = nanMeanStdErr(splitPhases{s},2);
    [meanMonPosHAlong(:,s,:),~,meanMonPosHAlong_err(:,s,:),~] = nanMeanStdErr(splitMonPosH{s},2);
    [meanMonPosVAlong(:,s,:),~,meanMonPosVAlong_err(:,s,:),~] = nanMeanStdErr(splitMonPosV{s},2);
    [meanBPMHAlong(:,s,:),~,meanBPMHAlong_err(:,s,:),~] = nanMeanStdErr(splitBPMH{s},2);
    [meanBPMVAlong(:,s,:),~,meanBPMVAlong_err(:,s,:),~] = nanMeanStdErr(splitBPMV{s},2);
    [meanBPMSAlong(:,s,:),~,meanBPMSAlong_err(:,s,:),~] = nanMeanStdErr(splitBPMS{s},2);
    [meanCorrsH(:,s),~,meanCorrsH_err(:,s),~] = nanMeanStdErr(splitCorrsH{s},2);
    [meanCorrsV(:,s),~,meanCorrsV_err(:,s),~] = nanMeanStdErr(splitCorrsV{s},2);
end

[meanPhases,~,meanPhases_err,~] = nanMeanStdErr(meanPhasesAlong(:,:,phaseSampleRange),3);
[meanMonPosH,~,meanMonPosH_err,~] = nanMeanStdErr(meanMonPosHAlong(:,:,monPosSampleRange),3);
[meanMonPosV,~,meanMonPosV_err,~] = nanMeanStdErr(meanMonPosVAlong(:,:,monPosSampleRange),3);
[meanBPMH,~,meanBPMH_err,~] = nanMeanStdErr(meanBPMHAlong(:,:,bpmSampleRange),3);
[meanBPMV,~,meanBPMV_err,~] = nanMeanStdErr(meanBPMVAlong(:,:,bpmSampleRange),3);
[meanBPMS,~,meanBPMS_err,~] = nanMeanStdErr(meanBPMSAlong(:,:,bpmSampleRange),3);

%% calculate corrector offsets, then subtract init orbit in each BPM

maxCorrsH = max(meanCorrsH,[],2);
minCorrsH = min(meanCorrsH,[],2);
maxCorrsV = max(meanCorrsV,[],2);
minCorrsV = min(meanCorrsV,[],2);

% use difference between max/min corrector setting to calculate initial
% values (assumes scan is initValue +/- maxOffset)
offsetCorrsH = (abs(maxCorrsH)-abs(minCorrsH))/2;
if (abs(maxCorrsH)>abs(minCorrsH))
    offsetCorrsH = abs(offsetCorrsH);
else
    offsetCorrsH = -abs(offsetCorrsH);
end

offsetCorrsV = (abs(maxCorrsV)-abs(minCorrsV))/2;
if (maxCorrsV>minCorrsV)
    offsetCorrsV = abs(offsetCorrsV);
else
    offsetCorrsV = -abs(offsetCorrsV);
end

nominalDiff = NaN(nCorrs,nCorrSettings);
for c=1:nCorrs
    if (isHorizontal)
        nominalDiff(c,:) = meanCorrsH(c,:)-offsetCorrsH(c);
    else
        nominalDiff(c,:) = meanCorrsV(c,:)-offsetCorrsV(c);    
    end
end
[~,refInd] = min(sum(abs(nominalDiff),1)); % find index to use as reference orbit (min difference to calculated initial corrector offsets)

% subtract init bpm offsets
for b=1:nBPMs
    meanBPMH(b,:) = meanBPMH(b,:)-meanBPMH(b,refInd);
    meanBPMV(b,:) = meanBPMV(b,:)-meanBPMV(b,refInd);
    meanBPMHAlong(b,:,:) = meanBPMHAlong(b,:,:)-meanBPMH(b,refInd);
    meanBPMVAlong(b,:,:) = meanBPMVAlong(b,:,:)-meanBPMV(b,refInd);
end

%% calculate position on monitors by back propagating from bpms
% assuming ballistic beam

tmp = (meanBPMH(bpmIndsPosCalc(2),:)-meanBPMH(bpmIndsPosCalc(1),:))./(bpmPositions(bpmIndsPosCalc(2))-bpmPositions(bpmIndsPosCalc(1)));
tmp = (monPositions(2)-bpmPositions(bpmIndsPosCalc(1))).*tmp;
tmp = tmp + meanBPMH(bpmIndsPosCalc(1),:);

tmp = (meanBPMH(bpmIndsPosCalc(2),:)-meanBPMH(bpmIndsPosCalc(1),:))./(bpmPositions(bpmIndsPosCalc(2))-bpmPositions(bpmIndsPosCalc(1)));
tmp = (monPositions(1)-bpmPositions(bpmIndsPosCalc(1))).*tmp;
tmp = tmp + meanBPMH(bpmIndsPosCalc(1),:);

%% plots
fprintf('Plotting...\n');

% for b=1:nBPMs
%     figure;
%     errorbar(meanBPMH(b,goodCorrSettings),meanPhases(1,goodCorrSettings),meanPhases_err(1,goodCorrSettings),'b.')
%     hold all
%     errorbar(meanBPMH(b,goodCorrSettings),meanPhases(2,goodCorrSettings),meanPhases_err(2,goodCorrSettings),'r.')
%     herrorbar(meanBPMH(b,goodCorrSettings),meanPhases(1,goodCorrSettings),meanBPMH_err(b,goodCorrSettings),'b.')
%     herrorbar(meanBPMH(b,goodCorrSettings),meanPhases(2,goodCorrSettings),meanBPMH_err(b,goodCorrSettings),'r.')
%     xlabel('BPMH POS [mm]')
%     ylabel('PHASE [degrees]')
%     title(bpmNames{b})
%     legend('Mon1','Mon2')
% end
% 
% for b=1:nBPMs
%     figure;
%     errorbar(meanBPMV(b,goodCorrSettings),meanPhases(1,goodCorrSettings),meanPhases_err(1,goodCorrSettings),'b.')
%     hold all
%     errorbar(meanBPMV(b,goodCorrSettings),meanPhases(2,goodCorrSettings),meanPhases_err(2,goodCorrSettings),'r.')
%     herrorbar(meanBPMV(b,goodCorrSettings),meanPhases(1,goodCorrSettings),meanBPMV_err(b,goodCorrSettings),'b.')
%     herrorbar(meanBPMV(b,goodCorrSettings),meanPhases(2,goodCorrSettings),meanBPMV_err(b,goodCorrSettings),'r.')
%     xlabel('BPMV POS [mm]')
%     ylabel('PHASE [degrees]')
%     title(bpmNames{b})
%     legend('Mon1','Mon2')
% end
% 
for b=1:nBPMs
    figure;
    errorbar(meanBPMH(b,goodCorrSettings),meanMonPosH(1,goodCorrSettings),meanMonPosH_err(1,goodCorrSettings),'b.')
    hold all
    errorbar(meanBPMH(b,goodCorrSettings),meanMonPosH(2,goodCorrSettings),meanMonPosH_err(2,goodCorrSettings),'r.')
    herrorbar(meanBPMH(b,goodCorrSettings),meanMonPosH(1,goodCorrSettings),meanBPMH_err(b,goodCorrSettings),'b.')
    herrorbar(meanBPMH(b,goodCorrSettings),meanMonPosH(2,goodCorrSettings),meanBPMH_err(b,goodCorrSettings),'r.')
    xlabel('BPMH POS [mm]')
    ylabel('MON H POS [V]')
    title(bpmNames{b})
    legend('Mon1','Mon2')
end

% for b=1:nBPMs
%     figure;
%     errorbar(meanBPMV(b,goodCorrSettings),meanMonPosV(1,goodCorrSettings),meanMonPosV_err(1,goodCorrSettings),'b.')
%     hold all
%     errorbar(meanBPMV(b,goodCorrSettings),meanMonPosV(2,goodCorrSettings),meanMonPosV_err(2,goodCorrSettings),'r.')
%     herrorbar(meanBPMV(b,goodCorrSettings),meanMonPosV(1,goodCorrSettings),meanBPMV_err(b,goodCorrSettings),'b.')
%     herrorbar(meanBPMV(b,goodCorrSettings),meanMonPosV(2,goodCorrSettings),meanBPMV_err(b,goodCorrSettings),'r.')
%     xlabel('BPMV POS [mm]')
%     ylabel('MON V POS [V]')
%     title(bpmNames{b})
%     legend('Mon1','Mon2')
% end

