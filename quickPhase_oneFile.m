%calibConstFileName = '/home/jack/Documents/CTFData/PFF/201411/frascatiCalibrations/frascatiCalibration_20141113_1455_calResults_mean.mat'; % file containing calibration results


dataDirectory = 'data/20141215_1303_ElectronicsSwapped/'; % directory containing files saved from matlab monitor


%monitorIndices = [1, 2, 3]; % monitors included in calibration/data files
monitorNames = {'Mon1','Mon2','Mon3'};

calFile = 'frascatiCalibrationConstants';%'data/frascatiCalibrations/frascatiCalibrationConstants_20141210_1844'; %'frascatiCalibrationConstants'; % 
allCalVals = loadFrascatiCalibrationConstants(calFile);
calibrationConstants = allCalVals(:,1);
calibrationOffsets = allCalVals(:,4);
calibrationConstants = [1.5996, 0.0648, 0.00121];
calibrationOffsets = [0.3513, 0.1036, -0.0009];

useMixerOverSqrtDiode = false;


d = 80;
startSamples = [160-d, 160-d, 320-d];
endSamples = [230-d, 230-d, 390-d];
nMonitors = 3;

petsStartSample = 300;
petsEndSample = 700;


removeBadPulses = false;%true; % if true remove bad pulses
removeThreshold = 0.01; % remove pulses that are removeThreshold % away from max on diode
removeMonIndices = [1, 2]; % look at diode for these monitors when deciding which pulses to remove (require good transmission for all monitors with given indices)

%% Load data, remove bad pulses

% calFile = open(calibConstFileName);
% calibrationConstants = calFile.meanCalibrationFactors;
% calStartSample = calFile.calStartSample;
% calEndSample = calFile.calEndSample;

CTFData = mergeMatMonData(dataDirectory);
[mixers,diodes] = extractMixerDiode(CTFData);
[~,nPulses,nSamples] = size(mixers);

plot(squeeze(mixers(1,1,:)))


if (removeBadPulses==true && nPulses>=20) % require at least 20 pulses just so we actually have something to take a decent mean from
    transmission = cell(1,nMonitors);
    isGoodPulseMon = cell(1,nMonitors);
    for mon = 1:nMonitors
        transmission{mon} = squeeze(nanmean(diodes(mon,:,startSamples(mon):endSamples(mon)),3));
        isGoodPulseMon{mon} = (transmission{mon}-min(transmission{mon})) < removeThreshold;
    end
    
    isGoodPulseAll = ones(1,nPulses);
    for mon = removeMonIndices
        isGoodPulseAll = isGoodPulseAll & isGoodPulseMon{mon};
    end
    
    CTFData = CTFData(isGoodPulseAll);
    mixers = mixers(:,isGoodPulseAll,:);
    diodes = diodes(:,isGoodPulseAll,:);
    nPulses = length(CTFData);
    
end

%% Process Frascati data

phases = NaN*ones(nMonitors,nPulses,nSamples);
for mon = 1:nMonitors
    if (useMixerOverSqrtDiode)
        phases(mon,:,:) = getPhaseMixerDiode(mixers(mon,:,:), diodes(mon,:,:), calibrationConstants(mon), calibrationOffsets(mon)); 
    else
        phases(mon,:,:) = getPhaseMixerDiode(mixers(mon,:,:), [], calibrationConstants(mon), calibrationOffsets(mon)); 
    end
end

meanDiodes = squeeze(nanmean(diodes,2));
meanMixers = squeeze(nanmean(mixers,2));
meanPhases = squeeze(nanmean(phases,2)); % mean per sample
stdPhases = squeeze(nanstd(phases,0,2)); % std per sample
meanStdPhase = NaN*ones(1,nMonitors);
for mon=1:nMonitors
    meanStdPhase(mon) = squeeze(nanmean(stdPhases(mon,startSamples(mon):endSamples(mon)),2)); % mean sample std
end

nAvg = 10;
meanStdPhase_nAvgPulses = NaN*ones(nMonitors,ceil(nPulses/nAvg));
for p=1:10:nPulses
    if p+nAvg-1 > nPulses
        pulseRange = p:nPulses;
    else
        pulseRange = p:p+nAvg-1;
    end
    
    for mon=1:nMonitors
        tmpStdPhase = squeeze(nanstd(phases(mon,pulseRange,:),0,2));
        meanStdPhase_nAvgPulses(mon,ceil(p/nAvg)) = nanmean(tmpStdPhase(startSamples(mon):endSamples(mon)));
    end
end

diffPhases12 = phases(1,:,startSamples(1):endSamples(1)) - phases(2,:,startSamples(2):endSamples(2));
diffPhases13 = phases(1,:,startSamples(1):endSamples(1)) - phases(3,:,startSamples(3):endSamples(3));
diffPhases23 = phases(2,:,startSamples(2):endSamples(2)) - phases(3,:,startSamples(3):endSamples(3));

resolution12 = squeeze(nanstd(diffPhases12,0,3))./sqrt(2);
resolution13 = squeeze(nanstd(diffPhases13,0,3))./sqrt(2);
resolution23 = squeeze(nanstd(diffPhases23,0,3))./sqrt(2);

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

corr12_meanPulsePhase = corrcoef(meanPulsePhases(1,:),meanPulsePhases(2,:));
corr13_meanPulsePhase = corrcoef(meanPulsePhases(1,:),meanPulsePhases(3,:));
corr23_meanPulsePhase = corrcoef(meanPulsePhases(2,:),meanPulsePhases(3,:));
corr12_meanPulsePhase = corr12_meanPulsePhase(1,2);
corr13_meanPulsePhase = corr13_meanPulsePhase(1,2);
corr23_meanPulsePhase = corr23_meanPulsePhase(1,2);

meanPhases_offsetSubtract = NaN*ones(nMonitors,nSamples);
for mon=1:nMonitors
    meanPhases_offsetSubtract(mon,:) = meanPhases(mon,:)-nanmean(meanPhases(mon,startSamples(mon):endSamples(mon)));
end
meanDiffPhase12 = meanPhases_offsetSubtract(1,startSamples(1):endSamples(1)) - meanPhases_offsetSubtract(2,startSamples(2):endSamples(2));
meanDiffPhase13 = meanPhases_offsetSubtract(1,startSamples(1):endSamples(1)) - meanPhases_offsetSubtract(3,startSamples(3):endSamples(3));
meanDiffPhase23 = meanPhases_offsetSubtract(2,startSamples(2):endSamples(2)) - meanPhases_offsetSubtract(3,startSamples(3):endSamples(3));

%% Process PETS data

pets_I = extractCTFSignalFromMergedData('CE.SCOPE03.CH01.Acquisition.value.value',CTFData);
pets_Q = extractCTFSignalFromMergedData('CE.SCOPE03.CH02.Acquisition.value.value',CTFData);
[pets_Phase,pets_Power] = getPhaseIQ(pets_I, pets_Q);

figure;
plot(pets_Phase(1,:));

pets_meanPulsePhases = nanmean(pets_Phase(:,petsStartSample:petsEndSample),2); % mean per pulse
pets_meanPhases = nanmean(pets_Phase); % mean per sample
pets_stdPhases = nanstd(pets_Phase); % std per sample
pets_meanStdPhase = nanmean(pets_stdPhases(petsStartSample:petsEndSample)); % mean sample std

corr_PetsFrascati3_meanPulsePhase = corrcoef(meanPulsePhases(3,:),pets_meanPulsePhases);
corr_PetsFrascati3_meanPulsePhase = corr_PetsFrascati3_meanPulsePhase(1,2);

%% Plots

figure;
plot(resolution12);
hold all;
plot(resolution13);
plot(resolution23);
title('Resolution vs. Sample No.');
legend('Mon1-Mon2', 'Mon1-Mon3', 'Mon2-Mon3');
xlabel('Sample No,');
ylabel('Resolution [12GHz Degrees]');

% figure;
% for i=1:nPulses
%     plot(squeeze(phases(3,i,:)),'b');
%     hold all
%     %plot(squeeze(phases(2,i,:)),'r');
%     %plot(squeeze(phases(3,i,:)),'g');
%     %hold off
%     %pause;
% end

figure;
for mon=1:nMonitors
    plot(stdPhases(mon,startSamples(mon):endSamples(mon)))
    hold all
end
plot(pets_stdPhases(petsStartSample:petsEndSample));
title('Standard Deviation Phase vs. Sample No.');
legend(monitorNames);
xlabel('Sample No,');
ylabel('Phase Jitter [12GHz Degrees]');
legStr = cell(1,4);
for mon=1:nMonitors
    legStr{mon} = sprintf('%s: %.2f degrees',monitorNames{mon},meanStdPhase(mon));
end
legStr{4} = sprintf('PETS: %.2f degrees',pets_meanStdPhase);
legend(legStr);

figure;
plot(meanPulsePhases(1,:),meanPulsePhases(2,:),'o');
xlabel('Mon1 Phase [12GHz Degrees]');
ylabel('Mon2 Phase [12GHz Degrees]');
title(sprintf('Correlation Mon1 - Mon2: %.2f',corr12_meanPulsePhase));

figure;
plot(meanPulsePhases(1,:),meanPulsePhases(3,:),'o');
hold all
xlabel('Mon1 Phase [12GHz Degrees]');
ylabel('Mon3 Phase [12GHz Degrees]');
title(sprintf('Correlation Mon1 - Mon3: %.2f',corr13_meanPulsePhase));


figure;
plot(meanPulsePhases(3,:),pets_meanPulsePhases,'o');
xlabel('Mon3 Phase [12GHz Degrees]');
ylabel('PETS Phase [12GHz Degrees]');
title(sprintf('Correlation Mon3 - PETS: %.2f',corr_PetsFrascati3_meanPulsePhase));

figure;
plot(meanPulsePhases(1,:));
hold all;
plot(meanPulsePhases(2,:))
plot(meanPulsePhases(3,:))
plot(pets_meanPulsePhases);
title('Mean Phase vs. Time');
legend([monitorNames 'PETS']);
xlabel('Time [Pulse No.]');
ylabel('Phase [12GHz Degrees]');

% figure;
% plot(meanDiffPhase12);
% hold all;
% plot(meanDiffPhase13);
% plot(meanDiffPhase23);
% title('Mean Diff Phase');
% legend('1-2','1-3','2-3');
% xlabel('Time [Pulse No.]');
% ylabel('Phase [12GHz Degrees]');

% figure;
% for mon = 1:nMonitors
%     plot(squeeze(meanDiodes(mon,:)));
%     hold all;
% end
% title('Mean Diode');
% xlabel('Sample No.');
% ylabel('Output [V]');
% legend(monitorNames);

% figure;
% for mon = 1:nMonitors
%     plot(squeeze(meanMixers(mon,:)));
%     hold all;
% end
% title('Mean Mixer');
% xlabel('Sample No.');
% ylabel('Output [V]');
% legend(monitorNames);

figure;
for mon = 1:nMonitors
    plot(squeeze(meanPhases_offsetSubtract(mon,:)));
    hold all;
end
plot(pets_meanPhases);
title('Mean Phase vs. Sample No.');
xlabel('Sample No.');
ylabel('Phase [12 GHz Degrees]');
legend([monitorNames 'PETS']);

figure;
plot(squeeze(phases(1,nPulses,:)));
hold all
plot(squeeze(phases(2,nPulses,:)));
plot(squeeze(phases(3,nPulses,:)));
plot(pets_Phase(nPulses,:));
title('Last Phase vs. Sample No.');
xlabel('Sample No.');
ylabel('Phase [12 GHz Degrees]');
legend([monitorNames 'PETS']);

%% BPMs

CT.BPM0285H = extractCTFSignalFromMergedData('CT.SVBPM0285H.MeanAtCursor.mean.value',CTFData);
CT.BPI0608H = extractCTFSignalFromMergedData('CT.SVBPI0608H.MeanAtCursor.mean.value',CTFData);
CC.BPI0645H = extractCTFSignalFromMergedData('CC.SVBPI0645H.MeanAtCursor.mean.value',CTFData);

figure;
subplot(2,1,1)
plot(CT.BPM0285H);
hold all;
plot(CT.BPI0608H);
plot(CC.BPI0645H);
legend('CT.BPM0285H','CT.BPI0608H','CC.BPI0645H');

subplot(2,1,2)
plot(squeeze(meanPulsePhases(1,:)))
hold all
plot(squeeze(meanPulsePhases(2,:)))
plot(squeeze(meanPulsePhases(3,:)))
plot(pets_meanPulsePhases)
legend([monitorNames 'PETS'])

figure;
plot(squeeze(meanPulsePhases(3,:)),CT.BPM0285H,'o');
hold all
plot(squeeze(meanPulsePhases(3,:)),CT.BPI0608H,'o');
plot(squeeze(meanPulsePhases(3,:)),CC.BPI0645H,'o');
legend('CT.BPM0285H','CT.BPI0608H','CC.BPI0645H');
title('Mon3 Correlation With BPMs');

figure;
plot(squeeze(meanPulsePhases(3,:)),CT.BPM0285H,'o');
hold all
plot(squeeze(meanPulsePhases(3,:)),CT.BPI0608H,'o');
plot(squeeze(meanPulsePhases(3,:)),CC.BPI0645H,'o');
legend('CT.BPM0285H','CT.BPI0608H','CC.BPI0645H');
title('Mon3 Correlation With BPMs');