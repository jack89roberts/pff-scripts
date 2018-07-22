dataSetName = '20150430_1731_CalibrationMon1Mix1FONT_Mon2Mix3_Mon3Mix2_2Mechanical';
sampleRange = 60:120;%40:70;%[]; % if this is left empty, will prompt to ask for sample range
frascatiCalTimeStamp= [];

%% load data
dataDir = ['data/' dataSetName];

CTFData = mergeMatMonData(dataDir);

if (isempty(frascatiCalTimeStamp))
    [frascatiCalTimeStamp,devNames,devValues] = loadInitSettings([dataDir '/initSettings.dat']);
end
[calibrationConstants, useMixerOverSqrtDiode] = loadFrascatiCalibrationConstants(sprintf('FrascatiCalibrations/data/frascatiCalibrationConstants_%s', frascatiCalTimeStamp));

% JACK - calibration constants for mechanical shifter on mon 2
calibrationConstants(2,1) = 0.325;
calibrationConstants(2,4) = 0.025;

%% signal processing

[ mixers, diodes ] = extractMixerDiode(CTFData);
[nMons,nPulses,nSamples] = size(mixers);

sampInterval = extractCTFSignalFromMergedData('CT_SCOPE01_CH01.Acquisition.sampleInterval.value',CTFData);
sampInterval = sampInterval(1);
timeAxis = sampInterval.*(0:(nSamples-1));

% align mixers, diodes
for mon=1:nMons
    Diode = squeeze(diodes(mon,:,:));
    Mixer = squeeze(mixers(mon,:,:));
    [Diode,Mixer] = getAligned2(Diode,20,nSamples,1,{Mixer});
    diodes(mon,:,:) = Diode;
    mixers(mon,:,:) = Mixer{1};
end
% prompt for sample range to use if not given
if (isempty(sampleRange))
    a = figure();
    
    for mon=1:nMons
        subplot(1,2,1)
        plot(squeeze(mixers(mon,1,:)));
        hold all;
        subplot(1,2,2)
        plot(squeeze(diodes(mon,1,:)));
        hold all;
    end
    
    subplot(1,2,1);
    title('MIXERS');
    subplot(1,2,2);
    title('DIODES')
    
    startSample = input('Start sample: ');
    endSample = input('End sample: ');
    sampleRange = startSample:endSample;
    
    close(a);
end

% remove duplicate pulses
for mon=1:nMons
    Diode = squeeze(diodes(mon,:,:));
    Mixer = squeeze(mixers(mon,:,:));

    diodes(mon,:,:) = removeDuplicatePulses(Diode);
    mixers(mon,:,:) = removeDuplicatePulses(Mixer);
end

% remove bad pulses
for mon=1:nMons
    Diode = squeeze(diodes(mon,:,:));
    Mixer = squeeze(mixers(mon,:,:));
    
    [Diode,Mixer] = removeBadPulses(Diode,sampleRange,{Mixer});
    diodes(mon,:,:) = Diode;
    mixers(mon,:,:) = Mixer{1};

end

% calculate phases
phases = NaN(nMons,nPulses,nSamples);
for mon=1:nMons
    CalibrationConstant = calibrationConstants(mon,1);
    CalibrationOffset = calibrationConstants(mon,4);
    Mixer =  squeeze(mixers(mon,:,:));
    Diode = squeeze(diodes(mon,:,:));

    if (useMixerOverSqrtDiode)
        phases(mon,:,:) = getPhaseMixerDiode( Mixer, Diode, CalibrationConstant, CalibrationOffset );
    else
        phases(mon,:,:) = getPhaseMixerDiode( Mixer, [], CalibrationConstant, CalibrationOffset );
    end
end

%% mean/std phases: calculate and plot
meanPhaseAlongPulse = squeeze(nanmean(phases,2));
stdPhaseAlongPulse = squeeze(nanstd(phases,0,2));

figure;
plot(timeAxis,meanPhaseAlongPulse(1,:));
hold all;
plot(timeAxis,meanPhaseAlongPulse(2,:));
plot(timeAxis,meanPhaseAlongPulse(3,:));
title('Mean Phase Along Pulse');
legend('Mix1','Mix2','Mix3');

figure;
plot(timeAxis,stdPhaseAlongPulse(1,:));
hold all;
plot(timeAxis,stdPhaseAlongPulse(2,:));
plot(timeAxis,stdPhaseAlongPulse(3,:));
title('Std Phase Along Pulse');
legend('Mix1','Mix2','Mix3');

meanPulsePhase = squeeze(nanmean(phases(:,:,sampleRange),3));
stdPulsePhase = squeeze(nanstd(phases(:,:,sampleRange),0,3));
stdMeanPulsePhase = squeeze(nanstd(meanPulsePhase,0,2))

figure;
plot(meanPulsePhase(1,:));
hold all;
plot(meanPulsePhase(2,:));
plot(meanPulsePhase(3,:));
title('Mean Pulse Phase')
legend('Mix1','Mix2','Mix3');

figure;
plot(stdPulsePhase(1,:));
hold all;
plot(stdPulsePhase(2,:));
plot(stdPulsePhase(3,:));
title('Std Pulse Phase')
legend('Mix1','Mix2','Mix3');


%% resolutions: calculate and plot
diffPhase12 = squeeze(phases(1,:,:)-phases(2,:,:));
diffPhase13 = squeeze(phases(1,:,:)-phases(3,:,:));
diffPhase23 = squeeze(phases(2,:,:)-phases(3,:,:));

resolution12 = squeeze(nanstd(diffPhase12,0,1))./sqrt(2);
resolution13 = squeeze(nanstd(diffPhase13,0,1))./sqrt(2);
resolution23 = squeeze(nanstd(diffPhase23,0,1))./sqrt(2);

figure;
plot(timeAxis,resolution12);
hold all;
plot(timeAxis,resolution13);
plot(timeAxis,resolution23);
legend('Mix1 Mix2','Mix1 Mix3','Mix2 Mix3');
title('RESOLUTION')

%% BPR
bpr475Cal = 53.845515;
bpr475Offset = 0.238898;
rawBPR475 = extractCTFSignalFromMergedData('CL_STBPR0475S.Samples.samples.value',CTFData);
phaseBPR475 = 4.*getPhaseMixerDiode(rawBPR475,[],bpr475Cal,bpr475Offset); % inc. 3GHz -> 12GHz
meanPulsePhaseBPR475 = nanmean(phaseBPR475);
meanPulsePhaseBPR475 = meanPulsePhaseBPR475 - nanmean(meanPulsePhaseBPR475(110:120));
stdPulsePhaseBPR475 = nanstd(phaseBPR475);
[~,nSamplesBPR475] = size(phaseBPR475);
timeAxBPR475 = (1000/96).*(0:(nSamplesBPR475-1));

figure;
subplot(1,2,1)
plot(meanPulsePhaseBPR475);
title('MEAN BPR CL 475')
subplot(1,2,2)
plot(stdPulsePhaseBPR475);
title('STD BPR CL 475')

% 
% bpr505Cal = 0.856447;
% bpr505Offset = -1.086647;
% rawBPR505 = extractCTFSignalFromMergedData('CR_STBPR0505S.Samples.samples.value',CTFData);
% phaseBPR505 = 4.*getPhaseMixerDiode(rawBPR505,[],bpr505Cal,bpr505Offset); % inc. 3GHz -> 12GHz
% meanPulsePhaseBPR505 = nanmean(phaseBPR505);
% meanPulsePhaseBPR505 = meanPulsePhaseBPR505 - nanmean(meanPulsePhaseBPR505(110:120));
% stdPulsePhaseBPR505 = nanstd(phaseBPR505);
% [~,nSamplesBPR505] = size(phaseBPR505);
% timeAxBPR505 = (1000/96).*(0:(nSamplesBPR505-1));
% 
% figure;
% subplot(1,2,1)
% plot(meanPulsePhaseBPR505);
% title('MEAN BPR CR 505')
% subplot(1,2,2)
% plot(stdPulsePhaseBPR505);
% title('STD BPR CR 505')
% 

bprCCCal = 3.493476;
bprCCOffset = 1.526760;
bprCCSampleRange = 390:410;
rawBPRCC = extractCTFSignalFromMergedData('CC_STBPR0915S.Samples.samples.value',CTFData);
phaseBPRCC = 4.*getPhaseMixerDiode(rawBPRCC,[],bprCCCal,bprCCOffset); % inc. 3GHz -> 12GHz
phaseBPRCC = -phaseBPRCC;
meanPulsePhaseBPRCC = nanmean(phaseBPRCC);
meanPulsePhaseBPRCC = meanPulsePhaseBPRCC; %- nanmean(meanPulsePhaseBPRCC(bprCCSampleRange));
stdPulsePhaseBPRCC = nanstd(phaseBPRCC);
[~,nSamplesBPRCC] = size(phaseBPRCC);
timeAxBPRCC = (1000/(96*2)).*(0:(nSamplesBPRCC-1));
meanVsTimeBPRCC = nanmean(phaseBPRCC(:,bprCCSampleRange),2);

figure;
subplot(1,2,1)
plot(meanPulsePhaseBPRCC);
title('MEAN BPR CC 915')
subplot(1,2,2)
plot(stdPulsePhaseBPRCC);
title('STD BPR CC 915')

figure;
scatter(meanVsTimeBPRCC,meanPulsePhase(3,:));

% figure;
% plot(timeAxis,stdPhaseAlongPulse(1,:));
% hold all;
% plot(timeAxis,stdPhaseAlongPulse(2,:));
% plot(timeAxis,stdPhaseAlongPulse(3,:));
% plot(timeAxBPR475-950,stdPulsePhaseBPR475);
% % plot(timeAxBPR505,stdPulsePhaseBPR505);
% title('Std Phase Along Pulse');
% legend('Mix1','Mix2','Mix3','CL.BPR0475');%,'CR.BPR505');

%%
petsI = double(extractCTFSignalFromMergedData('CE_SCOPE03_CH01.Acquisition.value.value',CTFData));
petsQ = double(extractCTFSignalFromMergedData('CE_SCOPE03_CH02.Acquisition.value.value',CTFData));

[petsPhase,petsPower] = getPhaseIQ(petsI, petsQ);

petsSampleRange = 190:320;
meanVsTimePets = nanmean(petsPhase(:,petsSampleRange),2);

figure;
plot(petsPhase')

figure;
scatter(meanVsTimePets,meanPulsePhase(3,:))
nancorrcoef(meanVsTimePets,meanPulsePhase(3,:)')

figure;
scatter(meanVsTimeBPRCC,meanVsTimePets)

%%
figure;
plot(meanPulsePhase(3,:)-nanmean(meanPulsePhase(3,:)))
hold all;
plot(meanVsTimePets-nanmean(meanVsTimePets));
plot(meanVsTimeBPRCC-nanmean(meanVsTimeBPRCC));
