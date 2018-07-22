%% BPR
return;
%% signal processing

[bprCals,bprCalNames]=loadBPRCalibrationConstants('/home/jack/PhaseFeedforward/CTFData/201511/BPRCalibrations/bprCalibrationConstants_20151123_1644');

ctCal = bprCals(3,:);
ctBPRRaw = extractCTFSignalFromMergedData('CT_STBPR0532S.Samples.samples.value',CTFData);

% ctCal = bprCals(4,:);
% ctBPRRaw = extractCTFSignalFromMergedData('CR_STBPR0505S.Samples.samples.value',CTFData);


ctBPRPhase = 4.*real(asind((ctBPRRaw-ctCal(4))./ctCal(1)));

figure;
plot(ctBPRPhase')

nanstd(ctBPRPhase(:,250))

meanBPRAlongPulse = nanmean(ctBPRPhase,1);
meanBPRAlongPulse = meanBPRAlongPulse-meanBPRAlongPulse(240);

figure;
plot(meanBPRAlongPulse)

figure;
plot(meanPhaseAlongPulse(2,:)-meanPhaseAlongPulse(2,654))

bprTimeAxis = 10.4*((1:384)-156);

figure;
plot(timeAxis,meanPhaseAlongPulse(1,:)-meanPhaseAlongPulse(1,654));
hold all;
plot(timeAxis,0.93*(delaySignal(meanPhaseAlongPulse(2,:),-5)-meanPhaseAlongPulse(2,654)));
plot(bprTimeAxis,meanBPRAlongPulse);