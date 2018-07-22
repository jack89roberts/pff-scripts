dataFile = '/home/jack/PhaseFeedforward/CTFData/201412/allforstability/20141215_1411/processed/20141215_1411.mat';
load(dataFile)
bpmH = processedData.bpmData.CT.SVBPM0285H.Samples.samples;
bpmS = processedData.bpmData.CT.SVBPM0285S.Samples.samples;
bpmH = removeDuplicatePulses(bpmH);
bpmS = removeDuplicatePulses(bpmS);
[bpmS,bpmH] = getAligned2(bpmS,20,300,1,{bpmH});
bpmH = bpmH{1};
[bpmS,bpmH] = removeBadPulses(bpmS,63:75,{bpmH});
bpmH = bpmH{1};
stdBPMH = nanstd(bpmH);
figure;
plot(stdBPMH);
xlabel('Sample No.')
ylabel('Position Jitter [mm]')

bpmH = processedData.bpmData.CT.SVBPI0608H.Samples.samples;
bpmS = processedData.bpmData.CT.SVBPI0608S.Samples.samples;
bpmH = removeDuplicatePulses(bpmH);
bpmS = removeDuplicatePulses(bpmS);
[bpmS,bpmH] = getAligned2(bpmS,20,300,1,{bpmH});
bpmH = bpmH{1};
[bpmS,bpmH] = removeBadPulses(bpmS,112:131,{bpmH});
bpmH = bpmH{1};
stdBPMH = nanstd(bpmH);
figure;
plot(stdBPMH);
xlabel('Sample No.')
ylabel('Position Jitter [mm]')

bpmH = processedData.bpmData.CC.SVBPI0535H.Samples.samples;
bpmS = processedData.bpmData.CC.SVBPI0535S.Samples.samples;
bpmH = removeDuplicatePulses(bpmH);
bpmS = removeDuplicatePulses(bpmS);
[bpmS,bpmH] = getAligned2(bpmS,20,300,1,{bpmH});
bpmH = bpmH{1};
[bpmS,bpmH] = removeBadPulses(bpmS,112:131,{bpmH});
bpmH = bpmH{1};
stdBPMH = nanstd(bpmH);
figure;
plot(stdBPMH);
xlabel('Sample No.')
ylabel('Position Jitter [mm]')

