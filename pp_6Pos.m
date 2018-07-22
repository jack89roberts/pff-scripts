clf;
close all;

bpmNameS = 'CC.SVBPI0685S.Samples.samples.value';
bpmNameH = 'CC.SVBPI0685H.Samples.samples.value';

bpmTraceS = cell(1,nDataSets);
bpmTraceH = cell(1,nDataSets);

for i=1:nDataSets

    allProcessedPhaseData{i} = extractProcessedPhaseData(mergedData{i},calibrationConstants);
    
    bpmTraceS{i} = extractCTFSignalFromMergedData(bpmNameS,mergedData{i});
    bpmTraceH{i} = extractCTFSignalFromMergedData(bpmNameH,mergedData{i});
    
%    allMeanPhases{i} = allProcessedPhaseData{i}.meanPhases;
%    allMeanPulsePhases{i} = allProcessedPhaseData{i}.meanPulsePhases;
%    allStdPhases{i} = allProcessedPhaseData{i}.stdPhases;
    
%    allMeanMixers{i} = allProcessedPhaseData{i}.meanMixers;
%    allMeanPulseMixers{i} = allProcessedPhaseData{i}.meanPulseMixers;
    
end



bpmTraceS1 = bpmTraceS{1};
bpmTraceS2 = bpmTraceS{2};
bpmTraceS3 = bpmTraceS{3};

bpmTraceH1 = bpmTraceH{1};
bpmTraceH2 = bpmTraceH{2};
bpmTraceH3 = bpmTraceH{3};

np = 200;

alignedH1 = getAligned(bpmTraceH1,bpmTraceS1,np);
[n,m] = size(alignedH1);

for i=1:n
   figure(101);
   plot(alignedH1(i,:),'-k');
   hold on;
   
%   figure(102)
%   plot(bpmTraceH1(i,:));
%   hold on;
   
end



figure(101);
alignedH2 = getAligned(bpmTraceH2,bpmTraceS2,np);
[n,m] = size(alignedH2);
for i=1:n
   figure(101);
   plot(alignedH2(i,:),'-b');
   hold on;
   
%   figure(102)
%   plot(bpmTraceH1(i,:));
%   hold on;
end

alignedH3 = getAligned(bpmTraceH3,bpmTraceS3,np);
[n,m] = size(alignedH3);

for i=1:n
   figure(101);
   plot(alignedH3(i,:),'-r');
   hold on;
   
%   figure(102)
%   plot(bpmTraceH1(i,:));
%   hold on;
   
end

alignedH4 = getAligned(bpmTraceH{4},bpmTraceS{4},np);
alignedH5 = getAligned(bpmTraceH{5},bpmTraceS{5},np);
alignedH6 = getAligned(bpmTraceH{6},bpmTraceS{6},np);




figure(104)

plot(mean(alignedH1))
hold on;
plot(mean(alignedH2))
plot(mean(alignedH3))


figure(105)

plot(nanmean(alignedH1),'LineWidth',2)
hold on;
plot(nanmean(alignedH2),'LineWidth',2)
plot(nanmean(alignedH3),'LineWidth',2)


%figure(106)
plot(nanmean(alignedH4),'--','LineWidth',2)
hold on;
plot(nanmean(alignedH5),'--','LineWidth',2)
plot(nanmean(alignedH6),'--','LineWidth',2)

legend(dataSetNames);

figure(106)

plot(nanmean(alignedH2) - nanmean(alignedH1),'LineWidth',2)
hold on;
plot(nanmean(alignedH3) - nanmean(alignedH1),'LineWidth',2)
plot(nanmean(alignedH3) - nanmean(alignedH2),'LineWidth',2)

figure(107)

plot(nanmean(alignedH5) - nanmean(alignedH4),'LineWidth',2)
hold on;
plot(nanmean(alignedH6) - nanmean(alignedH4),'LineWidth',2)
plot(nanmean(alignedH6) - nanmean(alignedH5),'LineWidth',2)
