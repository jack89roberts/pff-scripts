close all;
%%
% CTFData = loadMergedData('20141204_1617_FF_NewFMW_Gate365to465_GainK1_-0_GainK2_0');
% CTFData = loadMergedData('20150708_1222_FFGain_63_Interleaved');
CTFData = loadMergedData('20151117_1520_Gain_700_R56_0.1');

sampRange = 230:260;
sampRange = 190:210;
sampRange = 180:250;

startSamp = 107;
endSamp = 324;
%%
try
    ct608 = extractCTFSignalFromMergedData('CT.SVBPI0608H.Samples.samples.value',CTFData);
    ct608S = extractCTFSignalFromMergedData('CT.SVBPI0608S.Samples.samples.value',CTFData);
catch
   ct608 = extractCTFSignalFromMergedData('CT_SVBPI0608H.Samples.samples.value',CTFData); 
   ct608S = extractCTFSignalFromMergedData('CT_SVBPI0608S.Samples.samples.value',CTFData);
end
ct608 = removeDuplicatePulses(ct608);
ct608 = removeBadPulses(ct608,sampRange);
ct608 = ct608/(1000*-0.605);

timeAxis = 5.2*((1:length(meanAlong608))-startSamp);
%%
meanAlong608 = nanmean(ct608);
meanAlong608 = meanAlong608-nanmean(meanAlong608(sampRange));
stdAlong608 = nanstd(ct608);
mean608 = nanmean(ct608(:,sampRange),2);

mean608 = mean608-nanmean(mean608);
stdMean608 = nanstd(mean608);

figure;
plot(timeAxis,meanAlong608,'b','LineWidth',2);
title('Energy Along Pulse')
xlabel('Time [ns]')
ylabel('Relative Energy Offset')
xlim([timeAxis(startSamp+8) timeAxis(endSamp-1)])
format_plots;

figure;
plot(timeAxis,stdAlong608,'b','LineWidth',2);
title('Energy Jitter Along Pulse')
xlabel('Time [ns]')
ylabel('Relative Energy Jitter')
xlim([timeAxis(startSamp+8) timeAxis(endSamp-1)])
format_plots;

figure;
plot(mean608,'b','LineWidth',2)
xlabel('Time [Pulse No.]');
ylabel('Relative Energy Offset')
title(sprintf('Mean Energy Variations (std = %.1d)',stdMean608))
format_plots;

meanS = nanmean(ct608S);
figure;
plot(meanS);

return;
%%
datDir = '/home/jack/PhaseFeedforward/CTFData/201511';
sampRange = 180:250;

datSets = dir(datDir);
datSets = {datSets.name};
datSets = datSets(3:end);
nDatSets = length(datSets);
jit608 = NaN(1,nDatSets);
for i=1:nDatSets
    fprintf('Dat set %d of %d: %s....\n',i,nDatSets,datSets{i})
    
    CTFData = loadMergedData(datSets{i});
    try
        ct608 = extractCTFSignalFromMergedData('CT.SVBPI0608H.Samples.samples.value',CTFData);
        ct608S = extractCTFSignalFromMergedData('CT.SVBPI0608S.Samples.samples.value',CTFData);
    catch
        try
            ct608 = extractCTFSignalFromMergedData('CT_SVBPI0608H.Samples.samples.value',CTFData);
            ct608S = extractCTFSignalFromMergedData('CT_SVBPI0608S.Samples.samples.value',CTFData);
        catch
            continue;
        end
    end
    ct608 = removeDuplicatePulses(ct608);
    ct608 = removeBadPulses(ct608,sampRange);
    
    mean608 = nanmean(ct608(:,sampRange),2);
    mean608 = mean608-nanmean(mean608);
    mean608 = mean608./(1000*-0.605);
    
    jit608(i) = nanstd(mean608);
    fprintf('%s: %.2d\n',datSets{i},jit608(i))
    fprintf('--------------------------------------------------\n');

end
fprintf('\n')
fprintf('--------------------------------------------------\n');
fprintf('\n')
for i=1:nDatSets
    if ~isnan(jit608(i))
        fprintf('%s: %.2d\n',datSets{i},jit608(i))
    end 
end

figure;
plot(jit608);

%%
datDir = '/home/jack/PhaseFeedforward/Analysis/201511';
datSets = dir([datDir '/*.mat']);
datSets = {datSets.name};
nDatSets = length(datSets);
corrMon2En = NaN(1,nDatSets);
corrMon2En_err = NaN(1,nDatSets);
corrMon3En =  NaN(1,nDatSets);
corrMon3En_err = NaN(1,nDatSets);
corrMon2Mon3 = NaN(1,nDatSets);
corrMon2Mon3_err = NaN(1,nDatSets);
nPulses = NaN(1,nDatSets);
for i=1:nDatSets
    tmpDat = load([datDir '/' datSets{i}]);
    nPulses(i) = tmpDat.nPulses;
    
    try
        corrMon2Mon3(i) = tmpDat.corrMeanMix2Mix3;
        corrMon2Mon3_err(i) = tmpDat.corrMeanMix2Mix3_err;
    catch
    end
    
    try
        corrMon2En(i) = tmpDat.corrMon2_BPMH(14);
        corrMon2En_err(i) = tmpDat.corrMon2_BPMH_err(14);
    catch
    end
    
    try
        corrMon3En(i) = tmpDat.corrMon3_BPMH(14);
        corrMon3En_err(i) = tmpDat.corrMon3_BPMH_err(14);
    catch
    end
    
    fprintf('%d) %s: %d pulses, 2-E: %.2f, 3-E: %.2f, 2-3: %.2f\n',i,datSets{i},nPulses(i),corrMon2En(i),corrMon3En(i),corrMon2Mon3(i))
end


%71) 20151117_1754_Gain_750_R56_0.1_Ilvd_Even.mat: 300 pulses, 2-E: 0.19, 3-E: 0.35, 2-3: 0.81
%72) 20151117_1754_Gain_750_R56_0.1_Ilvd_Odd.mat: 300 pulses, 2-E: 0.19, 3-E: 0.39, 2-3: 0.32
i=71;
tmpDat = load([datDir '/' datSets{i}]);
figure;
ct608 = tmpDat.meanBPMH{14};
ct608 = ct608 - nanmean(ct608);
ct608 = ct608/(1000*0.605);
mon2 = tmpDat.meanPulsePhase(2,:);
mon2 = mon2-nanmean(mon2);
plot(ct608,mon2,'bo','MarkerFaceColor','b')
xlabel('Relative Energy Offset')
ylabel('Upstream Phase [degrees]')
legend(sprintf('correlation: %.2f%c%.2f',corrMon2En(i),char(177),corrMon2En_err(i)))
title('Upstream Phase vs. Energy')
format_plots;