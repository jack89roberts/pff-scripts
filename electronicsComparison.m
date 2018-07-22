dataSets = {'20141215_1323_ElectronicsNominal', '20141215_1117_ElectronicsSwapped'};
nRMS = 20; % used for alignment

%%
close all;

processedData = loadProcessedData(dataSets);
nDataSets = length(dataSets);
pulseLength = 620;
usePortion = 0.5;
removeThreshold = 0.02; % remove pulses that are removeThreshold % away from max on diode

pulsesToUse = {[1:197 385:502], [1:502]};

mix1Phase = cell(1,nDataSets);
mix2Phase = cell(1,nDataSets);
mix1Diode = cell(1,nDataSets);
mix2Diode = cell(1,nDataSets);
mix1Mixer = cell(1,nDataSets);
mix2Mixer = cell(1,nDataSets);
meanMix1Phase = cell(1,nDataSets);
meanMix2Phase = cell(1,nDataSets);
meanMix1Diode = cell(1,nDataSets);
meanMix2Diode = cell(1,nDataSets);
meanMix1Mixer = cell(1,nDataSets);
meanMix2Mixer = cell(1,nDataSets);
stdMix1Phase = cell(1,nDataSets);
stdMix2Phase = cell(1,nDataSets);
stdMix1Diode = cell(1,nDataSets);
stdMix2Diode = cell(1,nDataSets);
stdMix1Mixer = cell(1,nDataSets);
stdMix2Mixer = cell(1,nDataSets);

for i=1:nDataSets

    %% extract, align signals

    % extract signals from processed data
    frascatiPhases = processedData{i}.phaseData.frascatiPhases;
    frascatiDiodes = processedData{i}.phaseData.frascatiDiodes;
    frascatiMixers = processedData{i}.phaseData.frascatiMixers;
    frascatiTimePerSample = processedData{i}.phaseData.frascatiTimePerSample;

    mix2Phase{i} = squeeze(frascatiPhases(2,pulsesToUse{i},:));
    mix1Phase{i} = squeeze(frascatiPhases(1,pulsesToUse{i},:));
    mix2Diode{i} = squeeze(frascatiDiodes(2,pulsesToUse{i},:));
    mix1Diode{i} = squeeze(frascatiDiodes(1,pulsesToUse{i},:));
    mix2Mixer{i} = squeeze(frascatiMixers(2,pulsesToUse{i},:));
    mix1Mixer{i} = squeeze(frascatiMixers(1,pulsesToUse{i},:));
    
    % align signals
    nPoints = round((pulseLength./frascatiTimePerSample)+ 4.*nRMS);
    [mix2Diode{i},mix2Others] = getAligned2(mix2Diode{i},nRMS,nPoints,1,{mix2Mixer{i},mix2Phase{i}});
    [mix1Diode{i},mix1Others] = getAligned2(mix1Diode{i},nRMS,nPoints,1,{mix1Mixer{i},mix1Phase{i}});
    mix2Mixer{i} = mix2Others{1};
    mix1Mixer{i} = mix1Others{1};
    mix2Phase{i} = mix2Others{2};
    mix1Phase{i} = mix1Others{2};

    %% find sample ranges to use to calculate means etc. (pulse starts at nRMS)
    pulseSamplesFrascati = pulseLength./frascatiTimePerSample;
    startFrascati = (nRMS+1) + ( (1-usePortion)./2).*pulseSamplesFrascati;
    endFrascati = (nRMS+1) + ( (1+usePortion)./2).*pulseSamplesFrascati;
    rangeFrascati = round(startFrascati):round(endFrascati);

    %% remove bad pulses
    transmission = nanmean(mix2Diode{i}(:,rangeFrascati),2);
    isGoodPulse = transmission-min(transmission) < removeThreshold;

    mix2Phase{i} = mix2Phase{i}(isGoodPulse,:);
    mix2Diode{i} = mix2Diode{i}(isGoodPulse,:);
    mix2Mixer{i} = mix2Mixer{i}(isGoodPulse,:);
    mix1Phase{i} = mix1Phase{i}(isGoodPulse,:);
    mix1Diode{i} = mix1Diode{i}(isGoodPulse,:);
    mix1Mixer{i} = mix1Mixer{i}(isGoodPulse,:);

    %% calculate means
    meanMix1Phase{i} = nanmean(mix1Phase{i}(:,rangeFrascati),2);
    meanMix2Phase{i} = nanmean(mix2Phase{i}(:,rangeFrascati),2);
    meanMix1Diode{i} = nanmean(mix1Diode{i}(:,rangeFrascati),2);
    meanMix2Diode{i} = nanmean(mix2Diode{i}(:,rangeFrascati),2);
    meanMix1Mixer{i} = nanmean(mix1Mixer{i}(:,rangeFrascati),2);
    meanMix2Mixer{i} = nanmean(mix2Mixer{i}(:,rangeFrascati),2);
    
    stdMix1Phase{i} = nanstd(mix1Phase{i}(:,rangeFrascati),0,1);
    stdMix2Phase{i} = nanstd(mix2Phase{i}(:,rangeFrascati),0,1);
    stdMix1Diode{i} = nanstd(mix1Diode{i}(:,rangeFrascati),0,1);
    stdMix2Diode{i} = nanstd(mix2Diode{i}(:,rangeFrascati),0,1);
    stdMix1Mixer{i} = nanstd(mix1Mixer{i}(:,rangeFrascati),0,1);
    stdMix2Mixer{i} = nanstd(mix2Mixer{i}(:,rangeFrascati),0,1);

end

%% drift subtract
nAvg = 5;
driftMix1Phase = cell(1,nDataSets);
driftMix2Phase = cell(1,nDataSets);
driftMix1Diode = cell(1,nDataSets);
driftMix2Diode = cell(1,nDataSets);
driftMix1Mixer = cell(1,nDataSets);
driftMix2Mixer = cell(1,nDataSets);

for i=1:nDataSets
    nPulses = length(meanMix1Phase{i});
    
    driftMix1Phase{i} = NaN*ones(1,nPulses);
    driftMix2Phase{i} = NaN*ones(1,nPulses);
    driftMix1Diode{i} = NaN*ones(1,nPulses);
    driftMix2Diode{i} = NaN*ones(1,nPulses);
    driftMix1Mixer{i} = NaN*ones(1,nPulses);
    driftMix2Mixer{i} = NaN*ones(1,nPulses);
    
    for n = nAvg:nPulses
    
        driftMix1Phase{i}(n) = meanMix1Phase{i}(n)-nanmean(meanMix1Phase{i}((n-nAvg+1):n));
        driftMix2Phase{i}(n) = meanMix2Phase{i}(n)-nanmean(meanMix2Phase{i}((n-nAvg+1):n));
        driftMix1Diode{i}(n) = meanMix1Diode{i}(n)-nanmean(meanMix1Diode{i}((n-nAvg+1):n));
        driftMix2Diode{i}(n) = meanMix2Diode{i}(n)-nanmean(meanMix2Diode{i}((n-nAvg+1):n));
        driftMix1Mixer{i}(n) = meanMix1Mixer{i}(n)-nanmean(meanMix1Mixer{i}((n-nAvg+1):n));
        driftMix2Mixer{i}(n) = meanMix2Mixer{i}(n)-nanmean(meanMix2Mixer{i}((n-nAvg+1):n));
            
    end
    
end

%%

figure;
subplot(1,3,1);
plot(meanMix1Diode{1})
hold all;
plot(meanMix2Diode{1})
plot(meanMix1Diode{2})
plot(meanMix2Diode{2})
legend('MIX1 - MON1','MIX2 - MON2','MIX1 - MON2','MIX2 - MON1')
title('MEAN DIODE')

subplot(1,3,2);
plot(meanMix1Mixer{1})
hold all;
plot(meanMix2Mixer{1})
plot(meanMix1Mixer{2})
plot(meanMix2Mixer{2})
legend('MIX1 - MON1','MIX2 - MON2','MIX1 - MON2','MIX2 - MON1')
title('MEAN MIXER')

subplot(1,3,3);
plot(meanMix1Phase{1})
hold all;
plot(meanMix2Phase{1})
plot(meanMix1Phase{2})
plot(meanMix2Phase{2})
legend('MIX1 - MON1','MIX2 - MON2','MIX1 - MON2','MIX2 - MON1')
title('MEAN PHASE')


figure;
subplot(1,3,1);
plot(stdMix1Diode{1})
hold all;
plot(stdMix2Diode{1})
plot(stdMix1Diode{2})
plot(stdMix2Diode{2})
legend('MIX1 - MON1','MIX2 - MON2','MIX1 - MON2','MIX2 - MON1')
title('STD DIODE')

subplot(1,3,2);
plot(stdMix1Mixer{1})
hold all;
plot(stdMix2Mixer{1})
plot(stdMix1Mixer{2})
plot(stdMix2Mixer{2})
legend('MIX1 - MON1','MIX2 - MON2','MIX1 - MON2','MIX2 - MON1')
title('STD MIXER')

subplot(1,3,3);
plot(stdMix1Phase{1})
hold all;
plot(stdMix2Phase{1})
plot(stdMix1Phase{2})
plot(stdMix2Phase{2})
legend('MIX1 - MON1','MIX2 - MON2','MIX1 - MON2','MIX2 - MON1')
title('STD PHASE')

figure;
subplot(1,3,1);
plot(driftMix1Diode{1})
hold all;
plot(driftMix2Diode{1})
plot(driftMix1Diode{2})
plot(driftMix2Diode{2})
legend('MIX1 - MON1','MIX2 - MON2','MIX1 - MON2','MIX2 - MON1')
title('MEAN DIODE (DRIFT SUBTRACT)')

subplot(1,3,2);
plot(driftMix1Mixer{1})
hold all;
plot(driftMix2Mixer{1})
plot(driftMix1Mixer{2})
plot(driftMix2Mixer{2})
legend('MIX1 - MON1','MIX2 - MON2','MIX1 - MON2','MIX2 - MON1')
title('MEAN MIXER (DRIFT SUBTRACT)')

subplot(1,3,3);
plot(driftMix1Phase{1})
hold all;
plot(driftMix2Phase{1})
plot(driftMix1Phase{2})
plot(driftMix2Phase{2})
legend('MIX1 - MON1','MIX2 - MON2','MIX1 - MON2','MIX2 - MON1')
title('MEAN PHASE (DRIFT SUBTRACT)')
