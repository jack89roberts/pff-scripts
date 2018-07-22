close all;
%% Inputs

dataFile = '/home/jack/PhaseFeedforward/CTFData/201412/allforstability/20141215_1411/processed/20141215_1411.mat';
saveDir = '/home/jack/PhaseFeedforward/Analysis/CorrelationAnalysis/tmp/';

% delays of crates (must be >= 0)
frascatiDelay = 0; % cfv-2001-ccaos5
petsDelay = 0; % cfi-2010-ccaos3
bpm1Delay = 1;%2; % dctfumal
bpm2Delay = 1; % dctfumad
bpm3Delay = 1; % dctfumar
bpm4Delay = 1; % cfc-2010-bpm

pulseLength = 620;
usePortion = 0.5; % fraction of pulse to use when calculating correlations etc. (taken in middle)
pulsesToUse = 265:606;%-1;%265:606; % set to -1 to use all pulses

savePlots = false;

%% Load data

load(dataFile);

if (savePlots && ~exist(saveDir,'dir'))
    mkdir(saveDir);
end

%% extract signals

% some parameters, shouldn't need to change them
nRMS = 20; % used for alignment
bpmCCTimePerSample = 5.208;
bpmCTTimePerSample = 10.417;
%removeThreshold = 0.02; % remove pulses that are removeThreshold % away from max on diode

% extract signals from processed data
frascatiPhases = processedData.phaseData.frascatiPhases;
frascatiDiodes = processedData.phaseData.frascatiDiodes;
frascatiTimePerSample = processedData.phaseData.frascatiTimePerSample;
petsPhase = processedData.phaseData.petsPhase;
petsQ = processedData.phaseData.petsQ;
petsTimePerSample = processedData.phaseData.petsTimePerSample;
bpmNames = processedData.bpmData.bpmNames;


mon2Phase = squeeze(frascatiPhases(2,:,:));
mon3Phase = squeeze(frascatiPhases(3,:,:));
mon2Diode = squeeze(frascatiDiodes(2,:,:));
mon3Diode = squeeze(frascatiDiodes(3,:,:));

try
    eval(sprintf('processedData.bpmData.%sV.Samples.samples;',bpmNames{1}));
    useVertical = 1;
catch
    useVertical = 0;
end

nBPMs = length(bpmNames);
bpmS = cell(1,nBPMs);
bpmH = cell(1,nBPMs);
bpmV = cell(1,nBPMs);
bpmNSamples = NaN*ones(1,nBPMs);
for i=1:nBPMs
    bpmS{i} = eval(sprintf('processedData.bpmData.%sS.Samples.samples;',bpmNames{i}));
    bpmH{i} = eval(sprintf('processedData.bpmData.%sH.Samples.samples;',bpmNames{i}));
    if (useVertical)
        bpmV{i} = eval(sprintf('processedData.bpmData.%sV.Samples.samples;',bpmNames{i}));
    end
    [~,bpmNSamples(i)] = size(bpmS{i});
end

bpmCL502Index = find(ismember(bpmNames,'CL.SVBPM0502')); 
bpmCT285Index = find(ismember(bpmNames,'CT.SVBPM0285'));
useCL502 = (bpmCL502Index > 0);
useCT285 = (bpmCT285Index > 0);

%% remove (NaN) any duplicate pulses

mon2Phase = removeDuplicatePulses(mon2Phase);
mon3Phase = removeDuplicatePulses(mon3Phase);
mon2Diode = removeDuplicatePulses(mon2Diode);
mon3Diode = removeDuplicatePulses(mon3Diode);
petsPhase = removeDuplicatePulses(petsPhase);
petsQ = removeDuplicatePulses(petsQ);

for i=1:nBPMs
    bpmS{i} = removeDuplicatePulses(bpmS{i});
    bpmH{i} = removeDuplicatePulses(bpmH{i});
    if (useVertical)
        bpmV{i} = removeDuplicatePulses(bpmV{i});
    end
end

%% shift oasis signals to line them up with BPMs 
% deletes appropriate number of entries at start/end of arrays to align
% pulses.
maxDelay = max([frascatiDelay petsDelay bpm1Delay bpm2Delay bpm3Delay bpm4Delay]);

% Frascati
if (frascatiDelay > 0)
    mon2Phase(1:frascatiDelay, :) = [];
    mon3Phase(1:frascatiDelay, :) = [];
    mon2Diode(1:frascatiDelay, :) = [];
    mon3Diode(1:frascatiDelay, :) = [];
end    
if (frascatiDelay<maxDelay && maxDelay>0)
    mon2Phase((end-maxDelay+frascatiDelay+1):end, :) = [];
    mon3Phase((end-maxDelay+frascatiDelay+1):end, :) = [];
    mon2Diode((end-maxDelay+frascatiDelay+1):end, :) = [];
    mon3Diode((end-maxDelay+frascatiDelay+1):end, :) = [];
end

% PETS
if (petsDelay > 0)
    petsPhase(1:petsDelay, :) = [];
    petsQ(1:petsDelay, :) = [];
end    
if (petsDelay<maxDelay && maxDelay>0)
    petsPhase((end-maxDelay+petsDelay+1):end, :) = [];
    petsQ((end-maxDelay+petsDelay+1):end, :) = [];
end

% BPMs
% load file with device crate information
filename = 'devices/deviceCrates.dat';
delimiter = ',';
formatSpec = '%s%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
devNames = dataArray{:, 1};
crateNames = dataArray{:, 2};

for i=1:nBPMs
    tmpBPMCrateIndex = find(ismember(devNames,bpmNames{i}));
    if (tmpBPMCrateIndex > 0)
        tmpCrateName = crateNames{tmpBPMCrateIndex};
    else
        error('Could not find %s in list of device crates.\n',bpmNames{i});
    end
    switch (tmpCrateName)
        case 'dctfumal'
            tmpDelay = bpm1Delay;
        case 'dctfumad'
            tmpDelay = bpm2Delay;
        case 'dctfumar'
            tmpDelay = bpm3Delay;
        case 'cfc-2010-bpm'
            tmpDelay = bpm4Delay;
        otherwise
            error('No case exists for BPM %s on crate %s\n',bpmNames{i},tmpCrateName);
    end
    
    if (tmpDelay > 0)
        bpmS{i}(1:tmpDelay, :) = [];
        bpmH{i}(1:tmpDelay, :) = [];
        if (useVertical)
            bpmV{i}(1:tmpDelay, :) = [];
        end
    end    
    if (tmpDelay<maxDelay && maxDelay>0)
        bpmS{i}((end-maxDelay+tmpDelay+1):end, :) = [];
        bpmH{i}((end-maxDelay+tmpDelay+1):end, :) = [];
        if (useVertical)
            bpmV{i}((end-maxDelay+tmpDelay+1):end, :) = [];
        end
    end

end

%% align signals
nPoints = round((pulseLength./frascatiTimePerSample)+ 4.*nRMS);
[mon2Diode,mon2Phase] = getAligned2(mon2Diode,nRMS,nPoints,1,{mon2Phase});
[mon3Diode,mon3Phase] = getAligned2(mon3Diode,nRMS,nPoints,1,{mon3Phase});
mon2Phase = mon2Phase{1};
mon3Phase = mon3Phase{1};
nPoints = round((pulseLength./petsTimePerSample)+ 4.*nRMS);
[petsQ,petsPhase] = getAligned2(petsQ,nRMS,nPoints,1,{petsPhase});
petsPhase = petsPhase{1};


for i=1:nBPMs
    tmpBPMS = bpmS{i};
    tmpBPMH = bpmH{i};
    tmpBPMV = bpmV{i};
    
    if (bpmNSamples(i) > 400)
        tmpTimePerSample = bpmCCTimePerSample;
    else
        tmpTimePerSample = bpmCTTimePerSample;
    end
    
    nPoints = round((pulseLength./tmpTimePerSample) + 4.*nRMS);
    if (useVertical)
        [tmpBPMS,tmpOthers] = getAligned2(tmpBPMS,nRMS,nPoints,1,{tmpBPMH, tmpBPMV});
        tmpBPMH = tmpOthers{1};
        tmpBPMV = tmpOthers{2};
    else
        [tmpBPMS,tmpOthers] = getAligned2(tmpBPMS,nRMS,nPoints,1,{tmpBPMH});
        tmpBPMH = tmpOthers{1};
    end
    bpmS{i} = tmpBPMS;
    bpmH{i} = tmpBPMH;
    if (useVertical)
        bpmV{i} = tmpBPMV;
    end
end

%% find sample ranges to use to calculate means etc. (pulse starts at nRMS)
pulseSamplesFrascati = pulseLength./frascatiTimePerSample;
startFrascati = (nRMS+1) + ( (1-usePortion)./2).*pulseSamplesFrascati;
endFrascati = (nRMS+1) + ( (1+usePortion)./2).*pulseSamplesFrascati;
rangeFrascati = round(startFrascati):round(endFrascati);

pulseSamplesPETS = pulseLength./petsTimePerSample;
startPETS = (nRMS+1) + ( (1-usePortion)./2).*pulseSamplesPETS;
endPETS = (nRMS+1) + ( (1+usePortion)./2).*pulseSamplesPETS;
rangePETS = round(startPETS):round(endPETS);

pulseSamplesBPMCC = pulseLength./bpmCCTimePerSample;
startBPMCC = (nRMS+1) + ( (1-usePortion)./2).*pulseSamplesBPMCC;
endBPMCC = (nRMS+1) + ( (1+usePortion)./2).*pulseSamplesBPMCC;
rangeBPMCC = round(startBPMCC):round(endBPMCC);

pulseSamplesBPMCT = pulseLength./bpmCTTimePerSample;
startBPMCT = (nRMS+1) + ( (1-usePortion)./2).*pulseSamplesBPMCT;
endBPMCT = (nRMS+1) + ( (1+usePortion)./2).*pulseSamplesBPMCT;
rangeBPMCT = round(startBPMCT):round(endBPMCT);

% check ranges
figure;
plot(mon2Diode(1,:));
hold all;
plot(mon3Diode(1,:));
minY = min([mon2Diode(1,:) mon3Diode(1,:)]);
maxY = max([mon2Diode(1,:) mon3Diode(1,:)]);
plot([startFrascati startFrascati], [minY maxY],'k','LineWidth',2);
plot([endFrascati endFrascati], [minY maxY],'k','LineWidth',2);
title('SELECTED RANGE - FRASCATI MONITORS');
xlabel('Sample No.');
ylabel('Output');

figure;
plot(petsQ(1,:));
hold all;
minY = min(petsQ(1,:));
maxY = max(petsQ(1,:));
plot([startPETS startPETS], [minY maxY],'k','LineWidth',2);
plot([endPETS endPETS], [minY maxY],'k','LineWidth',2);
title('SELECTED RANGE - PETS');
xlabel('Sample No.');
ylabel('Output');

figBPMCT = figure();
figBPMCC = figure();
minCC = 1000000;
maxCC = -1000000;
minCT = 1000000;
maxCT = -1000000;
for i=1:nBPMs
    tmpMin = min(bpmS{i}(1,:));
    tmpMax = max(bpmS{i}(1,:));
    if (bpmNSamples(i) > 400)
        figure(figBPMCC);
        if (tmpMin < minCC)
            minCC = tmpMin;
        end
        if (tmpMax > maxCC)
            maxCC = tmpMax;
        end
    else
        figure(figBPMCT);
        if (tmpMin < minCT)
            minCT = tmpMin;
        end
        if (tmpMax > maxCT)
            maxCT = tmpMax;
        end
    end
    
    plot(bpmS{i}(1,:));
    hold all;
end
figure(figBPMCC);
plot([startBPMCC startBPMCC], [minCC maxCC],'k','LineWidth',2);
plot([endBPMCC endBPMCC], [minCC maxCC],'k','LineWidth',2);
title('SELECTED RANGE - BPM CC');
xlabel('Sample No.');
ylabel('Output');
figure(figBPMCT);
plot([startBPMCT startBPMCT], [minCT maxCT],'k','LineWidth',2);
plot([endBPMCT endBPMCT], [minCT maxCT],'k','LineWidth',2);
title('SELECTED RANGE - BPM CT');
xlabel('Sample No.');
ylabel('Output');


%% remove bad pulses

% 3 sigma removal
mon2Diode = removeBadPulses(mon2Diode,rangeFrascati);
mon2Phase = removeBadPulses(mon2Phase,rangeFrascati);
mon3Diode = removeBadPulses(mon3Diode,rangeFrascati);
mon3Phase = removeBadPulses(mon3Phase,rangeFrascati);
petsQ = removeBadPulses(petsQ,rangePETS);
petsPhase = removeBadPulses(petsPhase,rangePETS);

for i=1:nBPMs
    if (bpmNSamples(i) > 400)
        tmpRange = rangeBPMCC;
    else
        tmpRange = rangeBPMCT;
    end
    
    bpmS{i} = removeBadPulses(bpmS{i},tmpRange);
    bpmH{i} = removeBadPulses(bpmH{i},tmpRange);
    if (useVertical)
        bpmV{i} = removeBadPulses(bpmV{i},tmpRange);
    end
end

% % based on Mon2 diode only (remove all pulses removeThreshold away from max)
% transmission = nanmean(mon2Diode(:,rangeFrascati),2);
% isGoodPulse = transmission-min(transmission) < removeThreshold;
% 
% mon2Phase = mon2Phase(isGoodPulse,:);
% mon2Diode = mon2Diode(isGoodPulse,:);
% mon3Phase = mon3Phase(isGoodPulse,:);
% mon3Diode = mon3Diode(isGoodPulse,:);
% petsPhase = petsPhase(isGoodPulse,:);
% petsQ = petsQ(isGoodPulse,:);
% 
% for i=1:nBPMs
%     bpmS{i} = bpmS{i}(isGoodPulse,:);
%     bpmH{i} = bpmH{i}(isGoodPulse,:);
%     if (useVertical)
%         bpmV{i} = bpmV{i}(isGoodPulse,:);
%     end
% end


%% calculate means
if (length(pulsesToUse)==1 && (pulsesToUse==-1))
    [nPulses,~] = size(mon2Phase);
    pulsesToUse = 1:nPulses;
end
nPulses = length(pulsesToUse);

meanMon2Phase = nanmean(mon2Phase(pulsesToUse,rangeFrascati),2);
meanMon3Phase = nanmean(mon3Phase(pulsesToUse,rangeFrascati),2);
meanPETSPhase = nanmean(petsPhase(pulsesToUse,rangeFrascati),2);
meanMon2Diode = nanmean(mon2Diode(pulsesToUse,rangeFrascati),2);
meanMon3Diode = nanmean(mon3Diode(pulsesToUse,rangeFrascati),2);
meanPETSQ = nanmean(petsQ(pulsesToUse,rangeFrascati),2);

meanBPMS = NaN*ones(nBPMs,nPulses);
meanBPMH = NaN*ones(nBPMs,nPulses);
meanBPMV = NaN*ones(nBPMs,nPulses);
for i=1:nBPMs
    if (bpmNSamples(i) > 400)
        tmpRange = rangeBPMCC;
    else
        tmpRange = rangeBPMCT;
    end
    
    meanBPMS(i,:) = nanmean(bpmS{i}(pulsesToUse,tmpRange),2);
    meanBPMH(i,:) = nanmean(bpmH{i}(pulsesToUse,tmpRange),2);
    if (useVertical)
        meanBPMV(i,:) = nanmean(bpmV{i}(pulsesToUse,tmpRange),2);
    end
end

if (useCT285)
    meanBPMCT285S = meanBPMS(bpmCT285Index,:);
    meanBPMCT285H = meanBPMH(bpmCT285Index,:);
    if (useVertical)
        meanBPMCT285V = meanBPMV(bpmCT285Index,:);
    end
end
if (useCL502)
    meanBPMCL502S = meanBPMS(bpmCL502Index,:);
    meanBPMCL502H = meanBPMH(bpmCL502Index,:);
    if (useVertical)
        meanBPMCL502V = meanBPMV(bpmCL502Index,:);
    end
end

figure;
plot(meanMon2Diode);
hold all;
plot(meanMon3Diode);

figure;
plot(meanBPMS')

%% calculate correlations

isNaN2 = isnan(meanMon2Phase);
isNaN3 = isnan(meanMon3Phase);
isGoodPulse = ~(isNaN2 | isNaN3);
goodMon2 = meanMon2Phase(isGoodPulse);
goodMon3 = meanMon3Phase(isGoodPulse);
corrMon2Mon3 = corrcoef(goodMon2,goodMon3);
corrMon2Mon3 = corrMon2Mon3(1,2);
figure;
scatter(goodMon2,goodMon3);
title(sprintf('Mon2-Mon3: %.2f',corrMon2Mon3))
xlabel('Mon2 Phase [degrees]');
ylabel('Mon3 Phase [degrees]');
if (savePlots)
    saveStr = sprintf('%scorr_Mon2_Mon3',saveDir);
    print([saveStr '.png'],'-dpng');
    savefig([saveStr '.fig']);
end

isNaNPETS = isnan(meanPETSPhase);
isGoodPulse = ~(isNaNPETS | isNaN3);
goodPETS = meanPETSPhase(isGoodPulse);
goodMon3 = meanMon3Phase(isGoodPulse);
corrMon3PETS = corrcoef(goodMon3,goodPETS);
corrMon3PETS = corrMon3PETS(1,2);
figure;
scatter(goodMon3,goodPETS);
title(sprintf('Mon3-PETS: %.2f',corrMon3PETS))
xlabel('Mon3 Phase [degrees]');
ylabel('PETS Phase [degrees]');
if (savePlots)
    saveStr = sprintf('%scorr_Mon3_PETS',saveDir);
    print([saveStr '.png'],'-dpng');
    savefig([saveStr '.fig']);
end

isNaNPhase = isnan(meanMon2Phase);
isNaNDiode = isnan(meanMon2Diode);
isGoodPulse = ~(isNaNPhase | isNaNDiode);
goodPhase = meanMon2Phase(isGoodPulse);
goodDiode = meanMon2Diode(isGoodPulse);
corrMon2Diode = corrcoef(goodPhase,goodDiode);
corrMon2Diode = corrMon2Diode(1,2);
figure;
scatter(goodPhase,goodDiode);
title(sprintf('Mon2 Phase-Diode: %.2f',corrMon2Diode))
xlabel('Mon2 Phase [degrees]');
ylabel('Mon2 Diode [V]');
if (savePlots)
    saveStr = sprintf('%scorr_Mon2_Diode',saveDir);
    print([saveStr '.png'],'-dpng');
    savefig([saveStr '.fig']);
end

isNaNPhase = isnan(meanMon3Phase);
isNaNDiode = isnan(meanMon3Diode);
isGoodPulse = ~(isNaNPhase | isNaNDiode);
goodPhase = meanMon3Phase(isGoodPulse);
goodDiode = meanMon3Diode(isGoodPulse);
corrMon3Diode = corrcoef(goodPhase,goodDiode);
corrMon3Diode = corrMon3Diode(1,2);
figure;
scatter(goodPhase,goodDiode);
title(sprintf('Mon3 Phase-Diode: %.2f',corrMon3Diode))
xlabel('Mon3 Phase [degrees]');
ylabel('Mon3 Diode [V]');
if (savePlots)
    saveStr = sprintf('%scorr_Mon3_Diode',saveDir);
    print([saveStr '.png'],'-dpng');
    savefig([saveStr '.fig']);
end

if (useCT285)
    isNaNBPMCT285H = isnan(meanBPMCT285H);
    isNaNBPMCT285S = isnan(meanBPMCT285S);
    if (useVertical)
        isNaNBPMCT285V = isnan(meanBPMCT285V);
    end
end
if (useCL502)
    isNaNBPMCL502H = isnan(meanBPMCL502H);
    isNaNBPMCL502S = isnan(meanBPMCL502S);
    if (useVertical)
        isNaNBPMCL502V = isnan(meanBPMCL502V);
    end
end

corrMon3_BPMS = NaN*ones(1,nBPMs);
corrMon3_BPMH = NaN*ones(1,nBPMs);
corrMon3_BPMV = NaN*ones(1,nBPMs);
corrMon2_BPMS = NaN*ones(1,nBPMs);
corrMon2_BPMH = NaN*ones(1,nBPMs);
corrMon2_BPMV = NaN*ones(1,nBPMs);

corrBPMCT285H_BPMH = NaN*ones(1,nBPMs);
corrBPMCT285H_BPMS = NaN*ones(1,nBPMs);
corrBPMCT285H_BPMV = NaN*ones(1,nBPMs);
corrBPMCT285V_BPMH = NaN*ones(1,nBPMs);
corrBPMCT285V_BPMS = NaN*ones(1,nBPMs);
corrBPMCT285V_BPMV = NaN*ones(1,nBPMs);
corrBPMCT285S_BPMH = NaN*ones(1,nBPMs);
corrBPMCT285S_BPMS = NaN*ones(1,nBPMs);
corrBPMCT285S_BPMV = NaN*ones(1,nBPMs);
corrBPMCL502H_BPMH = NaN*ones(1,nBPMs);
corrBPMCL502H_BPMS = NaN*ones(1,nBPMs);
corrBPMCL502H_BPMV = NaN*ones(1,nBPMs);
corrBPMCL502V_BPMH = NaN*ones(1,nBPMs);
corrBPMCL502V_BPMS = NaN*ones(1,nBPMs);
corrBPMCL502V_BPMV = NaN*ones(1,nBPMs);
corrBPMCL502S_BPMH = NaN*ones(1,nBPMs);
corrBPMCL502S_BPMS = NaN*ones(1,nBPMs);
corrBPMCL502S_BPMV = NaN*ones(1,nBPMs);

for i=1:nBPMs
    isNaNBPMH = isnan(meanBPMH(i,:));
    isNaNBPMS = isnan(meanBPMS(i,:));

    % Mon3 - BPMH
    isGoodPulse = ~(isNaNBPMH | isNaN3');
    goodBPM = meanBPMH(i,isGoodPulse);
    goodMon3 = meanMon3Phase(isGoodPulse);
    tmpCorr = corrcoef(goodMon3,goodBPM);
    corrMon3_BPMH(i) = tmpCorr(1,2);

    % Mon3 - BPMS
    isGoodPulse = ~(isNaNBPMS | isNaN3');
    goodBPM = meanBPMS(i,isGoodPulse);
    goodMon3 = meanMon3Phase(isGoodPulse);
    tmpCorr = corrcoef(goodMon3,goodBPM);
    corrMon3_BPMS(i) = tmpCorr(1,2);

    % Mon2 - BPMH
    isGoodPulse = ~(isNaNBPMH | isNaN2');
    goodBPM = meanBPMH(i,isGoodPulse);
    goodMon2 = meanMon2Phase(isGoodPulse);
    tmpCorr = corrcoef(goodMon2,goodBPM);
    corrMon2_BPMH(i) = tmpCorr(1,2);
    
    % Mon2 - BPMS
    isGoodPulse = ~(isNaNBPMS | isNaN2');
    goodBPM = meanBPMS(i,isGoodPulse);
    goodMon2 = meanMon2Phase(isGoodPulse);
    tmpCorr = corrcoef(goodMon2,goodBPM);
    corrMon2_BPMS(i) = tmpCorr(1,2);
        
    if (useVertical)
        isNaNBPMV = isnan(meanBPMV(i,:));

        % Mon3 - BPMV
        isGoodPulse = ~(isNaNBPMV | isNaN3');
        goodBPM = meanBPMV(i,isGoodPulse);
        goodMon3 = meanMon3Phase(isGoodPulse);
        tmpCorr = corrcoef(goodMon3,goodBPM);
        corrMon3_BPMV(i) = tmpCorr(1,2);

        % Mon2 - BPMV
        isGoodPulse = ~(isNaNBPMV | isNaN2');
        goodBPM = meanBPMV(i,isGoodPulse);
        goodMon2 = meanMon2Phase(isGoodPulse);
        tmpCorr = corrcoef(goodMon2,goodBPM);
        corrMon2_BPMV(i) = tmpCorr(1,2);   
    end
    
    if (useCT285)
        % CT285H - BPMH
        isGoodPulse = ~(isNaNBPMH | isNaNBPMCT285H);
        goodBPM = meanBPMH(i,isGoodPulse);
        goodBPMCT285H = meanBPMCT285H(isGoodPulse);
        tmpCorr = corrcoef(goodBPMCT285H,goodBPM);
        corrBPMCT285H_BPMH(i) = tmpCorr(1,2);

        % CT285H - BPMS
        isGoodPulse = ~(isNaNBPMS | isNaNBPMCT285H);
        goodBPM = meanBPMS(i,isGoodPulse);
        goodBPMCT285H = meanBPMCT285H(isGoodPulse);
        tmpCorr = corrcoef(goodBPMCT285H,goodBPM);
        corrBPMCT285H_BPMS(i) = tmpCorr(1,2);

         % CT285S - BPMH
        isGoodPulse = ~(isNaNBPMH | isNaNBPMCT285S);
        goodBPM = meanBPMH(i,isGoodPulse);
        goodBPMCT285S = meanBPMCT285S(isGoodPulse);
        tmpCorr = corrcoef(goodBPMCT285S,goodBPM);
        corrBPMCT285S_BPMH(i) = tmpCorr(1,2);

        % CT285S - BPMS
        isGoodPulse = ~(isNaNBPMS | isNaNBPMCT285S);
        goodBPM = meanBPMS(i,isGoodPulse);
        goodBPMCT285S = meanBPMCT285S(isGoodPulse);
        tmpCorr = corrcoef(goodBPMCT285S,goodBPM);
        corrBPMCT285S_BPMS(i) = tmpCorr(1,2);

        
        if (useVertical)
            isNaNBPMV = isnan(meanBPMV(i,:));
            
            % CT285H - BPMV
            isGoodPulse = ~(isNaNBPMV | isNaNBPMCT285H);
            goodBPM = meanBPMV(i,isGoodPulse);
            goodBPMCT285H = meanBPMCT285H(isGoodPulse);
            tmpCorr = corrcoef(goodBPMCT285H,goodBPM);
            corrBPMCT285H_BPMV(i) = tmpCorr(1,2);

            % CT285V - BPMH
            isGoodPulse = ~(isNaNBPMH | isNaNBPMCT285V);
            goodBPM = meanBPMH(i,isGoodPulse);
            goodBPMCT285V = meanBPMCT285V(isGoodPulse);
            tmpCorr = corrcoef(goodBPMCT285V,goodBPM);
            corrBPMCT285V_BPMH(i) = tmpCorr(1,2);

            % CT285V - BPMS
            isGoodPulse = ~(isNaNBPMS | isNaNBPMCT285V);
            goodBPM = meanBPMS(i,isGoodPulse);
            goodBPMCT285V = meanBPMCT285V(isGoodPulse);
            tmpCorr = corrcoef(goodBPMCT285V,goodBPM);
            corrBPMCT285V_BPMS(i) = tmpCorr(1,2);

            % CT285V - BPMV
            isGoodPulse = ~(isNaNBPMV | isNaNBPMCT285V);
            goodBPM = meanBPMV(i,isGoodPulse);
            goodBPMCT285V = meanBPMCT285V(isGoodPulse);
            tmpCorr = corrcoef(goodBPMCT285V,goodBPM);
            corrBPMCT285V_BPMV(i) = tmpCorr(1,2);

            % CT285S - BPMV
            isGoodPulse = ~(isNaNBPMV | isNaNBPMCT285S);
            goodBPM = meanBPMV(i,isGoodPulse);
            goodBPMCT285S = meanBPMCT285S(isGoodPulse);
            tmpCorr = corrcoef(goodBPMCT285S,goodBPM);
            corrBPMCT285S_BPMV(i) = tmpCorr(1,2);
        end
        
    end
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (useCL502)
        % CL502H - BPMH
        isGoodPulse = ~(isNaNBPMH | isNaNBPMCL502H);
        goodBPM = meanBPMH(i,isGoodPulse);
        goodBPMCL502H = meanBPMCL502H(isGoodPulse);
        tmpCorr = corrcoef(goodBPMCL502H,goodBPM);
        corrBPMCL502H_BPMH(i) = tmpCorr(1,2);

        % CL502H - BPMS
        isGoodPulse = ~(isNaNBPMS | isNaNBPMCL502H);
        goodBPM = meanBPMS(i,isGoodPulse);
        goodBPMCL502H = meanBPMCL502H(isGoodPulse);
        tmpCorr = corrcoef(goodBPMCL502H,goodBPM);
        corrBPMCL502H_BPMS(i) = tmpCorr(1,2);

        % CL502S - BPMH
        isGoodPulse = ~(isNaNBPMH | isNaNBPMCL502S);
        goodBPM = meanBPMH(i,isGoodPulse);
        goodBPMCL502S = meanBPMCL502S(isGoodPulse);
        tmpCorr = corrcoef(goodBPMCL502S,goodBPM);
        corrBPMCL502S_BPMH(i) = tmpCorr(1,2);

        % CL502S - BPMS
        isGoodPulse = ~(isNaNBPMS | isNaNBPMCL502S);
        goodBPM = meanBPMS(i,isGoodPulse);
        goodBPMCL502S = meanBPMCL502S(isGoodPulse);
        tmpCorr = corrcoef(goodBPMCL502S,goodBPM);
        corrBPMCL502S_BPMS(i) = tmpCorr(1,2);

        if (useVertical)
            % CL502H - BPMV
            isGoodPulse = ~(isNaNBPMV | isNaNBPMCL502H);
            goodBPM = meanBPMV(i,isGoodPulse);
            goodBPMCL502H = meanBPMCL502H(isGoodPulse);
            tmpCorr = corrcoef(goodBPMCL502H,goodBPM);
            corrBPMCL502H_BPMV(i) = tmpCorr(1,2);

            % CL502V - BPMH
            isGoodPulse = ~(isNaNBPMH | isNaNBPMCL502V);
            goodBPM = meanBPMH(i,isGoodPulse);
            goodBPMCL502V = meanBPMCL502V(isGoodPulse);
            tmpCorr = corrcoef(goodBPMCL502V,goodBPM);
            corrBPMCL502V_BPMH(i) = tmpCorr(1,2);

            % CL502V - BPMS
            isGoodPulse = ~(isNaNBPMS | isNaNBPMCL502V);
            goodBPM = meanBPMS(i,isGoodPulse);
            goodBPMCL502V = meanBPMCL502V(isGoodPulse);
            tmpCorr = corrcoef(goodBPMCL502V,goodBPM);
            corrBPMCL502V_BPMS(i) = tmpCorr(1,2);

            % CL502V - BPMV
            isGoodPulse = ~(isNaNBPMV | isNaNBPMCL502V);
            goodBPM = meanBPMV(i,isGoodPulse);
            goodBPMCL502V = meanBPMCL502V(isGoodPulse);
            tmpCorr = corrcoef(goodBPMCL502V,goodBPM);
            corrBPMCL502V_BPMV(i) = tmpCorr(1,2);

            % CL502S - BPMV
            isGoodPulse = ~(isNaNBPMV | isNaNBPMCL502S);
            goodBPM = meanBPMV(i,isGoodPulse);
            goodBPMCL502S = meanBPMCL502S(isGoodPulse);
            tmpCorr = corrcoef(goodBPMCL502S,goodBPM);
            corrBPMCL502S_BPMV(i) = tmpCorr(1,2);
        end
         
    end
    
end

%% Correlation between phase correlation and bpm correlation

if (useCT285)
    
    bpmH_CORR_Mon3_CT285H = corrcoef(corrMon3_BPMH,corrBPMCT285H_BPMH);
    bpmS_CORR_Mon3_CT285H = corrcoef(corrMon3_BPMS,corrBPMCT285H_BPMS);
    bpmH_CORR_Mon3_CT285S = corrcoef(corrMon3_BPMH,corrBPMCT285S_BPMH);
    bpmS_CORR_Mon3_CT285S = corrcoef(corrMon3_BPMS,corrBPMCT285S_BPMS);
    bpmH_CORR_Mon2_CT285H = corrcoef(corrMon2_BPMH,corrBPMCT285H_BPMH);
    bpmS_CORR_Mon2_CT285H = corrcoef(corrMon2_BPMS,corrBPMCT285H_BPMS);
    bpmH_CORR_Mon2_CT285S = corrcoef(corrMon2_BPMH,corrBPMCT285S_BPMH);
    bpmS_CORR_Mon2_CT285S = corrcoef(corrMon2_BPMS,corrBPMCT285S_BPMS);

    bpmH_CORR_Mon3_CT285H = bpmH_CORR_Mon3_CT285H(1,2);
    bpmS_CORR_Mon3_CT285H = bpmS_CORR_Mon3_CT285H(1,2);
    bpmH_CORR_Mon3_CT285S = bpmH_CORR_Mon3_CT285S(1,2);
    bpmS_CORR_Mon3_CT285S = bpmS_CORR_Mon3_CT285S(1,2);
    bpmH_CORR_Mon2_CT285H = bpmH_CORR_Mon2_CT285H(1,2);
    bpmS_CORR_Mon2_CT285H = bpmS_CORR_Mon2_CT285H(1,2);
    bpmH_CORR_Mon2_CT285S = bpmH_CORR_Mon2_CT285S(1,2);
    bpmS_CORR_Mon2_CT285S = bpmS_CORR_Mon2_CT285S(1,2);
    
    if (useVertical)
        bpmH_CORR_Mon3_CT285V = corrcoef(corrMon3_BPMH,corrBPMCT285V_BPMH);
        bpmS_CORR_Mon3_CT285V = corrcoef(corrMon3_BPMS,corrBPMCT285V_BPMS);
        bpmV_CORR_Mon3_CT285V = corrcoef(corrMon3_BPMV,corrBPMCT285V_BPMV);
        bpmV_CORR_Mon3_CT285H = corrcoef(corrMon3_BPMV,corrBPMCT285H_BPMV);
        bpmV_CORR_Mon3_CT285S = corrcoef(corrMon3_BPMV,corrBPMCT285S_BPMV);
        bpmV_CORR_Mon2_CT285H = corrcoef(corrMon2_BPMV,corrBPMCT285H_BPMV);
        bpmH_CORR_Mon2_CT285V = corrcoef(corrMon2_BPMH,corrBPMCT285V_BPMH);
        bpmS_CORR_Mon2_CT285V = corrcoef(corrMon2_BPMS,corrBPMCT285V_BPMS);
        bpmV_CORR_Mon2_CT285V = corrcoef(corrMon2_BPMV,corrBPMCT285V_BPMV);
        bpmV_CORR_Mon2_CT285S = corrcoef(corrMon2_BPMV,corrBPMCT285S_BPMV);

        bpmH_CORR_Mon3_CT285V = bpmH_CORR_Mon3_CT285V(1,2);
        bpmS_CORR_Mon3_CT285V = bpmS_CORR_Mon3_CT285V(1,2);
        bpmV_CORR_Mon3_CT285V = bpmV_CORR_Mon3_CT285V(1,2);
        bpmV_CORR_Mon3_CT285H = bpmV_CORR_Mon3_CT285H(1,2);
        bpmV_CORR_Mon3_CT285S = bpmV_CORR_Mon3_CT285S(1,2);    
        bpmV_CORR_Mon2_CT285H = bpmV_CORR_Mon2_CT285H(1,2);
        bpmH_CORR_Mon2_CT285V = bpmH_CORR_Mon2_CT285V(1,2);
        bpmS_CORR_Mon2_CT285V = bpmS_CORR_Mon2_CT285V(1,2);
        bpmV_CORR_Mon2_CT285V = bpmV_CORR_Mon2_CT285V(1,2);
        bpmV_CORR_Mon2_CT285S = bpmV_CORR_Mon2_CT285S(1,2);
    end
    
end

if (useCL502)
    
    bpmH_CORR_Mon3_CL502H = corrcoef(corrMon3_BPMH,corrBPMCL502H_BPMH);
    bpmS_CORR_Mon3_CL502H = corrcoef(corrMon3_BPMS,corrBPMCL502H_BPMS);
    bpmH_CORR_Mon3_CL502S = corrcoef(corrMon3_BPMH,corrBPMCL502S_BPMH);
    bpmS_CORR_Mon3_CL502S = corrcoef(corrMon3_BPMS,corrBPMCL502S_BPMS);
    bpmH_CORR_Mon3_CL502H = bpmH_CORR_Mon3_CL502H(1,2);
    bpmS_CORR_Mon3_CL502H = bpmS_CORR_Mon3_CL502H(1,2);
    bpmH_CORR_Mon3_CL502S = bpmH_CORR_Mon3_CL502S(1,2);
    bpmS_CORR_Mon3_CL502S = bpmS_CORR_Mon3_CL502S(1,2);

    bpmH_CORR_Mon2_CL502H = corrcoef(corrMon2_BPMH,corrBPMCL502H_BPMH);
    bpmS_CORR_Mon2_CL502H = corrcoef(corrMon2_BPMS,corrBPMCL502H_BPMS);
    bpmH_CORR_Mon2_CL502S = corrcoef(corrMon2_BPMH,corrBPMCL502S_BPMH);
    bpmS_CORR_Mon2_CL502S = corrcoef(corrMon2_BPMS,corrBPMCL502S_BPMS);
    bpmH_CORR_Mon2_CL502H = bpmH_CORR_Mon2_CL502H(1,2);
    bpmS_CORR_Mon2_CL502H = bpmS_CORR_Mon2_CL502H(1,2);
    bpmH_CORR_Mon2_CL502S = bpmH_CORR_Mon2_CL502S(1,2);
    bpmS_CORR_Mon2_CL502S = bpmS_CORR_Mon2_CL502S(1,2);

    if (useVertical)
        bpmV_CORR_Mon3_CL502H = corrcoef(corrMon3_BPMV,corrBPMCL502H_BPMV);
        bpmH_CORR_Mon3_CL502V = corrcoef(corrMon3_BPMH,corrBPMCL502V_BPMH);
        bpmS_CORR_Mon3_CL502V = corrcoef(corrMon3_BPMS,corrBPMCL502V_BPMS);
        bpmV_CORR_Mon3_CL502V = corrcoef(corrMon3_BPMV,corrBPMCL502V_BPMV);
        bpmV_CORR_Mon3_CL502S = corrcoef(corrMon3_BPMV,corrBPMCL502S_BPMV);
        bpmV_CORR_Mon3_CL502H = bpmV_CORR_Mon3_CL502H(1,2);
        bpmH_CORR_Mon3_CL502V = bpmH_CORR_Mon3_CL502V(1,2);
        bpmS_CORR_Mon3_CL502V = bpmS_CORR_Mon3_CL502V(1,2);
        bpmV_CORR_Mon3_CL502V = bpmV_CORR_Mon3_CL502V(1,2);
        bpmV_CORR_Mon3_CL502S = bpmV_CORR_Mon3_CL502S(1,2);

        bpmV_CORR_Mon2_CL502H = corrcoef(corrMon2_BPMV,corrBPMCL502H_BPMV);
        bpmH_CORR_Mon2_CL502V = corrcoef(corrMon2_BPMH,corrBPMCL502V_BPMH);
        bpmS_CORR_Mon2_CL502V = corrcoef(corrMon2_BPMS,corrBPMCL502V_BPMS);
        bpmV_CORR_Mon2_CL502V = corrcoef(corrMon2_BPMV,corrBPMCL502V_BPMV);
        bpmV_CORR_Mon2_CL502S = corrcoef(corrMon2_BPMV,corrBPMCL502S_BPMV);
        bpmV_CORR_Mon2_CL502H = bpmV_CORR_Mon2_CL502H(1,2);
        bpmH_CORR_Mon2_CL502V = bpmH_CORR_Mon2_CL502V(1,2);
        bpmS_CORR_Mon2_CL502V = bpmS_CORR_Mon2_CL502V(1,2);
        bpmV_CORR_Mon2_CL502V = bpmV_CORR_Mon2_CL502V(1,2);
        bpmV_CORR_Mon2_CL502S = bpmV_CORR_Mon2_CL502S(1,2);  
    end
    
end

%% plots

%%%%%%%%%% Scatter plots for max correlation (H,S,V) acros BPMs %%%%%%%%%%
% Mon3 - BPMH
[~,i] = max(abs(corrMon3_BPMH));
isNaNBPM = isnan(meanBPMH(i,:));
isGoodPulse = ~(isNaNBPM | isNaN3');
goodBPM = meanBPMH(i,isGoodPulse);
goodMon3 = meanMon3Phase(isGoodPulse);
tmpCorr = corrcoef(goodMon3,goodBPM);
corrMon3_BPMH(i) = tmpCorr(1,2);
figure;
scatter(goodMon3,goodBPM);
title(sprintf('Mon3-%sH: %.2f',bpmNames{i},corrMon3_BPMH(i)) ) 
xlabel('Mon3 Phase [degrees]');
ylabel('BPM Position [mm]');
if (savePlots)
    saveStr = sprintf('%scorr_Mon3_BPMH_%s',saveDir,bpmNames{i});
    print([saveStr '.png'],'-dpng');
    savefig([saveStr '.fig']);
end

% Mon3 - BPMV
if (useVertical)
    [~,i] = max(abs(corrMon3_BPMV));
    isNaNBPM = isnan(meanBPMV(i,:));
    isGoodPulse = ~(isNaNBPM | isNaN3');
    goodBPM = meanBPMV(i,isGoodPulse);
    goodMon3 = meanMon3Phase(isGoodPulse);
    tmpCorr = corrcoef(goodMon3,goodBPM);
    corrMon3_BPMV(i) = tmpCorr(1,2);
    figure;
    scatter(goodMon3,goodBPM);
    title(sprintf('Mon3-%sV: %.2f',bpmNames{i},corrMon3_BPMV(i)) )
    xlabel('Mon3 Phase [degrees]');
    ylabel('BPM Position [mm]');
    if (savePlots)
        saveStr = sprintf('%scorr_Mon3_BPMV_%s',saveDir,bpmNames{i});
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end
end

% Mon3 - BPMS
[~,i] = max(abs(corrMon3_BPMS));
isNaNBPM = isnan(meanBPMS(i,:));
isGoodPulse = ~(isNaNBPM | isNaN3');
goodBPM = meanBPMS(i,isGoodPulse);
goodMon3 = meanMon3Phase(isGoodPulse);
tmpCorr = corrcoef(goodMon3,goodBPM);
corrMon3_BPMS(i) = tmpCorr(1,2);
figure;
scatter(goodMon3,goodBPM);
title(sprintf('Mon3-%sS: %.2f',bpmNames{i},corrMon3_BPMS(i)) )
xlabel('Mon3 Phase [degrees]');
ylabel('BPM Transmission [V]');
if (savePlots)
    saveStr = sprintf('%scorr_Mon3_BPMS_%s',saveDir,bpmNames{i});
    print([saveStr '.png'],'-dpng');
    savefig([saveStr '.fig']);
end

% Mon2 - BPMH
[~,i] = max(abs(corrMon2_BPMH));
isNaNBPM = isnan(meanBPMH(i,:));
isGoodPulse = ~(isNaNBPM | isNaN2');
goodBPM = meanBPMH(i,isGoodPulse);
goodMon2 = meanMon2Phase(isGoodPulse);
tmpCorr = corrcoef(goodMon2,goodBPM);
corrMon2_BPMH(i) = tmpCorr(1,2);
figure;
scatter(goodMon2,goodBPM);
title(sprintf('Mon2-%sH: %.2f',bpmNames{i},corrMon2_BPMH(i)) )
xlabel('Mon2 Phase [degrees]');
ylabel('BPM Position [mm]');
if (savePlots)
    saveStr = sprintf('%scorr_Mon2_BPMH_%s',saveDir,bpmNames{i});
    print([saveStr '.png'],'-dpng');
    savefig([saveStr '.fig']);
end

% Mon2 - BPMV
if (useVertical)
    [~,i] = max(abs(corrMon2_BPMV));
    isNaNBPM = isnan(meanBPMV(i,:));
    isGoodPulse = ~(isNaNBPM | isNaN2');
    goodBPM = meanBPMV(i,isGoodPulse);
    goodMon2 = meanMon2Phase(isGoodPulse);
    tmpCorr = corrcoef(goodMon2,goodBPM);
    corrMon2_BPMV(i) = tmpCorr(1,2);
    figure;
    scatter(goodMon2,goodBPM);
    title(sprintf('Mon2-%sV: %.2f',bpmNames{i},corrMon2_BPMV(i)) )
    xlabel('Mon2 Phase [degrees]');
    ylabel('BPM Position [mm]');
    if (savePlots)
        saveStr = sprintf('%scorr_Mon2_BPMV_%s',saveDir,bpmNames{i});
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end
end

% Mon2 - BPMS
[~,i] = max(abs(corrMon2_BPMS));
isNaNBPM = isnan(meanBPMS(i,:));
isGoodPulse = ~(isNaNBPM | isNaN2');
goodBPM = meanBPMS(i,isGoodPulse);
goodMon2 = meanMon2Phase(isGoodPulse);
tmpCorr = corrcoef(goodMon2,goodBPM);
corrMon2_BPMS(i) = tmpCorr(1,2);
figure;
scatter(goodMon2,goodBPM);
title(sprintf('Mon2-%sS: %.2f',bpmNames{i},corrMon2_BPMS(i)) )
xlabel('Mon2 Phase [degrees]');
ylabel('BPM Transmission [V]');
if (savePlots)
    saveStr = sprintf('%scorr_Mon2_BPMS_%s',saveDir,bpmNames{i});
    print([saveStr '.png'],'-dpng');
    savefig([saveStr '.fig']);
end

%%%%%%%%%% Correlation vs. BPM %%%%%%%%%%

yMin = -1;
yMax = 1;

if (useVertical)
    figure;
    plot(corrMon2_BPMV,'b','LineWidth',2);
    hold all;
    plot(corrMon3_BPMV,'r','LineWidth',2);
    title('VERTICAL POSITION: Correlation with Phase')
    legend('Mon2 (CT)','Mon3 (CB)');
    correlationAnalysis_plotFormat_corrVsBPM;
    if (savePlots)
        saveStr = sprintf('%scorrBPMV_Phase',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end
end

figure;
plot(corrMon2_BPMH,'b','LineWidth',2);
hold all;
plot(corrMon3_BPMH,'r','LineWidth',2);
title('HORIZONTAL POSITION: Correlation with Phase')
legend('Mon2 (CT)','Mon3 (CB)');xlabel('BPM Index')
correlationAnalysis_plotFormat_corrVsBPM;
if (savePlots)
    saveStr = sprintf('%scorrBPMH_Phase',saveDir);
    print([saveStr '.png'],'-dpng');
    savefig([saveStr '.fig']);
end

figure;
plot(corrMon2_BPMS,'b','LineWidth',2);
hold all;
plot(corrMon3_BPMS,'r','LineWidth',2);
title('BPM TRANSMISSION: Correlation with Phase')
legend('Mon2 (CT)','Mon3 (CB)');
correlationAnalysis_plotFormat_corrVsBPM;
if (savePlots)
    saveStr = sprintf('%scorrBPMS_Phase',saveDir);
    print([saveStr '.png'],'-dpng');
    savefig([saveStr '.fig']);
end

%%%%%% Correlations: BPM CL502/CT285 signals and VERTICAL POSITION %%%%%%
if (useVertical)
    green = [0, 0.8, 0];
    yMin = -1;
    yMax = 1;

    if (useCL502)    
        figure;
        plot(corrMon2_BPMV,'b','LineWidth',2);
        hold all;
        plot(corrMon3_BPMV,'r','LineWidth',2);
        plot(corrBPMCL502H_BPMV, 'Color',green,'LineWidth',2)
        plot([bpmCL502Index bpmCL502Index],[yMin yMax],'--','Color',green,'LineWidth',2);
        title('Correlation with Vertical Position in all BPMs')
        legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CL502H');
        correlationAnalysis_plotFormat_corrVsBPM;
        if (savePlots)
            saveStr = sprintf('%scorrBPMV_CL502H',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end
        
        figure;
        plot(corrMon2_BPMV,'b','LineWidth',2);
        hold all;
        plot(corrMon3_BPMV,'r','LineWidth',2);
        plot(corrBPMCL502S_BPMV, 'Color',green,'LineWidth',2)
        plot([bpmCL502Index bpmCL502Index],[yMin yMax],'--','Color',green,'LineWidth',2);
        title('Correlation with Vertical Position in all BPMs')
        legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CL502S');
        correlationAnalysis_plotFormat_corrVsBPM;
        if (savePlots)
            saveStr = sprintf('%scorrBPMV_CL502S',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        figure;
        plot(corrMon2_BPMV,'b','LineWidth',2);
        hold all;
        plot(corrMon3_BPMV,'r','LineWidth',2);
        plot(corrBPMCL502V_BPMV, 'Color',green,'LineWidth',2)
        plot([bpmCL502Index bpmCL502Index],[yMin yMax],'--','Color',green,'LineWidth',2);
        title('Correlation with Vertical Position in all BPMs')
        legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CL502V');
        correlationAnalysis_plotFormat_corrVsBPM;
        if (savePlots)
            saveStr = sprintf('%scorrBPMV_CL502V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end
    end
    
    if (useCT285)
        figure;
        plot(corrMon2_BPMV,'b','LineWidth',2);
        hold all;
        plot(corrMon3_BPMV,'r','LineWidth',2);
        plot(corrBPMCT285H_BPMV, 'Color',green,'LineWidth',2)
        plot([bpmCT285Index bpmCT285Index],[yMin yMax],'--','Color',green,'LineWidth',2);
        title('Correlation with Vertical Position in all BPMs')
        grid on;
        legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CT285H');
        correlationAnalysis_plotFormat_corrVsBPM;
        if (savePlots)
            saveStr = sprintf('%scorrBPMV_CT285H',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        figure;
        plot(corrMon2_BPMV,'b','LineWidth',2);
        hold all;
        plot(corrMon3_BPMV,'r','LineWidth',2);
        plot(corrBPMCT285S_BPMV, 'Color',green,'LineWidth',2)
        plot([bpmCT285Index bpmCT285Index],[yMin yMax],'--','Color',green,'LineWidth',2);
        title('Correlation with Vertical Position in all BPMs')
        legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CT285S');
        correlationAnalysis_plotFormat_corrVsBPM;
        if (savePlots)
            saveStr = sprintf('%scorrBPMV_CT285S',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end

        figure;
        plot(corrMon2_BPMV,'b','LineWidth',2);
        hold all;
        plot(corrMon3_BPMV,'r','LineWidth',2);
        plot(corrBPMCT285V_BPMV, 'Color',green,'LineWidth',2)
        plot([bpmCT285Index bpmCT285Index],[yMin yMax],'--','Color',green,'LineWidth',2);
        title('Correlation with Vertical Position in all BPMs')
        legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CT285V');
        correlationAnalysis_plotFormat_corrVsBPM;
        if (savePlots)
            saveStr = sprintf('%scorrBPMV_CT285V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end 
    end
end

%%%%%% Correlations: BPM CL502/CT285 signals and HORIZONTAL POSITION %%%%%%
green = [0, 0.8, 0];
yMin = -1;
yMax = 1;

if (useCL502)    
    figure;
    plot(corrMon2_BPMH,'b','LineWidth',2);
    hold all;
    plot(corrMon3_BPMH,'r','LineWidth',2);
    plot(corrBPMCL502H_BPMH, 'Color',green,'LineWidth',2)
    plot([bpmCL502Index bpmCL502Index],[yMin yMax],'--','Color',green,'LineWidth',2);
    title('Correlation with Horizontal Position in all BPMs')
    legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CL502H');
    correlationAnalysis_plotFormat_corrVsBPM;
    if (savePlots)
        saveStr = sprintf('%scorrBPMH_CL502H',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    figure;
    plot(corrMon2_BPMH,'b','LineWidth',2);
    hold all;
    plot(corrMon3_BPMH,'r','LineWidth',2);
    plot(corrBPMCL502S_BPMH, 'Color',green,'LineWidth',2)
    plot([bpmCL502Index bpmCL502Index],[yMin yMax],'--','Color',green,'LineWidth',2);
    title('Correlation with Horizontal Position in all BPMs')
    legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CL502S');
    correlationAnalysis_plotFormat_corrVsBPM;
    if (savePlots)
        saveStr = sprintf('%scorrBPMH_CL502S',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    if (useVertical)
        figure;
        plot(corrMon2_BPMH,'b','LineWidth',2);
        hold all;
        plot(corrMon3_BPMH,'r','LineWidth',2);
        plot(corrBPMCL502V_BPMH, 'Color',green,'LineWidth',2)
        plot([bpmCL502Index bpmCL502Index],[yMin yMax],'--','Color',green,'LineWidth',2);
        title('Correlation with Horizontal Position in all BPMs')
        legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CL502V');
        correlationAnalysis_plotFormat_corrVsBPM;
        if (savePlots)
            saveStr = sprintf('%scorrBPMH_CL502V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end
    end
end

if (useCT285)
    figure;
    plot(corrMon2_BPMH,'b','LineWidth',2);
    hold all;
    plot(corrMon3_BPMH,'r','LineWidth',2);
    plot(corrBPMCT285H_BPMH, 'Color',green,'LineWidth',2)
    plot([bpmCT285Index bpmCT285Index],[yMin yMax],'--','Color',green,'LineWidth',2);
    title('Correlation with Horizontal Position in all BPMs')
    legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CT285H');
    correlationAnalysis_plotFormat_corrVsBPM;
    if (savePlots)
        saveStr = sprintf('%scorrBPMH_CT285H',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    figure;
    plot(corrMon2_BPMH,'b','LineWidth',2);
    hold all;
    plot(corrMon3_BPMH,'r','LineWidth',2);
    plot(corrBPMCT285S_BPMH, 'Color',green,'LineWidth',2)
    plot([bpmCT285Index bpmCT285Index],[yMin yMax],'--','Color',green,'LineWidth',2);
    title('Correlation with Horizontal Position in all BPMs')
    legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CT285S');
    correlationAnalysis_plotFormat_corrVsBPM;
    if (savePlots)
        saveStr = sprintf('%scorrBPMH_CT285S',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    if (useVertical)
        figure;
        plot(corrMon2_BPMH,'b','LineWidth',2);
        hold all;
        plot(corrMon3_BPMH,'r','LineWidth',2);
        plot(corrBPMCT285V_BPMH, 'Color',green,'LineWidth',2)
        plot([bpmCT285Index bpmCT285Index],[yMin yMax],'--','Color',green,'LineWidth',2);
        title('Correlation with Horizontal Position in all BPMs')
        legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CT285V');
        correlationAnalysis_plotFormat_corrVsBPM;
        if (savePlots)
            saveStr = sprintf('%scorrBPMH_CT285V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end
    end
end

%%%%%% Correlations: BPM CL502/CT285 signals and TRANSMISSION %%%%%%
green = [0, 0.8, 0];
yMin = -1;
yMax = 1;

if (useCL502)
    figure;
    plot(corrMon2_BPMS,'b','LineWidth',2);
    hold all;
    plot(corrMon3_BPMS,'r','LineWidth',2);
    plot(corrBPMCL502H_BPMS, 'Color',green,'LineWidth',2)
    plot([bpmCL502Index bpmCL502Index],[yMin yMax],'--','Color',green,'LineWidth',2);
    title('Correlation with Transmission in all BPMs')
    legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CL502H');
    correlationAnalysis_plotFormat_corrVsBPM;
    if (savePlots)
        saveStr = sprintf('%scorrBPMS_CL502H',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    figure;
    plot(corrMon2_BPMS,'b','LineWidth',2);
    hold all;
    plot(corrMon3_BPMS,'r','LineWidth',2);
    plot(corrBPMCL502S_BPMS, 'Color',green,'LineWidth',2)
    plot([bpmCL502Index bpmCL502Index],[yMin yMax],'--','Color',green,'LineWidth',2);
    title('Correlation with Transmission in all BPMs')
    legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CL502S');
    correlationAnalysis_plotFormat_corrVsBPM;
    if (savePlots)
        saveStr = sprintf('%scorrBPMS_CL502S',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    if (useVertical)
        figure;
        plot(corrMon2_BPMS,'b','LineWidth',2);
        hold all;
        plot(corrMon3_BPMS,'r','LineWidth',2);
        plot(corrBPMCL502V_BPMS, 'Color',green,'LineWidth',2)
        plot([bpmCL502Index bpmCL502Index],[yMin yMax],'--','Color',green,'LineWidth',2);
        title('Correlation with Transmission in all BPMs')
        legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CL502V');
        correlationAnalysis_plotFormat_corrVsBPM;
        if (savePlots)
            saveStr = sprintf('%scorrBPMS_CL502V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end
    end
end

if (useCT285)
    figure;
    plot(corrMon2_BPMS,'b','LineWidth',2);
    hold all;
    plot(corrMon3_BPMS,'r','LineWidth',2);
    plot(corrBPMCT285H_BPMS, 'Color',green,'LineWidth',2)
    plot([bpmCT285Index bpmCT285Index],[yMin yMax],'--','Color',green,'LineWidth',2);
    title('Correlation with Transmission in all BPMs')
    legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CT285H');
    correlationAnalysis_plotFormat_corrVsBPM;
    if (savePlots)
        saveStr = sprintf('%scorrBPMS_CT285H',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    figure;
    plot(corrMon2_BPMS,'b','LineWidth',2);
    hold all;
    plot(corrMon3_BPMS,'r','LineWidth',2);
    plot(corrBPMCT285S_BPMS, 'Color',green,'LineWidth',2)
    plot([bpmCT285Index bpmCT285Index],[yMin yMax],'--','Color',green,'LineWidth',2);
    title('Correlation with Transmission in all BPMs')
    legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CT285S');
    correlationAnalysis_plotFormat_corrVsBPM;
    if (savePlots)
        saveStr = sprintf('%scorrBPMS_CT285S',saveDir);
        print([saveStr '.png'],'-dpng');
        savefig([saveStr '.fig']);
    end

    if (useVertical)
        figure;
        plot(corrMon2_BPMS,'b','LineWidth',2);
        hold all;
        plot(corrMon3_BPMS,'r','LineWidth',2);
        plot(corrBPMCT285V_BPMS, 'Color',green,'LineWidth',2)
        plot([bpmCT285Index bpmCT285Index],[yMin yMax],'--','Color',green,'LineWidth',2);
        title('Correlation with Transmission in all BPMs')
        legend('Phase Mon2 (CT)','Phase Mon3 (CB)', 'BPM CT285V');
        correlationAnalysis_plotFormat_corrVsBPM;
        if (savePlots)
            saveStr = sprintf('%scorrBPMS_CT285V',saveDir);
            print([saveStr '.png'],'-dpng');
            savefig([saveStr '.fig']);
        end
    end
end

%%% Correlation between phase correlation and bpm correlation
figure;
correlationAnalysis_plots_PhaseCorrVsBPMCorr;

%% save data
if (savePlots)
   save([saveDir 'correlationAnalysis.mat']); 
end