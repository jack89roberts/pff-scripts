close all;
% load processed data first
%%
% Crates that each name corresponds to:
% frascati = cfv-2001-ccaos5
% pets = cfi-2010-ccaos3
% bpm1 = dctfumal (CL, CT 155,215,242,258,285,335,430)
% bpm2 = dctfumad (CT 105,487,495,515,532)
% bpm3 = dctfumar (CT 608,622,645,665,692,722,758, CR, CC 130)
% bpm4 = cfc-2010-bpm (CC, CB)

% load file with device crate information
filename = 'devices/deviceCrates.dat';
delimiter = ',';
formatSpec = '%s%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
devNames = dataArray{:, 1};
crateNames = dataArray{:, 2};

bpmNames = processedData.bpmData.bpmNames;
nBPMs = length(bpmNames);

figAllBPM1 = figure();
figAllBPM2 = figure();
figAllBPM3 = figure();
figAllBPM4 = figure();
for i=1:nBPMs
    tmpBPMCrateIndex = find(ismember(devNames,bpmNames{i}));
    if (tmpBPMCrateIndex > 0)
        tmpCrateName = crateNames{tmpBPMCrateIndex};
    else
        error('Could not find %s in list of device crates.\n',bpmNames{i});
    end
    switch (tmpCrateName)
        case 'dctfumal'
            figure(figAllBPM1);
        case 'dctfumad'
            figure(figAllBPM2);
        case 'dctfumar'
            figure(figAllBPM3);
        case 'cfc-2010-bpm'
            figure(figAllBPM4);
        otherwise
            error('No case exists for BPM %s on crate %s\n',bpmNames{i},tmpCrateName);
    end
    
    plot(getAligned2(eval(sprintf('processedData.bpmData.%sS.Samples.samples',bpmNames{i})),20,300,1)');
    hold all;
end

figure(figAllBPM1)
title(sprintf('All BPMs, BPM Crate 1 (dctfumal)'));
figure(figAllBPM2)
title(sprintf('All BPMs, BPM Crate 2 (dctfumad)'));
figure(figAllBPM3)
title(sprintf('All BPMs, BPM Crate 3 (dctfumar)'));
figure(figAllBPM4)
title(sprintf('All BPMs, BPM Crate 4 (cfc-2010-bpm)'));
%%

% extract signals to use from each crate
bpm1S = processedData.bpmData.CT.SVBPM0430S.Samples.samples;
bpm2S = processedData.bpmData.CT.SVBPI0487S.Samples.samples;
bpm3S = processedData.bpmData.CT.SVBPI0608S.Samples.samples;
bpm4S = processedData.bpmData.CC.SVBPM0235S.Samples.samples;
frascatiD = squeeze(processedData.phaseData.frascatiDiodes(2,:,:));
petsI = processedData.phaseData.petsI;
% pki1A = extractCTFSignalFromMergedData('CK.SVPKI11A.Samples.samples.value',CTFData);
% pki2A = extractCTFSignalFromMergedData('CK.SVPKI06A.Samples.samples.value',CTFData);
% pki3A = extractCTFSignalFromMergedData('CK.SVPKI15A.Samples.samples.value',CTFData);

bpm1H = processedData.bpmData.CT.SVBPM0430H.Samples.samples;
bpm2H = processedData.bpmData.CT.SVBPI0487H.Samples.samples;
bpm3H = processedData.bpmData.CT.SVBPI0608H.Samples.samples;
bpm4H = processedData.bpmData.CC.SVBPM0235H.Samples.samples;
frascatiM = squeeze(processedData.phaseData.frascatiMixers(2,:,:));
petsQ = processedData.phaseData.petsQ;
% pki1P = extractCTFSignalFromMergedData('CK.SVPKI11P.Samples.samples.value',CTFData);
% pki2P = extractCTFSignalFromMergedData('CK.SVPKI06P.Samples.samples.value',CTFData);
% pki3P = extractCTFSignalFromMergedData('CK.SVPKI15P.Samples.samples.value',CTFData);


% remove any duplicate pulses coming from bad acquisition
bpm1S = removeDuplicatePulses(bpm1S);
bpm2S = removeDuplicatePulses(bpm2S);
bpm3S = removeDuplicatePulses(bpm3S);
bpm4S = removeDuplicatePulses(bpm4S);
frascatiD = removeDuplicatePulses(frascatiD);
petsI = removeDuplicatePulses(petsI);
bpm1H = removeDuplicatePulses(bpm1H);
bpm2H = removeDuplicatePulses(bpm2H);
bpm3H = removeDuplicatePulses(bpm3H);
bpm4H = removeDuplicatePulses(bpm4H);
frascatiM = removeDuplicatePulses(frascatiM);
petsQ = removeDuplicatePulses(petsQ);

% align signals
[bpm1SA,bpm1HA] = getAligned2(bpm1S,20,300,1,{bpm1H});
[bpm2SA,bpm2HA] = getAligned2(bpm2S,20,300,1,{bpm2H});
[bpm3SA,bpm3HA] = getAligned2(bpm3S,20,300,1,{bpm3H});
[bpm4SA,bpm4HA] = getAligned2(bpm4S,20,300,1,{bpm4H});
[frascatiDA,frascatiMA] = getAligned2(frascatiD,20,300,1,{frascatiM});
[petsQA,petsIA] = getAligned2(petsQ,20,300,1,{petsI});
% [pki1AA,pki1PA] = getAligned2(pki1A,20,300,1,{pki1P});
% [pki2AA,pki2PA] = getAligned2(pki2A,20,300,1,{pki2P});
% [pki3AA,pki3PA] = getAligned2(pki3A,20,300,1,{pki3P});
bpm1HA = bpm1HA{1};
bpm2HA = bpm2HA{1};
bpm3HA = bpm3HA{1};
bpm4HA = bpm4HA{1};
frascatiMA = frascatiMA{1};
petsIA = petsIA{1};
% pki1PA = pki1PA{1};
% pki2PA = pki2PA{1};
% pki3PA = pki3PA{1};

figure;
subplot(1,2,1);
plot(bpm1S');
title('BPM1S Raw')
subplot(1,2,2);
plot(bpm1SA');
title('BPM1S Aligned')

figure;
subplot(1,2,1);
plot(bpm2S');
title('BPM2S Raw')
subplot(1,2,2);
plot(bpm2SA');
title('BPM2S Aligned')

figure;
subplot(1,2,1);
plot(bpm3S');
title('BPM3S Raw')
subplot(1,2,2);
plot(bpm3SA');
title('BPM3S Aligned')

figure;
subplot(1,2,1);
plot(bpm4S');
title('BPM4S Raw')
subplot(1,2,2);
plot(bpm4SA');
title('BPM4S Aligned')

figure;
subplot(1,2,1);
plot(frascatiD');
title('FrascatiD Raw')
subplot(1,2,2);
plot(frascatiDA');
title('FrascatiD Aligned')

figure;
subplot(1,2,1);
plot(petsQ');
title('PETSQ Raw')
subplot(1,2,2);
plot(petsQA');
title('PETSQ Aligned')

sample = 50;
figure;
plot(bpm1SA(:,sample)./max(abs(bpm1SA(:,sample))));
hold all;
plot(bpm2SA(:,sample)./max(abs(bpm2SA(:,sample))));
plot(bpm3SA(:,sample*2)./max(abs(bpm3SA(:,sample*2))));
plot(bpm4SA(:,sample*2)./max(abs(bpm4SA(:,sample*2))));
plot(frascatiDA(:,round(sample*(500/384)))./max(abs(frascatiDA(:,round(sample*(500/384))))));
plot(petsQA(:,round(sample*(1000/384)))./max(abs(petsQA(:,round(sample*(1000/384))))));
% plot(pki1AA(:,round(sample*(1000/384)))./max(abs(pki1AA(:,round(sample*(1000/384))))));
% plot(pki2AA(:,round(sample*(1000/384)))./max(abs(pki2AA(:,round(sample*(1000/384))))));
% plot(pki3AA(:,round(sample*(1000/384)))./max(abs(pki3AA(:,round(sample*(1000/384))))));
legend('BPM1','BPM2','BPM3','BPM4','Frascati','PETS');
title('VALUE vs. TIME');

[xcfBPM2,lagsBPM2] = nancrosscorr(bpm1SA(:,sample),bpm2SA(:,sample));
[xcfBPM3,lagsBPM3] = nancrosscorr(bpm1SA(:,sample),bpm3SA(:,sample*2));
[xcfBPM4,lagsBPM4] = nancrosscorr(bpm1SA(:,sample),bpm4SA(:,sample*2));
[xcfFrascati,lagsFrascati] = nancrosscorr(bpm1SA(:,sample),frascatiDA(:,round(sample*(500/384))));
[xcfPETS,lagsPETS] = nancrosscorr(bpm1SA(:,sample),petsQA(:,round(sample*(1000/384))));
% [xcfPKI1,lagsPKI1,~] = crosscorr(bpm1SA(:,50),pki1AA(:,round(sample*(1000/384))));
% [xcfPKI2,lagsPKI2,~] = crosscorr(bpm1SA(:,50),pki2AA(:,round(sample*(1000/384))));
% [xcfPKI3,lagsPKI3,~] = crosscorr(bpm1SA(:,50),pki3AA(:,round(sample*(1000/384))));

[~,maxIndBPM2] = max(abs(xcfBPM2));
[~,maxIndBPM3] = max(abs(xcfBPM3));
[~,maxIndBPM4] = max(abs(xcfBPM4));
[~,maxIndFrascati] = max(abs(xcfFrascati));
[~,maxIndPETS] = max(abs(xcfPETS));
legStrBPM2 = sprintf('BPM2 (%.2f at %d)',xcfBPM2(maxIndBPM2),lagsBPM2(maxIndBPM2));
legStrBPM3 = sprintf('BPM3 (%.2f at %d)',xcfBPM3(maxIndBPM3),lagsBPM3(maxIndBPM3));
legStrBPM4 = sprintf('BPM4 (%.2f at %d)',xcfBPM4(maxIndBPM4),lagsBPM4(maxIndBPM4));
legStrFrascati = sprintf('Frascati (%.2f at %d)',xcfFrascati(maxIndFrascati),lagsFrascati(maxIndFrascati));
legStrPETS = sprintf('PETS (%.2f at %d)',xcfPETS(maxIndPETS),lagsPETS(maxIndPETS));

figure
plot(lagsBPM2,xcfBPM2)
hold all
plot(lagsBPM3,xcfBPM3)
plot(lagsBPM4,xcfBPM4)
plot(lagsFrascati,xcfFrascati)
plot(lagsPETS,xcfPETS)
% plot(lagsPKI1,xcfPKI1)
% plot(lagsPKI2,xcfPKI2)
% plot(lagsPKI3,xcfPKI3)
legend(legStrBPM2, legStrBPM3, legStrBPM4, legStrFrascati, legStrPETS);
title('CROSS CORRELATION (Relative to BPM1)')

% figure;
% tmp1 = frascatiDA;
% tmp1(end,:) = [];
% tmp2 = bpm1SA;
% tmp2(1,:) = [];
% plot(tmp1(:,round(sample*(500/384)))./max(abs(tmp1(:,round(sample*(500/384))))));
% hold all;
% plot(tmp2(:,sample)./max(abs(tmp2(:,sample))));


% 
% figure;
% subplot(1,2,1);
% plot(bpm1H','b');
% hold all;
% plot(bpm2H','r');
% plot(bpm3H','g');
% plot(bpm4H','m');
% plot(frascatiM','k');
% title('RAW');
% subplot(1,2,2);
% plot(bpm1HA','b');
% hold all;
% plot(bpm2HA','r');
% plot(bpm3HA','g');
% plot(bpm4HA','m');
% plot(frascatiMA','k');
% title('ALIGNED')
% 
% sample = 50;
% figure;
% plot(bpm1HA(:,sample)./max(abs(bpm1HA(:,sample))),'b');
% hold all;
% plot(bpm2HA(:,sample)./max(abs(bpm2HA(:,sample))),'r');
% plot(bpm3HA(:,sample*2)./max(abs(bpm3HA(:,sample*2))),'g');
% plot(bpm4HA(:,sample*2)./max(abs(bpm4HA(:,sample*2))),'m');
% plot(frascatiMA(:,round(sample*(500/384)))./max(abs(frascatiMA(:,round(sample*(500/384))))),'k');


% 
% 
% filename = '/home/jack/Documents/MATLAB/PhaseFeedforward/ctfmod/PhaseFeedforward/devices/deviceCrates.dat';
% delimiter = ',';
% formatSpec = '%s%s%[^\n\r]';
% fileID = fopen(filename,'r');
% dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
% fclose(fileID);
% allDevs = dataArray{:, 1};
% allCrateNames = dataArray{:, 2};
% clearvars filename delimiter formatSpec fileID dataArray ans;
% 
% uniqueCrates = unique(allCrateNames);
% nDevs = length(allDevs);
% allCrateIndices = NaN*ones(1,nDevs);
% 
% for i=1:nDevs
%     [~,crateIndex] = ismember(allCrateNames{i},uniqueCrates);
%     allCrateIndices(i) = crateIndex;
% end