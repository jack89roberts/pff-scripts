%% inputs
dataSetNames = {...
'20151120_1526_Gain-950_R56_0.1_Interleaaved',...
'20151120_1538_Gain-800_R56_0.1_Interleaaved',...
'20151120_1551_Gain-750_R56_0.1_Interleaved',...
'20151120_1600_Gain950_R56_0.1_Interleaaved',...
'20151120_1604_Gain-775_R56_0.1_Interleaved',...
'20151120_1611_Gain-775_R56_0.1_Interleaved',...
'20151120_1617_Gain-825_R56_0.1_Interleaved',...
% '20151120_1622_Gain-825_R56_0.1_Interleaved',...
% '20151120_1639_Gain-790_R56_0.1_Interleaved',...
% '20151120_1644_Gain-790_R56_0.1_Interleaved',...
% '20151120_1650_Gain-810_R56_0.1_Interleaved',...
% '20151120_1654_Gain-810_R56_0.1_Interleaved',...
% '20151120_1700_Gain-850_R56_0.1_Interleaved',...
% '20151120_1721_Gain-850_R56_0.1_Interleaved',...
% '20151120_1800_Gain-850_R56_0.1_Interleaved'...
};

%% definition of klystron signal names
pkiANames = {...
    'CK_SVPKI02A'
    'CK_SVPKI03A'
    'CK_SVPKI05A'
    'CK_SVPKI06A'
    'CK_SVPKI07A'
    'CK_SVPKI11A'
    'CK_SVPKI12A'
    'CK_SVPKI13A'
    'CK_SVPKI14A'
    'CK_SVPKI15A'
};

pkiPNames = {...
    'CK_SVPKI02P'
    'CK_SVPKI03P'
    'CK_SVPKI05P'
    'CK_SVPKI06P'
    'CK_SVPKI07P'
    'CK_SVPKI11P'
    'CK_SVPKI12P'
    'CK_SVPKI13P'
    'CK_SVPKI14P'
    'CK_SVPKI15P'
};

psiANames = {...
    'CK_SVPSI03A'
    'CK_SVPSI05A'
    'CK_SVPSI06A'
    'CK_SVPSI07A'
    'CK_SVPSI11A'
    'CK_SVPSI12A'
    'CK_SVPSI13A'
    'CK_SVPSI14A'
    'CK_SVPSI15A'
};

psiPNames = {...
    'CK_SVPSI03P'
    'CK_SVPSI05P'
    'CK_SVPSI06P'
    'CK_SVPSI07P'
    'CK_SVPSI11P'
    'CK_SVPSI12P'
    'CK_SVPSI13P'
    'CK_SVPSI14P'
    'CK_SVPSI15P'
};

peiNames = {...
    'CK_SVPEI0305A'
    'CK_SVPEI0530A'
    'CK_SVPEI0630A'
    'CK_SVPEI0730A'
    'CK_SVPEI1130A'
    'CK_SVPEI1230A'
    'CK_SVPEI1330A'
    'CK_SVPEI1430A'
    'CK_SVPEI1530A'
};

sigToExtract = '.Samples.samples.value';

%% load data
nDataSets = length(dataSetNames);
CTFData = [];
for i=1:nDataSets
    CTFData = [CTFData, loadMergedData(dataSetNames{i})];
end
nPulses = length(CTFData);

nPKIA = length(pkiANames);
nPKIP = length(pkiPNames);
nPSIA = length(psiANames);
nPSIP = length(psiPNames);
nPEI = length(peiNames);
nSamplesKlystron = length(extractCTFSignalFromMergedData([pkiANames{1} sigToExtract],CTFData(1)));

rawPKIA = NaN(nPKIA,nPulses,nSamplesKlystron);
rawPKIP = NaN(nPKIP,nPulses,nSamplesKlystron);
rawPSIA = NaN(nPSIA,nPulses,nSamplesKlystron);
rawPSIP = NaN(nPSIP,nPulses,nSamplesKlystron);
rawPEI = NaN(nPEI,nPulses,nSamplesKlystron);

for i=1:nPKIA
    rawPKIA(i,:,:) = extractCTFSignalFromMergedData([pkiANames{i} sigToExtract],CTFData);
end
for i=1:nPKIP
    rawPKIP(i,:,:) = extractCTFSignalFromMergedData([pkiPNames{i} sigToExtract],CTFData);
end
for i=1:nPSIA
    rawPSIA(i,:,:) = extractCTFSignalFromMergedData([psiANames{i} sigToExtract],CTFData);
end
for i=1:nPSIP
    rawPSIP(i,:,:) = extractCTFSignalFromMergedData([psiPNames{i} sigToExtract],CTFData);
end
for i=1:nPEI
    rawPEI(i,:,:) = extractCTFSignalFromMergedData([peiNames{i} sigToExtract],CTFData);
end
%% choose sample range
klysSampleRange = 630;
% figure;
% for i=1:nPKIA
%     plot(squeeze(rawPKIA(i,1,:)))
%     hold all;
% end
% figure;
% for i=1:nPKIP
%     plot(squeeze(rawPKIP(i,1,:)))
%     hold all
% end
% figure;
% for i=1:nPSIA
%     plot(squeeze(rawPSIA(i,1,:)))
%     hold all;
% end
% figure;
% for i=1:nPSIP
%     plot(squeeze(rawPSIP(i,1,:)))
%     hold all;
% end
% figure;
% for i=1:nPEI
%     plot(squeeze(rawPEI(i,1,:)))
%     hold all;
% end

%% calculate means
meanPKIA = nanmean(rawPKIA(:,:,klysSampleRange),3);
meanPKIP = nanmean(rawPKIP(:,:,klysSampleRange),3);
meanPSIA = nanmean(rawPSIA(:,:,klysSampleRange),3);
meanPSIP = nanmean(rawPSIP(:,:,klysSampleRange),3);
meanPEI = nanmean(rawPEI(:,:,klysSampleRange),3);

figure;
plot(meanPKIA')
legend(pkiANames)

figure;
plot(meanPKIP')
legend(pkiPNames)

figure;
plot(meanPSIA')
legend(psiANames)

figure;
plot(meanPSIP')
legend(psiPNames)

figure;
plot(meanPEI')
legend(peiNames)