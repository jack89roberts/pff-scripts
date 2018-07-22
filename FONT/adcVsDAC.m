% To analyse data where constant DAC output sent to ADCs to check their
% linearity.
clearvars; close all;

%% INPUTS

dataDir = '/home/jack/PhaseFeedforward/FONTData/201510/20151026_1340_ADCCheck/Extracted';
ADCsToUse = [2, 4]; % ADCs where DAC out applied
dacValueString = 'DAC'; % substring in file name used to identify DAC out
sampleRange = 200:800;

%%

[rawADCs, dacValues] = packageFONTDataByString(dataDir,dacValueString);

ADCs = rawADCs(ADCsToUse,:,:,sampleRange);
meanADCs = nanmean(nanmean(ADCs,4),3);

figure;
plot(dacValues,meanADCs');