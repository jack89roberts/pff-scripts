close all;
%clear all;
%%

dataDir = '/home/jack/Documents/CTFData/PFF/201406/20140618_2125';

frascatiSampleRange = 200:300;
frascatiCalibrationConstant = 1.4212;

petsSampleRange = 200:300;

tblSampleRange = 140:220;
tl2SampleRange = 320:400;
%%

CTFData = mergeMatMonData(dataDir);
nPulses = length(CTFData);

frascatiDiodeString = 'CT.SCOPE01.CH01.Acquisition.value.value';
frascatiDiode = extractCTFSignalFromMergedData(frascatiDiodeString, CTFData);
frascatiDiodeSensitivityString = 'CT.SCOPE01.CH01.Acquisition.sensitivity.value';
frascatiDiodeSensitivity = extractCTFSignalFromMergedData(frascatiDiodeSensitivityString, CTFData);
frascatiMixerString = 'CT.SCOPE01.CH02.Acquisition.value.value';
frascatiMixer = extractCTFSignalFromMergedData(frascatiMixerString, CTFData);
frascatiMixerSensitivityString = 'CT.SCOPE01.CH02.Acquisition.sensitivity.value';
frascatiMixerSensitivity = extractCTFSignalFromMergedData(frascatiMixerSensitivityString, CTFData);

petsIString = 'CE.SCOPE03.CH01.Acquisition.value.value';
petsI = extractCTFSignalFromMergedData(petsIString, CTFData);
petsQString = 'CE.SCOPE03.CH02.Acquisition.value.value';
petsQ = extractCTFSignalFromMergedData(petsQString, CTFData);

tblTransmissionString = 'CB.SVBPM0150S.Samples.samples.value';
tblTransmission = extractCTFSignalFromMergedData(tblTransmissionString, CTFData);
tl2TransmissionString = 'CC.SVBPM0235S.Samples.samples.value';
tl2Transmission = extractCTFSignalFromMergedData(tl2TransmissionString, CTFData);


%%

[~,nSamplesFrascati] = size(frascatiDiode);
[~,nSamplesPETS] = size(petsI);
phasePETS = zeros(nPulses,nSamplesPETS);
phaseFrascati = zeros(nPulses,nSamplesFrascati);

for t=1:nPulses
    tmpFrasDiode = frascatiDiode(t,:);
    tmpFrasMixer = frascatiMixer(t,:);
    tmpFrasDiodeSens = frascatiDiodeSensitivity(t);
    tmpFrasMixerSens = frascatiMixerSensitivity(t);
    phaseFrascati(t,:) = getPhaseMixerDiode(tmpFrasMixer,tmpFrasDiode,tmpFrasMixerSens,tmpFrasDiodeSens,frascatiCalibrationConstant);
    
    tmpPETSI = petsI(t,:);
    tmpPETSQ = petsQ(t,:);
    phasePETS(t,:) = getPhaseIQ(tmpPETSI,tmpPETSQ);
  
end

meanPhasePETS = mean(phasePETS(:,petsSampleRange),2);
meanPhasePETS = meanPhasePETS - nanmean(meanPhasePETS);

meanPhaseFrascati = mean(phaseFrascati(:,frascatiSampleRange),2);
meanPhaseFrascati = meanPhaseFrascati - nanmean(meanPhaseFrascati);

stdMeanFrascati = std(meanPhaseFrascati)
stdMeanPETS = std(meanPhasePETS)

figure;
plot(meanPhaseFrascati,'b');
hold all
plot(meanPhasePETS,'r');
xlabel('Time [Pulse No.]');
ylabel('Mean Phase [degrees]');
legend('Frascati','PETS');
print([dataDir '/plots/meanPhaseVsTimeFrascatiPETS.eps'],'-depsc');
print([dataDir '/plots/meanPhaseVsTimeFrascatiPETS.png'],'-dpng');
savefig([dataDir '/plots/meanPhaseVsTimeFrascatiPETS.fig']);

figure;
minFrascati = min(meanPhaseFrascati);
maxFrascati = max(meanPhaseFrascati);
minPETS = min(meanPhasePETS);
maxPETS = max(meanPhasePETS);
minPhase = min([minFrascati minPETS]);
maxPhase = max([maxFrascati maxPETS]);

[nbf,xbf] = hist(meanPhaseFrascati,linspace(minPhase,maxPhase,50));
[nbp,xbp] = hist(meanPhasePETS,linspace(minPhase,maxPhase,50));
hold all
bar(xbf,nbf,1,'FaceColor','b')
bar(xbp,nbp,1,'FaceColor','r')
xlabel('Mean Phase [degrees]');
ylabel('Bin Count');
frasLegStr = ['Frascati: \sigma = ' num2str(stdMeanFrascati,'%.2f') '\circ'];
petsLegStr = ['PETS: \sigma = ' num2str(stdMeanPETS,'%.2f') '\circ'];
legend(frasLegStr,petsLegStr)
print([dataDir '/plots/histMeanPhaseFrascatiPETS.eps'],'-depsc');
print([dataDir '/plots/histMeanPhaseFrascatiPETS.png'],'-dpng');
savefig([dataDir '/plots/histMeanPhaseFrascatiPETS.fig']);



%%

% fitting/correlation doesn't work with NaN values so remove them
notNaNMeanFrascati = ~isnan(meanPhaseFrascati);
notNaNMeanPETS = ~isnan(meanPhasePETS);
notNaNBoth = notNaNMeanFrascati & notNaNMeanPETS;
validMeanPhaseFrascati = meanPhaseFrascati(notNaNBoth);
validMeanPhasePETS = meanPhasePETS(notNaNBoth);

fitParams = polyfit(validMeanPhaseFrascati,validMeanPhasePETS,1);
gradient = fitParams(1);
correlation = corrcoef(validMeanPhaseFrascati,validMeanPhasePETS);
correlation = correlation(1,2);

xFit = linspace(min(meanPhaseFrascati),max(meanPhaseFrascati),100);
yFit = gradient.*xFit + fitParams(2);

figure;
scatter(meanPhaseFrascati,meanPhasePETS);
hold all;
plot(xFit,yFit,'Color','k','LineWidth',2);

xlabel('Mean Phase Frascati [degrees]');
ylabel('Mean Phase PETS [degrees]');
datStr = sprintf('Data (Correlation = %.2f)',correlation);
fitStr = sprintf('Linear Fit (PETS = %.2f*Frascati + %.2f)', gradient, fitParams(2));
legend(datStr,fitStr);
print([dataDir '/plots/scatterFrascatiPETS.eps'],'-depsc');
print([dataDir '/plots/scatterFrascatiPETS.png'],'-dpng');
savefig([dataDir '/plots/scatterFrascatiPETS.fig']);


gradient
correlation

%%

meanTransmissionTL2 = nanmean(tl2Transmission(:,tl2SampleRange),2);
meanTransmissionTBL = nanmean(tblTransmission(:,tblSampleRange),2);

stdTransmissionTL2 = nanstd(tl2Transmission(:,tl2SampleRange),0,2);
stdTransmissionTBL = nanstd(tblTransmission(:,tblSampleRange),0,2);

figure;
plot(meanTransmissionTL2)
hold all
plot(meanTransmissionTBL);
xlabel('Time [Pulse No,]');
ylabel('Mean Transmission [A]');
legend('TL2 (CC.BPM0235)', 'TBL (CB.BPM0150)','Location','best')
print([dataDir '/plots/transmission.eps'],'-depsc');
print([dataDir '/plots/transmission.png'],'-dpng');
savefig([dataDir '/plots/transmission.fig']);

%% correlation in jitter (difference between two points)

diffFrascati = meanPhaseFrascati(2:length(meanPhaseFrascati))-meanPhaseFrascati(1:length(meanPhaseFrascati)-1);
diffPETS = meanPhasePETS(2:length(meanPhasePETS))-meanPhasePETS(1:length(meanPhasePETS)-1);
figure
plot(diffPETS,'g')
hold all
plot(diffFrascati,'b')
figure
scatter(diffFrascati,diffPETS)
hold off
scatter(diffFrascati,diffPETS)
corrcoef(diffFrascati,diffPETS)