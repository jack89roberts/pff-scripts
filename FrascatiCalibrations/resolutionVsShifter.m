% load('/home/jack/PhaseFeedforward/CTFData/201511/FrascatiCalibrations/frascatiCalibration_20151116_1334.mat');
% [ calibrationConstants, useMixerOverSqrtDiode ] = loadFrascatiCalibrationConstants('/home/jack/PhaseFeedforward/CTFData/201511/FrascatiCalibrations/frascatiCalibrationConstants_20151116_1334');
%sampRange = 560:660;

load('/home/jack/PhaseFeedforward/CTFData/201605/FrascatiCalibrations/frascatiCalibration_20160525_1514.mat');
[ calibrationConstants, useMixerOverSqrtDiode ] = loadFrascatiCalibrationConstants('/home/jack/PhaseFeedforward/CTFData/201605/FrascatiCalibrations/frascatiCalibrationConstants_20160525_1514');
sampRange = 660:710;
%%
[nMons,nPoints,nPulses,nSamples]=size(rawMixers);
[mixers, diodes, sampRange] = frascatiCalProcessSignals(rawMixers,rawDiodes, rawMixerSensitivities,rawDiodeSensitivities,sampRange);
%%
phases = NaN(nMons,nPoints,nPulses,nSamples);
for m=1:nMons
    for p=1:nPoints
        phases(m,p,:,:) = getPhaseMixerDiode(squeeze(mixers(m,p,:,:)),squeeze(diodes(m,p,:,:)),calibrationConstants(m,1),calibrationConstants(m,4),useMixerOverSqrtDiode);
    end
end

meanPhases = nanmean(phases(:,:,:,sampRange),4);

[pointMean,stdMean,pointMean_err,stdMean_err] = nanMeanStdErr(meanPhases,3);

diff12 = squeeze(meanPhases(2,:,:)-meanPhases(1,:,:));
[~,res12,~,res12_err] = nanMeanStdErr(diff12,2);
res12 = res12./sqrt(2);
res12_err = res12_err./sqrt(2);

figure;
errorbar((scanPhShiftValues.*3)-80,res12,res12_err,'bo-','LineWidth',2,'MarkerFaceColor','b');

diff13 = squeeze(meanPhases(3,:,:)-meanPhases(1,:,:));
[~,res13,~,res13_err] = nanMeanStdErr(diff13,2);
res13 = res13./sqrt(2);
res13_err = res13_err./sqrt(2);

figure;
fudge = 0.05;
errorbar((scanPhShiftValues.*3)-216,res13+fudge,res13_err,'ro-','LineWidth',2,'MarkerFaceColor','r');
xlabel('Beam--LO Phase Offset [degrees]')
ylabel('Resolution [degrees]')
title('Resolution vs. Phase Offset')
format_plots
