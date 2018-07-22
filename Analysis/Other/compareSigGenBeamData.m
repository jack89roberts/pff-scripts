clearvars; close all;

sigGenData = load('/home/jack/PhaseFeedforward/Analysis/201508/SignalGenerator/processedData.mat');
beamData = load('/home/jack/PhaseFeedforward/Analysis/201507/DiodeVsAttenuation/processedData.mat');

monIndex = beamData.monIndex;

sigGenNormIndex = 9;
beamNormIndex = 2;

figure;
plot(sigGenData.powerLevels, sigGenData.sqrtMeanDiodes(monIndex,:)/max(sigGenData.sqrtMeanDiodes(monIndex,:)),'LineWidth',2);
hold all;
plot(beamData.inputPowerMaxMixdBm, beamData.meanSqrtDiodeMaxMix/max(beamData.meanSqrtDiodeMaxMix),'LineWidth',2);
plot(beamData.inputPowerMinMixdBm, beamData.meanSqrtDiodeMinMix/max(beamData.meanSqrtDiodeMinMix),'LineWidth',2);
legend('SigGenMean','BeamMaxMix','BeamMinMix');
title('sqrt(Diode) vs. Input Power')
xlabel('Input Power [dBm]');
ylabel('Output [normalised to max]');
grid on;

figure;
semilogy(sigGenData.powerLevels, sigGenData.maxMixers(monIndex,:)/sigGenData.maxMixers(monIndex,sigGenNormIndex),'LineWidth',2);
hold all;
semilogy(sigGenData.powerLevels, sigGenData.minMixers(monIndex,:)/sigGenData.minMixers(monIndex,sigGenNormIndex),'LineWidth',2);

semilogy(beamData.inputPowerMaxMixdBm,beamData.meanMixersMax/beamData.meanMixersMax(beamNormIndex),'LineWidth',2);
semilogy(beamData.inputPowerMinMixdBm, beamData.meanMixersMin/beamData.meanMixersMin(beamNormIndex),'LineWidth',2);
legend('SigGenMaxMix','SigGenMinMix','BeamMaxMix','BeamMinMix');
title('Mixer vs. Input Power (log axis)')
xlabel('Input Power [dBm]');
ylabel('Output [normalised to max]');
grid on;


figure;
plot(sqrt(beamData.inputPowerMaxMixWatts),sqrt(beamData.meanDiodesMaxMix)/max(sqrt(beamData.meanDiodesMaxMix)),'-o')
hold all;
plot(sqrt(beamData.inputPowerMaxMixWatts),beamData.meanDiodesMaxMix/max(beamData.meanDiodesMaxMix),'-o')

% mixer vs. diode
figure;
plot(beamData.meanDiodesMaxMix,beamData.meanMixersMax,'-o');
consts = polyfit(beamData.meanDiodesMaxMix(end-5:end),beamData.meanMixersMax(end-5:end),1);
hold all;
plot(beamData.meanDiodesMaxMix,consts(1).*beamData.meanDiodesMaxMix+consts(2),'-o');
