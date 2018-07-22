addpath('/home/jack/Documents/MATLAB/PhaseFeedforward/ctfmod/PhaseFeedforward');

load('/home/jack/PhaseFeedforward/CTFData/201511/20151124_2017_BeamPickUp_DAC_200_K1Delay_0_K2Delay_1/Merged/20151124_2017_BeamPickUp_DAC_200_K1Delay_0_K2Delay_1.mat')
ampSigs = extractAmpMon(CTFData);
ampSigs = squeeze(nanmean(ampSigs,2));

figure;
plot(ampSigs(1,:),'b','LineWidth',2)
hold all;

load('/home/jack/PhaseFeedforward/CTFData/201511/20151124_2012_BeamPickUp_DAC_200_K1Delay_9_K2Delay_10/Merged/20151124_2012_BeamPickUp_DAC_200_K1Delay_9_K2Delay_10.mat')
ampSigs = extractAmpMon(CTFData);
ampSigs = squeeze(nanmean(ampSigs,2));
plot(ampSigs(1,:),'r','LineWidth',2)

load('/home/jack/PhaseFeedforward/CTFData/201511/20151124_2029_BeamPickUp_DAC_200_K1Delay_24_K2Delay_25/Merged/20151124_2029_BeamPickUp_DAC_200_K1Delay_24_K2Delay_25.mat')
ampSigs = extractAmpMon(CTFData);
ampSigs = squeeze(nanmean(ampSigs,2));
plot(ampSigs(1,:),'Color',[0 0.8 0],'LineWidth',2)

% load('/home/jack/PhaseFeedforward/CTFData/201511/20151124_2021_BeamPickUp_DAC_200_K1Delay_6_K2Delay_7/Merged/20151124_2021_BeamPickUp_DAC_200_K1Delay_6_K2Delay_7.mat');
% ampSigs = extractAmpMon(CTFData);
% ampSigs = squeeze(nanmean(ampSigs,2));
% plot(ampSigs(1,:))
 
% load('/home/jack/PhaseFeedforward/CTFData/201511/20151124_2025_BeamPickUp_DAC_200_K1Delay_18_K2Delay_19/Merged/20151124_2025_BeamPickUp_DAC_200_K1Delay_18_K2Delay_19.mat');
% ampSigs = extractAmpMon(CTFData);
% ampSigs = squeeze(nanmean(ampSigs,2));
% plot(ampSigs(1,:))

%legend('Delay 0ns','Delay 17ns','Delay 50ns');
legend('Delay 0ns','Delay 25ns','Delay 65ns');
xlabel('Sample No.')
ylabel('Output [V]')
format_plots;
plot([479 479],get(gca,'YLim'),'--k','LineWidth',1.5);
grid off;

load('/home/jack/PhaseFeedforward/CTFData/201511/20151124_1959_BeamPickUp_TrigOutDisabled/Merged/20151124_1959_BeamPickUp_TrigOutDisabled.mat')
ampSigs = extractAmpMon(CTFData);
ampSigs = squeeze(nanmean(ampSigs,2));
figure;
plot(ampSigs(1,:),'LineWidth',2)
format_plots;
xlabel('Sample No.')
ylabel('Output [V]')
title('Beam Pickup from Kicker')
