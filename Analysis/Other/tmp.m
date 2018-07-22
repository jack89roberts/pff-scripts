% load('/home/jack/PhaseFeedforward/CTFData/201511/20151124_1911_Nominal_R56_0.08_Gain710_Interleaved/Merged/20151124_1911_Nominal_R56_0.08_Gain710_Interleaved.mat');
% % load('/home/jack/PhaseFeedforward/CTFData/201511/20151124_1658_WiggleAll_Gauss_0.75_Gain700_Interleaved/Merged/20151124_1658_WiggleAll_Gauss_0.75_Gain700_Interleaved.mat');
% % load('/home/jack/PhaseFeedforward/CTFData/201511/20151124_1156_Gain700_NaturalFrascatiChicane_Interleaved/Merged/20151124_1156_Gain700_NaturalFrascatiChicane_Interleaved.mat');
% ampSigs = extractAmpMon(CTFData);
% %%
% LA_odd = squeeze(ampSigs(1,1:2:end,:));
% LB_odd = squeeze(ampSigs(2,1:2:end,:));
% RA_odd = squeeze(ampSigs(3,1:2:end,:));
% RB_odd = squeeze(ampSigs(4,1:2:end,:));
% 
% LA_even = squeeze(ampSigs(1,2:2:end,:));
% LB_even = squeeze(ampSigs(2,2:2:end,:));
% RA_even = squeeze(ampSigs(3,2:2:end,:));
% RB_even = squeeze(ampSigs(4,2:2:end,:));
% 
% [mixers,diodes] = extractMixerDiode(CTFData);
% mix_odd = squeeze(mixers(2,1:2:end,:));
% mix_even = squeeze(mixers(2,2:2:end,:));
% 
% [N,~] = size(mix_odd);
% 
% % figure;
% % for i=1:N
% %     subplot(1,2,1)
% %     plot(LA_odd(i,:))
% %     title('LA')
% %     subplot(1,2,2)
% %     plot(mix_odd(i,:))
% %     title('MIX')
% %     input(sprintf('i = %d...',i))
% % end
% 
% 
% figure;
% plot(8*(1:500),LB_odd(32,:))
% hold all;
% plot(5.2083*(1:768),delaySignal(mix_odd(32,:),93)*13+0.6)
% legend('LB','MIX')
% hold off
% 
% figure;
% plot(8*(1:500),LB_even(1,:)./(max(abs(LB_even(1,:)))))
% hold all;
% plot(5.2083*(1:768),delaySignal(mix_even(1,:),95)./(max(abs(mix_even(1,:)))))
% legend('LB','MIX')
% 
% 
% figure;
% plot(LB_odd(32,:)./(max(abs(LB_odd(32,:)))))
% hold all;
% plot(delaySignal(mix_odd(32,:),95)./(max(abs(mix_odd(32,:)))))
% legend('LB','MIX')

%%
load('/home/jack/PhaseFeedforward/CTFData/201512/20151202_1829_FF_Gain_800_Interleaved/Merged/20151202_1829_FF_Gain_800_Interleaved.mat');
load('/home/jack/PhaseFeedforward/CTFData/201512/20151202_1947_FFScan_Gain_-1200_Interleaved/Merged/20151202_1947_FFScan_Gain_-1200_Interleaved.mat');
calUpA = 1.083143;
calUpO = 0.049037;
calDownA = 0.993613;
calDownO = 0.023793;

ampMon = double(extractCTFSignalFromMergedData('CT_SCOPE02_CH07.Acquisition.value.value',CTFData));
[mix,dio] = extractMixerDiode(CTFData);

upPhase = getPhaseMixerDiode(squeeze(mix(2,:,:)),dio,calUpA,calUpO,0);
downPhase = getPhaseMixerDiode(squeeze(mix(3,:,:)),dio,calDownA,calDownO,0);

[up,~,up_err] = nanMeanStdErr(upPhase,1);
[evenDown,~,evenDown_err] = nanMeanStdErr(downPhase(2:2:end,:),1);
[oddDown,~,oddDown_err] = nanMeanStdErr(downPhase(1:2:end,:),1);
[oddAmp,~,oddAmp_err] = nanMeanStdErr(ampMon(1:2:end,:),1);
[evenAmp,~,evenAmp_err] = nanMeanStdErr(ampMon(2:2:end,:),1);

up = delaySignal(up,124);
up_err = delaySignal(up_err,124);

oddAmp = delaySignal(oddAmp,17);
oddAmp_err = delaySignal(oddAmp_err,17);


diffDown = evenDown-oddDown;
diffDown_err = sqrt(evenDown_err.^2 + oddDown_err.^2);

figure;
shadedErrorBar(1:768,up,up_err,'b');
hold all;
shadedErrorBar(1:768,diffDown-2,diffDown_err,'r');
% shadedErrorBar(1:768,oddAmp/3000-2,oddAmp_err/3000,'g')
shadedErrorBar(1:768,evenAmp/3000-2,evenAmp_err/3000,'g')
title('blue=up, red=diff down, green=amp')
ylabel('phase');
xlabel('sample');
format_plots;
% savePlot('/home/jack/Documents/pff-thesis/Figures/commissioning','pffKickShape')

figure;
plot(up);
hold all;
plot(diffDown-2);