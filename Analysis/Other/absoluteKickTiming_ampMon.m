dataSets = {...
    '20151124_2017_BeamPickUp_DAC_200_K1Delay_0_K2Delay_1',...
    '20151124_2019_BeamPickUp_DAC_200_K1Delay_3_K2Delay_4',...
    '20151124_2021_BeamPickUp_DAC_200_K1Delay_6_K2Delay_7',...    
    '20151124_2033_BeamPickUp_DAC_200_K1Delay_7_K2Delay_8',...
    '20151124_2034_BeamPickUp_DAC_200_K1Delay_8_K2Delay_9',...
    '20151124_2012_BeamPickUp_DAC_200_K1Delay_9_K2Delay_10',...
    '20151124_2035_BeamPickUp_DAC_200_K1Delay_10_K2Delay_11',...
    '20151124_2036_BeamPickUp_DAC_200_K1Delay_11_K2Delay_12',...
    '20151124_2022_BeamPickUp_DAC_200_K1Delay_12_K2Delay_13',...
    '20151124_2024_BeamPickUp_DAC_200_K1Delay_15_K2Delay_16',...
    '20151124_2025_BeamPickUp_DAC_200_K1Delay_18_K2Delay_19',...
    '20151124_2026_BeamPickUp_DAC_200_K1Delay_21_K2Delay_22',...
    '20151124_2029_BeamPickUp_DAC_200_K1Delay_24_K2Delay_25',...
    '20151124_2030_BeamPickUp_DAC_200_K1Delay_27_K2Delay_28',...
    '20151124_2031_BeamPickUp_DAC_200_K1Delay_30_K2Delay_31',...
};

saveDir = '/home/jack/PhaseFeedforward/Analysis/201511/Plots/20151124_2017_KickDelay';

delays = [0 3 6 7 8 9 10 11 12 15 18 21 24 27 30];

ampMonToUse = 4;

%%
nDataSets = length(dataSets);

for i=1:nDataSets
    CTFData = loadMergedData(dataSets{i});
    if i==1
        nPulses = length(CTFData);
        nSamples = length(CTFData(1).CT_SCOPE02_CH08.Acquisition.value.value);
        meanAmp = NaN(nDataSets,nSamples);
        meanAmp_err = NaN(nDataSets,nSamples);
        sampInt = CTFData(1).CT_SCOPE02_CH08.Acquisition.sampleInterval.value;
    end
    
    ampMon = extractAmpMon(CTFData);
    [meanAmp(i,:),~,meanAmp_err(i,:)] = nanMeanStdErr(squeeze(ampMon(ampMonToUse,:,:)),1);
end

%%
cols = varycolor(nDataSets);

figure;
for i=1:nDataSets
    plot(meanAmp(i,:),'Color',cols(i,:));
    hold all;
end
xlabel('Sample No.')
ylabel('Amplifier Output [V]')
title('Amplifier Output Delay Scan: Whole Pulse')
set(gcf, 'Colormap', cols);
figColBar = colorbar;
allTicks = linspace(0,1,nDataSets+1);
set(figColBar,'YTick',allTicks(1:nDataSets));
set(figColBar,'YTickLabel',delays);
ylabel(figColBar,'Delay [Clock Cycles]');
xlim([220 550])
ylim([-40 40])
% xlim([220 550])
% ylim([-70 30])
format_plots;
% savePlot(saveDir,'allTraces_fullPulse');

figure;
for i=1:nDataSets
    plot(meanAmp(i,:),'Color',cols(i,:),'LineWidth',1.5);
    hold all;
end
xlabel('Sample No.')
ylabel('Amplifier Output [V]')
title('Amplifier Output Delay Scan: End of Pulse')
set(gcf, 'Colormap', cols);
figColBar = colorbar;
allTicks = linspace(0,1,nDataSets+1);
set(figColBar,'YTick',allTicks(1:nDataSets));
set(figColBar,'YTickLabel',delays);
ylabel(figColBar,'Delay [Clock Cycles]');
xlim([485 510])
ylim([-40 40])
% xlim([470 495])
% ylim([-70 30])
format_plots;
% savePlot(saveDir,'allTraces_endPulse');


for datSetToPlot=1:nDataSets
    figure;
    plot(meanAmp(datSetToPlot,:),'Color',cols(datSetToPlot,:),'LineWidth',2);
    xlabel('Sample No.')
    ylabel('Amplifier Output [V]')
    title(sprintf('Output Delay: %d clock cycles',delays(datSetToPlot)))
    xlim([485 510])
    ylim([-40 40])
%     xlim([470 495])
%     ylim([-70 30])
    format_plots;
    % savePlot(saveDir,sprintf('endPulse_delay_%d',delays(datSetToPlot)));

    figure;
    plot(meanAmp(datSetToPlot,:),'Color',cols(datSetToPlot,:),'LineWidth',2);
    xlabel('Sample No.')
    ylabel('Amplifier Output [V]')
    title(sprintf('Output Delay: %d clock cycles',delays(datSetToPlot)))
    xlim([220 550])
    ylim([-40 40])
%     xlim([220 550])
%     ylim([-70 30])
    format_plots;
%     savePlot(saveDir,sprintf('fullPulse_delay_%d',delays(datSetToPlot)));
    pause(3)
end

% for i=1:nDataSets
%     figure;
%     plot(meanAmp(i,:))
%     title(delays(i))
%     xlim([485 510])
% %     input('next...')
% end


%%
% noKickDatSet = '20151124_1959_BeamPickUp_TrigOutDisabled';
% CTFData = loadMergedData(noKickDatSet);
% ampMon = extractAmpMon(CTFData);
% [meanNoKick,~,meanNoKick_err] = nanMeanStdErr(squeeze(ampMon(4,:,:)),1);
% 
% figure;
% plot(8*((1:500)-243),meanNoKick,'LineWidth',2);
% xlabel('Time [ns]')
% ylabel('Amplifier Output [V]')
% title('Beam Pickup On Amplifier Monitoring Signals')
% xlim([-200 1400])
% ylim([-30 30])
% set(gca,'XTick',-200:200:1328)
% format_plots;
% % savePlot(saveDir,'beamPickup_noKick');
