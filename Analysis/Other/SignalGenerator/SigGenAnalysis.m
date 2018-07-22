%%
close all;clearvars;
%% input parameters
mix1RefOffsetFile = '20150731_1345_SigGen_0dbM_RfOff';
mix1DataSetNames = {
    '20150731_1347_SigGen_0dbM_RfOn',...
    '20150731_1349_SigGen_3dbM_RfOn',...
    '20150731_1352_SigGen_6dbM_RfOn',...
    '20150731_1354_SigGen_9dbM_RfOn',...
    '20150731_1357_SigGen_12dbM_RfOn',...
    '20150731_1358_SigGen_15dbM_RfOn',...
    '20150731_1404_SigGen_18dbM_RfOn',...
    '20150731_1406_SigGen_21dbM_RfOn',...
    '20150731_1407_SigGen_24dbM_RfOn',...
    '20150731_1409_SigGen_27dbM_RfOn',...
    '20150731_1411_SigGen_30dbM_RfOn',...
    '20150731_1412_SigGen_33dbM_RfOn',...
};

mix2RefOffsetFile = '20150805_1727_SigGen_Mix2_0dbM_RfOff';
mix2DataSetNames = {
    '20150805_1730_SigGen_Mix2_0dbM_RfOn',...
    '20150805_1733_SigGen_Mix2_3dbM_RfOn',...
    '20150805_1735_SigGen_Mix2_6dbM_RfOn',...
    '20150805_1736_SigGen_Mix2_9dbM_RfOn',...
    '20150805_1737_SigGen_Mix2_12dbM_RfOn',...
    '20150805_1738_SigGen_Mix2_15dbM_RfOn',...
    '20150805_1740_SigGen_Mix2_18dbM_RfOn',...
    '20150805_1740_SigGen_Mix2_21dbM_RfOn',...
    '20150805_1742_SigGen_Mix2_24dbM_RfOn',...
    '20150805_1743_SigGen_Mix2_27dbM_RfOn',...
    '20150805_1745_SigGen_Mix2_30dbM_RfOn',...
    '20150805_1746_SigGen_Mix2_33dbM_RfOn',...
};

mix3RefOffsetFile = '20150805_1759_SigGen_Mix3_0dbM_RfOff';
mix3DataSetNames = {
    '20150805_1800_SigGen_Mix3_0dbM_RfOn',...
    '20150805_1801_SigGen_Mix3_3dbM_RfOn',...
    '20150805_1802_SigGen_Mix3_6dbM_RfOn',...
    '20150805_1804_SigGen_Mix3_9dbM_RfOn',...
    '20150805_1805_SigGen_Mix3_12dbM_RfOn',...
    '20150805_1806_SigGen_Mix3_15dbM_RfOn',...
    '20150805_1807_SigGen_Mix3_18dbM_RfOn',...
    '20150805_1808_SigGen_Mix3_21dbM_RfOn',...
    '20150805_1809_SigGen_Mix3_24dbM_RfOn',...
    '20150805_1810_SigGen_Mix3_27dbM_RfOn',...
    '20150805_1811_SigGen_Mix3_30dbM_RfOn',...
    '20150805_1812_SigGen_Mix3_33dbM_RfOn',...
};

nMixers = 3;
nPulses = 30;
nSamples = 384;
powerLevels = 0:3:33; % dBm

mixerPowerLevelsToFit = 3:7; % indices
diodePowerLevelsToFit = 3:6;

mixerGoodPlotRange = [1 9];  %indices
diodeGoodPlotRange = [1 6];

savePlots = 0;
savePlotDir = '/home/jack/PhaseFeedforward/Analysis/201508/SignalGenerator';

%% some path/arrays definitions
addpath('../../../'); % PhaseFeedforward directory

nPowerLevels = length(mix1DataSetNames);
dataSetNames = cell(nMixers,length(mix1DataSetNames));
dataSetNames(1,:) = mix1DataSetNames;
dataSetNames(2,:) = mix2DataSetNames;
dataSetNames(3,:) = mix3DataSetNames;

refOffsetFileNames = {mix1RefOffsetFile,mix2RefOffsetFile,mix3RefOffsetFile};

%% extract mixers and didoes from data files

powerInWatts = dBmToWatts(powerLevels);
powerInVolts = dBmToVolts(powerLevels);

mixers = NaN(nMixers,nPowerLevels,nPulses,nSamples);
diodes = NaN(nMixers,nPowerLevels,nPulses,nSamples);

for m=1:nMixers
    fprintf('------------------- MIXER %d -------------------\n',m);
    
    % load data
    [refData,~,~] = loadMergedData(refOffsetFileNames{m});
    data = cell(1,nPowerLevels);
    for d=1:nPowerLevels
        fprintf('Loading %s...\n',dataSetNames{m,d});
        [data{d},~,~] = loadMergedData(dataSetNames{m,d});
    end

    % extract mixers and diodes
    fprintf('Extracting mixers and diodes...\n');

    [mixerBaseline,diodeBaseline] = extractMixerDiode(refData,1);
    mixerBaseline = squeeze(mixerBaseline(m,:,:));
    diodeBaseline = squeeze(diodeBaseline(m,:,:));
    mixerBaseline = nanmean(mixerBaseline(:));
    diodeBaseline = nanmean(diodeBaseline(:));

    for d=1:nPowerLevels
        [tmpMixers,tmpDiodes] = extractMixerDiode(data{d},1);
        tmpMixers = squeeze(tmpMixers(m,:,:)) - mixerBaseline;
        tmpDiodes = squeeze(tmpDiodes(m,:,:)) - diodeBaseline;
        mixers(m,d,:,:) = tmpMixers(1:nPulses,:);
        diodes(m,d,:,:) = tmpDiodes(1:nPulses,:);
    end
end

%% align mixers/diodes and take mean

fprintf('Aligning signals...\n');

% align signals within each power level
for m=1:nMixers
    for d=1:nPowerLevels
        [mixers(m,d,:,:),diodes(m,d,:,:)] = getAlignedCont(squeeze(mixers(m,d,:,:)),squeeze(diodes(m,d,:,:)));
    end
end

mixers = squeeze(nanmean(mixers,3));
diodes = squeeze(nanmean(diodes,3));

% align mean signals between power levels
for m=1:nMixers
    [mixers(m,:,:),diodes(m,:,:)] = getAlignedCont(squeeze(mixers(m,:,:)),squeeze(diodes(m,:,:)));
end

% % align mean signals between mixers
% for d=1:nPowerLevels
%     [mixers(:,d,:),diodes(:,d,:)] = getAlignedCont(squeeze(mixers(:,d,:)),squeeze(diodes(:,d,:)));
% end

% want diode to be positive
mixers = -mixers;
diodes = -diodes;

meanDiodes = squeeze(nanmean(diodes,3));
sqrtMeanDiodes = sqrt(meanDiodes);
maxDiodes = max(diodes,[],3);
minDiodes = min(diodes,[],3);
diffMaxMinDiodes = (maxDiodes-minDiodes)/2;

maxMixers = max(mixers,[],3);
minMixers = min(mixers,[],3);
meanAmpMixers = (maxMixers-minMixers)/2;

%% Initial plots
fprintf('Making plots...\n');

mixerColours = [0, 0, 1; 1, 0, 0; 0, 0.5, 0];
powerColours = varycolor(nPowerLevels);

legLabels = cell(1,nPowerLevels);
for m=1:nMixers
    figure;
    for d=1:nPowerLevels
        plot(squeeze(diodes(m,d,:)),'Color',powerColours(d,:));
        hold all;
        legLabels{d} = [num2str(powerLevels(d)) ' dBm'];
    end
    xlim([round(nSamples/10) 9*round(nSamples/10)]);
    ylim([0 2]);
    title(sprintf('Diode %d',m));
    xlabel('Sample No.');
    ylabel('Output [V]');
    grid on;
    
    set(gcf, 'Colormap', powerColours);
    figColBar = colorbar;
    set(figColBar,'YTick',(1:nPowerLevels)/nPowerLevels);
    set(figColBar,'YTickLabel',powerLevels);   
    ylabel(figColBar,'Input Power [dBm]');
    
    if savePlots; savePlot(savePlotDir,sprintf('Diode%d_AllPowerLevels',m)); end;
end

legLabels = cell(1,nPowerLevels);
for m=1:nMixers
    figure;
    for d=1:nPowerLevels
        plot(squeeze(mixers(m,d,:)),'Color',powerColours(d,:));
        hold all;
        legLabels{d} = [num2str(powerLevels(d)) ' dBm'];
    end
    xlim([round(nSamples/10) 9*round(nSamples/10)]);
    ylim([-0.6 0.6]);
    title(sprintf('Mixer %d',m));
    xlabel('Sample No.');
    ylabel('Output [V]');
    grid on;
    
    set(gcf, 'Colormap', powerColours);
    figColBar = colorbar;
    set(figColBar,'YTick',(1:nPowerLevels)/nPowerLevels);
    set(figColBar,'YTickLabel',powerLevels);   
    ylabel(figColBar,'Input Power [dBm]');
    
    if savePlots; savePlot(savePlotDir,sprintf('Mixer%d_AllPowerLevels',m)); end;
end

figure; % diode vs. input volts
for m=1:nMixers
    plot(powerInVolts,meanDiodes(m,:),'-o','LineWidth',2,'Color',mixerColours(m,:));
    hold all;
end
legend('Diode 1','Diode 2','Diode 3');
title('Diode Output vs. Input Voltage')
xlabel('Input [V]');
ylabel('Output [V]');
grid on;
if savePlots; savePlot(savePlotDir,'MeanDiodeVsVolts'); end;

figure; % sqrt(diode) vs. input volts
for m=1:nMixers
    plot(powerInVolts,sqrtMeanDiodes(m,:),'-o','LineWidth',2,'Color',mixerColours(m,:));
    hold all;
end
legend('Diode 1','Diode 2','Diode 3');
title('sqrt(Diode) vs. Input Voltage')
xlabel('Input [V]');
ylabel('Output [sqrt(V)]');
grid on;
if savePlots; savePlot(savePlotDir,'SqrtDiodeVsVolts'); end;

figure; % mixer vs. volts
for m=1:nMixers
    plot(powerInVolts,meanAmpMixers(m,:),'-o','LineWidth',2,'Color',mixerColours(m,:));
    hold all;
    plot(powerInVolts,maxMixers(m,:),'--^','Color',mixerColours(m,:));
    plot(powerInVolts,-minMixers(m,:),'--v','Color',mixerColours(m,:));

end
legend( 'Mixer 1 (mean amplitude)',...
        'Mixer 1 (amplitude at max)',...
        'Mixer 1 (amplitude at min)',...
        'Mixer 2 (mean amplitude)',...
        'Mixer 2 (amplitude at max)',...
        'Mixer 2 (amplitude at min)',...
        'Mixer 3 (mean amplitude)',...
        'Mixer 3 (amplitude at max)',...
        'Mixer 3 (amplitude at min)',...
        'Location','best');
title('Mixer Output vs. Input Voltage')
xlabel('Input [V]');
ylabel('Output [V]');
grid on;
if savePlots; savePlot(savePlotDir,'MixerVsVolts'); end;

figure; % mixer vs. volts
for m=1:nMixers
    plot(powerInVolts,(abs(maxMixers(m,:))-abs(minMixers(m,:))),'o-','Color',mixerColours(m,:));
    hold all;
end
legend( 'Mixer 1',...
        'Mixer 2',...
        'Mixer 3',...
        'Location','best');
title('Mixer Max-Min vs. Input Voltage')
xlabel('Input [V]');
ylabel('Output [V]');
grid on;
if savePlots; savePlot(savePlotDir,'DiffMaxMinMixerVsVolts'); end;

figure; % mixer vs. sqrt(diode)
for m=1:nMixers
    plot(sqrtMeanDiodes(m,:),meanAmpMixers(m,:),'-o','LineWidth',2,'Color',mixerColours(m,:));
    hold all;
end
legend( 'Mixer 1',...
        'Mixer 2',...
        'Mixer 3',...
        'Location','best');
title('Mixer vs. sqrt(Diode)')
ylabel('Mixer [V]');
xlabel('sqrt(Diode) [sqrt(V)]');
grid on;
if savePlots; savePlot(savePlotDir,'MixerVsSqrtDopde'); end;

figure; % mixer vs. diode
for m=1:nMixers
    plot(meanDiodes(m,:),meanAmpMixers(m,:),'-o','LineWidth',2,'Color',mixerColours(m,:));
    hold all;
end
legend( 'Mixer 1',...
        'Mixer 2',...
        'Mixer 3',...
        'Location','best');
title('Mixer vs. Diode')
ylabel('Mixer [V]');
xlabel('Diode [V]');
grid on;
if savePlots; savePlot(savePlotDir,'MixerVsDopde'); end;

% figure; % dBm vs log(output)
% for m=1:nMixers
%     semilogy(powerLevels,meanDiodes(m,:),'LineWidth',2,'Color',mixerColours(m,:));
%     hold all;
% end
% legend('Diode 1','Diode 2','Diode 3');
% title('Diode Output vs. Input Power (Log Axis)')
% xlabel('Input Power [dBm]');
% ylabel('Output [V]');
% grid on;
% if savePlots; savePlot(savePlotDir,'MeanDiodeVsPower_Log'); end;

% figure; % dBm vs. output
% for m=1:nMixers
%     plot(powerLevels,meanDiodes(m,:),'LineWidth',2,'Color',mixerColours(m,:));
%     hold all;
% end
% legend('Diode 1','Diode 2','Diode 3');
% title('Diode Output vs. Input Power')
% xlabel('Input Power [dBm]');
% ylabel('Output [V]');
% grid on;
% if savePlots; savePlot(savePlotDir,'MeanDiodeVsPower'); end;

% figure; % dBm vs log(output)
% for m=1:nMixers
%     semilogy(powerLevels,meanDiodes(m,:),'LineWidth',2,'Color',mixerColours(m,:));
%     hold all;
% end
% legend('Diode 1','Diode 2','Diode 3');
% title('Diode Output vs. Input Power (Log Axis)')
% xlabel('Input Power [dBm]');
% ylabel('Output [V]');
% grid on;
% if savePlots; savePlot(savePlotDir,'MeanDiodeVsPower_Log'); end;

% figure; % W vs. output
% for m=1:nMixers
%     plot(powerInWatts,meanDiodes(m,:),'LineWidth',2,'Color',mixerColours(m,:));
%     hold all;
% end
% legend('Diode 1','Diode 2','Diode 3');
% title('Diode Output vs. Input Power')
% xlabel('Input Power [W]');
% ylabel('Output [V]');
% grid on;
% if savePlots; savePlot(savePlotDir,'MeanDiodeVsPowerWatts'); end;

% figure;
% for m=1:nMixers
%     plot(powerLevels,sqrtMeanDiodes(m,:),'LineWidth',2,'Color',mixerColours(m,:));
%     hold all;
% end
% legend('Diode 1','Diode 2','Diode 3');
% title('sqrt(Diode) Output vs. Input Power')
% xlabel('Input Power [dBm]');
% ylabel('Output [sqrt(V)]');
% grid on;
% if savePlots; savePlot(savePlotDir,'sqrtMeanDiodeVsPower'); end;

% figure;
% for m=1:nMixers
%     plot(powerInWatts,sqrtMeanDiodes(m,:),'LineWidth',2,'Color',mixerColours(m,:));
%     hold all;
% end
% legend('Diode 1','Diode 2','Diode 3');
% title('sqrt(Diode) Output vs. Input Power')
% xlabel('Input Power [W]');
% ylabel('Output [sqrt(V)]');
% grid on;
% if savePlots; savePlot(savePlotDir,'sqrtMeanDiodeVsPowerWatts'); end;

% figure;
% for m=1:nMixers
%     semilogy(powerLevels,sqrtMeanDiodes(m,:),'LineWidth',2,'Color',mixerColours(m,:));
%     hold all;
% end
% legend('Diode 1','Diode 2','Diode 3');
% title('sqrt(Diode) Output vs. Input Power (Log Axis)')
% xlabel('Input Power [dBm]');
% ylabel('Output [sqrt(V)]');
% grid on;
% if savePlots; savePlot(savePlotDir,'sqrtMeanDiodeVsPower_Log'); end;

% figure;
% for m=1:nMixers
%     plot(powerLevels,meanAmpMixers(m,:),'-o','LineWidth',2,'Color',mixerColours(m,:));
%     hold all;
%     plot(powerLevels,maxMixers(m,:),'--^','Color',mixerColours(m,:));
%     plot(powerLevels,-minMixers(m,:),'--v','Color',mixerColours(m,:));
% end
% legend( 'Mixer 1 (mean amplitude)',...
%         'Mixer 1 (amplitude at max)',...
%         'Mixer 1 (amplitude at min)',...
%         'Mixer 2 (mean amplitude)',...
%         'Mixer 2 (amplitude at max)',...
%         'Mixer 2 (amplitude at min)',...
%         'Mixer 3 (mean amplitude)',...
%         'Mixer 3 (amplitude at max)',...
%         'Mixer 3 (amplitude at min)',...
%         'Location','best');
% title('Mixer Output vs. Input Power')
% xlabel('Input Power [dBm]');
% ylabel('Output [V]');
% grid on;
% if savePlots; savePlot(savePlotDir,'MixerVsPower'); end;

% figure;
% for m=1:nMixers
%     plot(powerInWatts,meanAmpMixers(m,:),'-o','LineWidth',2,'Color',mixerColours(m,:));
%     hold all;
%     plot(powerInWatts,maxMixers(m,:),'--^','Color',mixerColours(m,:));
%     plot(powerInWatts,-minMixers(m,:),'--v','Color',mixerColours(m,:));
% end
% legend( 'Mixer 1 (mean amplitude)',...
%         'Mixer 1 (amplitude at max)',...
%         'Mixer 1 (amplitude at min)',...
%         'Mixer 2 (mean amplitude)',...
%         'Mixer 2 (amplitude at max)',...
%         'Mixer 2 (amplitude at min)',...
%         'Mixer 3 (mean amplitude)',...
%         'Mixer 3 (amplitude at max)',...
%         'Mixer 3 (amplitude at min)',...
%         'Location','best');
% title('Mixer Output vs. Input Power')
% xlabel('Input Power [W]');
% ylabel('Output [V]');
% grid on;
% if savePlots; savePlot(savePlotDir,'MixerVsPowerWatts'); end;

% figure;
% for m=1:nMixers
%     semilogy(powerLevels,meanAmpMixers(m,:),'-o','LineWidth',2,'Color',mixerColours(m,:));
%     hold all;
%     semilogy(powerLevels,maxMixers(m,:),'--^','Color',mixerColours(m,:));
%     semilogy(powerLevels,-minMixers(m,:),'--v','Color',mixerColours(m,:));
% end
% legend( 'Mixer 1 (mean amplitude)',...
%         'Mixer 1 (amplitude at max)',...
%         'Mixer 1 (amplitude at min)',...
%         'Mixer 2 (mean amplitude)',...
%         'Mixer 2 (amplitude at max)',...
%         'Mixer 2 (amplitude at min)',...
%         'Mixer 3 (mean amplitude)',...
%         'Mixer 3 (amplitude at max)',...
%         'Mixer 3 (amplitude at min)',...
%         'Location','best');
% title('Mixer Output vs. Input Power (Log Axis)')
% xlabel('Input Power [dBm]');
% ylabel('Output [V]');
% grid on;
% if savePlots; savePlot(savePlotDir,'MixerVsPower_Log'); end;

%% Calibrations - mixer
fprintf('Calculating Mixer calibrations...\n');

% estimate frequency (number of samples per 360 degrees)
estFreq = NaN(nMixers,nPowerLevels);
for m=1:nMixers
    for d=1:nPowerLevels
        [pks,locs] = findpeaks(squeeze(mixers(m,d,:)),'MinPeakHeight',0,'NPeaks',7,'SortStr','descend'); % returns 7 largest peaks sorted from largest to smallest. 7 picked to match no. of peaks in data.
        estFreq(m,d) = nanmean(diff(sort(locs))); % average of the difference in peak locations. Sort needed as peaks are returned in order from largest to smallest but we want difference between neighbouring peaks.
    end
end
estFreq = nanmean(removeBadPulses(estFreq(:))); % overall best estimate, mean removing any outliers.

mixerCalConsts = NaN(nMixers,nPowerLevels,4);
mixerCalRSquare = NaN(nMixers,nPowerLevels);

for m=1:nMixers
    fprintf('Mixer %d...\n',m);

    for d=1:nPowerLevels        
        tmpMixer = squeeze(mixers(m,d,:));
        isGood = ~isnan(tmpMixer);
        tmpMixer = tmpMixer(isGood);
        tmpSamples = 1:nSamples;
        tmpSamples = tmpSamples(isGood)';
        
        initA = (maxMixers(m,d) - minMixers(m,d))./2; % est amplitude
        initB = 2*pi/estFreq; % est frequency, s = no. samples between neighbouring max and min (180 degrees)
        initD = (maxMixers(m,d) + minMixers(m,d))./2; % est vertical offset
        
        sinArg = (tmpMixer-initD)./initA;
        sinArg(sinArg>1)= NaN; % avoid imaginary component if data point above amplitude guess
        sinArg(sinArg<-1)= NaN;
        sinArg = asin(sinArg);
        sinArg(diff(tmpMixer)<0) = NaN;
        initC = sinArg-(initB.*tmpSamples);
        initC = mod(initC,2*pi); % make everything between 0 and 2pi
        initC(initC>pi) = initC(initC>pi)-2*pi;
        initC(initC<-pi) = initC(initC<-pi)+2*pi;
        initC = nanmean(initC);

        initParams = [initA, initB, initC, initD];


            % need to select 1st quadrant where asin valid - positive
            % gradient in difference
            
            initC = mod(initC,2*pi); % make everything between 0 and 2pi
            initC(initC>pi) = initC(initC>pi)-2*pi;
            initC(initC<-pi) = initC(initC<-pi)+2*pi;
            initC = nanmean(initC);

        
        
        [calibrationFactors, fitRSquare, fitConfInt] = offsetSinFit(tmpSamples, tmpMixer, initParams);
        
        % correct negative amplitude fit
        mixerCalConsts(m,d,:) = calibrationFactors;           
        mixerCalRSquare(m,d) = fitRSquare;
    end
end

figure;
for m=1:nMixers
    semilogy(powerLevels, squeeze(abs(mixerCalConsts(m,:,1))),'LineWidth',2,'Color',mixerColours(m,:));
    hold all;
end
legend( 'Mixer 1','Mixer 2','Mixer 3','Location','best');
title('Fitted Mixer Output vs. Input Power (Log Axis)')
xlabel('Input Power [dBm]');
ylabel('Output [V]');
grid on;
if savePlots; savePlot(savePlotDir,'FittedMixerVsPower_Log'); end;

figure;
for m=1:nMixers
    plot(powerLevels, mixerCalRSquare(m,:),'LineWidth',2,'Color',mixerColours(m,:));
    hold all;
end
legend( 'Mixer 1','Mixer 2','Mixer 3','Location','best');
title('Fit R square: Mixer Output vs. Input Power')
xlabel('Input Power [dBm]');
ylabel('Rsquare');
grid on;
if savePlots; savePlot(savePlotDir,'FitRSquareMixerVsPower'); end;

mixerStdResidVolt = NaN(nMixers,nPowerLevels);
mixerStdResidPhas = NaN(nMixers,nPowerLevels);

figure;
for m=1:nMixers
    for d=1:nPowerLevels
        fitFunc = @(x) (mixerCalConsts(m,d,1).*sin( (mixerCalConsts(m,d,2).*x) + mixerCalConsts(m,d,3) )) + mixerCalConsts(m,d,4);
        fitX = 1:nSamples;
        fitY = fitFunc(fitX);

        mixerStdResidVolt(m,d) = nanstd(squeeze(mixers(m,d,:))'-fitY);        
        mixerStdResidPhas(m,d) = asin(mixerStdResidVolt(m,d)/mixerCalConsts(m,d,1));

        
        plot(fitX,fitY,'k','LineWidth',1.5);
        hold all;
        plot(squeeze(mixers(m,d,:)),'o','Color',mixerColours(m,:));
        title(sprintf('Mixer %d, Power %d dBm (Amplitude: %.2d V, Rsquare: %.3f)',m,powerLevels(d),abs(mixerCalConsts(m,d,1)),mixerCalRSquare(m,d)));
        xlabel('Sample No.');
        ylabel('Output [V]');
        grid on;
        legend('Fit','Data');
        hold off;
        
        if savePlots
            savePlot([savePlotDir '/Fits'],sprintf('Mixer%d_Power%d',m,powerLevels(d)));
        else
%             input(sprintf('Mixer %d, Power %d dBm',m,powerLevels(d)));
        end

    end
end

figure;
plot(powerLevels,mixerStdResidVolt','LineWidth',2);
xlabel('Input Power [dBm]')
ylabel('Std Residuals [V]')
title('Std Deviation of Residuals to Sin Fit')
legend('Mixer 1','Mixer 2','Mixer 3');
xlim([powerLevels(1) powerLevels(end)])
format_plots;

figure;
plot(powerLevels,mixerStdResidPhas','LineWidth',2);
xlabel('Input Power [dBm]')
ylabel('Std Residuals [degrees]')
title('Std Deviation of Residuals to Sin Fit')
legend('Mixer 1','Mixer 2','Mixer 3');
xlim([powerLevels(1) powerLevels(end)])
format_plots;

%% Calibrations - diode

diodeCalConsts = NaN(nMixers,nPowerLevels,4);
diodeCalRSquare = NaN(nMixers,nPowerLevels);
for m=1:nMixers
    fprintf('Diode %d...\n',m);

    for d=1:nPowerLevels        
        tmpDiode = squeeze(diodes(m,d,:));
        isGood = ~isnan(tmpDiode);
        tmpDiode = tmpDiode(isGood);
        tmpSamples = 1:nSamples;
        tmpSamples = tmpSamples(isGood)';
        
        initA =  (abs(maxDiodes(m,d)) - abs(minDiodes(m,d)))./2; % est amplitude
        initB = 2*pi/estFreq; % est frequency, s = no. samples between neighbouring max and min (180 degrees)
        initC = (tmpDiode(1)./initA)*(pi/2); % est phase shift
        initD = meanDiodes(m,d); % est vertical offset
        initParams = [initA, initB, initC, initD];

        [calibrationFactors, fitRSquare, fitConfInt] = offsetSinFit(tmpSamples, tmpDiode, initParams);
    
        % correct negative amplitude fit
        diodeCalConsts(m,d,:) = calibrationFactors;           
        diodeCalRSquare(m,d) = fitRSquare;
    end
end

relativeDiodeXTalk = abs(diodeCalConsts(:,:,1)./meanDiodes);

figure;
for m=1:nMixers
    plot(powerLevels, squeeze(abs(diodeCalConsts(m,:,1))),'LineWidth',2,'Color',mixerColours(m,:));
    hold all;
end
drawnow();
legend( 'Diode 1','Diode 2','Diode 3','Location','best');
title('Fitted Diode Cross-Talk Amplitude vs. Input Power')
xlabel('Input Power [dBm]');
ylabel('Amplitude [V]');
grid on;
if savePlots; savePlot(savePlotDir,'FittedDiodeXTalkVsPower'); end;

figure;
for m=1:nMixers
    plot(powerLevels, diodeCalRSquare(m,:),'LineWidth',2,'Color',mixerColours(m,:));
    hold all;
end
drawnow();
legend( 'Diode 1','Diode 2','Diode 3','Location','best');
title('Fit R square: Mixer Diode Cross-Talk Amplitude vs. Input Power')
xlabel('Input Power [dBm]');
ylabel('Rsquare');
grid on;
if savePlots; savePlot(savePlotDir,'FitRSquareDiodeXTalkVsPower'); end;

figure;
for m=1:nMixers
    plot(powerLevels, sqrt(squeeze(relativeDiodeXTalk(m,:))),'LineWidth',2,'Color',mixerColours(m,:));
    hold all;
end
drawnow();
legend( 'Diode 1','Diode 2','Diode 3','Location','best');
title('Relative Diode Cross-Talk vs. Input Power')
xlabel('Input Power [dBm]');
ylabel('Relative Cross-Talk');
grid on;
if savePlots; savePlot(savePlotDir,'RelativeDiodeXTalkVsPower'); end;

figure;
for m=1:nMixers
    for d=1:nPowerLevels
        fitFunc = @(x) (diodeCalConsts(m,d,1).*sin( (diodeCalConsts(m,d,2).*x) + diodeCalConsts(m,d,3) )) + diodeCalConsts(m,d,4);
        fitX = 1:nSamples;
        fitY = fitFunc(fitX);
        
        plot(fitX,fitY,'k','LineWidth',1.5);
        hold all;
        plot(squeeze(diodes(m,d,:)),'o','Color',mixerColours(m,:));
        title(sprintf('Diode %d, Power %d dBm (Amplitude: %.2d V, Rsquare: %.3f)',m,powerLevels(d),abs(diodeCalConsts(m,d,1)),diodeCalRSquare(m,d)));
        xlabel('Sample No.');
        ylabel('Output [V]');
        grid on;
        drawnow();
        legend('Fit','Data');
        hold off;
        
        if savePlots
            savePlot([savePlotDir '/Fits'],sprintf('Diode%d_Power%d',m,powerLevels(d)));
        else
%             input(sprintf('Diode %d, Power %d dBm',m,powerLevels(d)));
        end
    end
end

%% Linear fits to mixer and diode

linFitMixer = NaN(nMixers,2);
linFitDiode = NaN(nMixers,2);
quadFitDiode = NaN(nMixers,3);
linFitSqrtDiode = NaN(nMixers,2);
linFitMixerVsDiode = NaN(nMixers,2);

linFitMixerPoints = NaN(nMixers,nPowerLevels);
linFitDiodePoints = NaN(nMixers,nPowerLevels);
linFitSqrtDiodePoints = NaN(nMixers,nPowerLevels);
linFitMixerVsDiodePoints = NaN(nMixers,nPowerLevels);
for m=1:nMixers
    tmpMixer = meanAmpMixers(m,mixerPowerLevelsToFit);%mixerCalConsts(m,mixerPowerLevelsToUse,1);
    tmpDiode = meanDiodes(m,diodePowerLevelsToFit);
    tmpSqrtDiode = sqrtMeanDiodes(m,diodePowerLevelsToFit);
    
    linFitMixer(m,:) = polyfit(powerInVolts(mixerPowerLevelsToFit),tmpMixer,1);
    linFitDiode(m,:) = polyfit(powerInVolts(diodePowerLevelsToFit),tmpDiode,1);
    linFitSqrtDiode(m,:) = polyfit(powerInVolts(diodePowerLevelsToFit),tmpSqrtDiode,1);    
    linFitMixerVsDiode(m,:) = polyfit(tmpDiode,meanAmpMixers(m,diodePowerLevelsToFit),1);
    quadFitDiode(m,:) = polyfit(powerInVolts(diodePowerLevelsToFit),tmpDiode,2);
    
    %linFitMixerAmp(m,:) = 10.^(linFitMixer(m,1).*powerLevels + linFitMixer(m,2));
    linFitMixerPoints(m,:) = linFitMixer(m,1).*powerInVolts + linFitMixer(m,2);
    linFitDiodePoints(m,:) = linFitDiode(m,1).*powerInVolts + linFitDiode(m,2); 
    linFitSqrtDiodePoints(m,:) = linFitSqrtDiode(m,1).*powerInVolts + linFitSqrtDiode(m,2); 
    linFitMixerVsDiodePoints(m,:) = linFitMixerVsDiode(m,1)*meanDiodes(m,:) + linFitMixerVsDiode(m,2);
    
end

figure;
for m=1:nMixers
    plot(powerInVolts,meanDiodes(m,:),'-o','LineWidth',1,'Color',mixerColours(m,:));
    hold all;
    fitPowers = linspace(min(powerInVolts),max(powerInVolts),1000);
    fitValues = polyval(quadFitDiode(m,:),fitPowers);
    plot(fitPowers,fitValues,'--','LineWidth',2,'Color',mixerColours(m,:));
end
drawnow();
legend( 'Diode 1 (mean amplitude)',...
        'Diode 1 (quad fit)',...
        'Diode 2 (mean amplitude)',...
        'Diode 2 (quad fit)',...
        'Diode 3 (mean amplitude)',...
        'Diode 3 (quad fit)',...
        'Location','best');
plot([powerInVolts(diodePowerLevelsToFit(1)) powerInVolts(diodePowerLevelsToFit(1))],get(gca,'YLim'),'k');
plot([powerInVolts(diodePowerLevelsToFit(end)) powerInVolts(diodePowerLevelsToFit(end))],get(gca,'YLim'),'k');
title('Diode Output vs. Input Voltage')
xlabel('Input [V]');
ylabel('Output [V]');
grid on;
ylim([0 2]);
xlim([0 3])
if savePlots; savePlot(savePlotDir,'QuadFitDiodeVsVolts'); end;



figure;
for m=1:nMixers
    plot(powerInVolts,meanAmpMixers(m,:),'-o','LineWidth',1,'Color',mixerColours(m,:));
    hold all;
    plot(powerInVolts,linFitMixerPoints(m,:),'--','LineWidth',2,'Color',mixerColours(m,:));
end
drawnow();
legend( 'Mixer 1 (mean amplitude)',...
        'Mixer 1 (linear fit)',...
        'Mixer 2 (mean amplitude)',...
        'Mixer 2 (linear fit)',...
        'Mixer 3 (mean amplitude)',...
        'Mixer 3 (linear fit)',...
        'Location','best');
plot([powerInVolts(mixerPowerLevelsToFit(1)) powerInVolts(mixerPowerLevelsToFit(1))],get(gca,'YLim'),'k');
plot([powerInVolts(mixerPowerLevelsToFit(end)) powerInVolts(mixerPowerLevelsToFit(end))],get(gca,'YLim'),'k');
title('Mixer Output vs. Input Voltage')
xlabel('Input [V]');
ylabel('Output [V]');
grid on;
ylim([0 0.7]);
if savePlots; savePlot(savePlotDir,'LinFitMixerVsVolts'); end;

mixerGoodXLim = powerInVolts(mixerGoodPlotRange);
mixerGoodXLim(2) = mixerGoodXLim(2) + (powerInVolts(mixerGoodPlotRange(2)+1)-powerInVolts(mixerGoodPlotRange(2)))/2;
xlim(mixerGoodXLim);
legend( 'Mixer 1 (mean amplitude)',...
        'Mixer 1 (linear fit)',...
        'Mixer 2 (mean amplitude)',...
        'Mixer 2 (linear fit)',...
        'Mixer 3 (mean amplitude)',...
        'Mixer 3 (linear fit)',...
        'Location','best');
if savePlots; savePlot(savePlotDir,'LinFitMixerVsVolts_zoom'); end;

figure;
for m=1:nMixers
    plot(powerInVolts,meanDiodes(m,:),'o-','LineWidth',1,'Color',mixerColours(m,:));
    hold all;
    plot(powerInVolts,linFitDiodePoints(m,:),'--','LineWidth',2,'Color',mixerColours(m,:));
end
drawnow();
legend( 'Diode 1 (mean amplitude)',...
        'Diode 1 (linear fit)',...
        'Diode 2 (mean amplitude)',...
        'Diode 2 (linear fit)',...
        'Diode 3 (mean amplitude)',...
        'Diode 3 (linear fit)',...
        'Location','best');
plot([powerInVolts(diodePowerLevelsToFit(1)) powerInVolts(diodePowerLevelsToFit(1))],get(gca,'YLim'),'k');
plot([powerInVolts(diodePowerLevelsToFit(end)) powerInVolts(diodePowerLevelsToFit(end))],get(gca,'YLim'),'k');
title('Diode Output vs. Input Voltage')
xlabel('Input [V]');
ylabel('Output [V]');
grid on;
ylim([0 2]);
if savePlots; savePlot(savePlotDir,'LinFitDiodeVsVolts'); end;

diodeGoodXLim = powerInVolts(diodeGoodPlotRange);
diodeGoodXLim(2) = diodeGoodXLim(2) + (powerInVolts(diodeGoodPlotRange(2)+1)-powerInVolts(diodeGoodPlotRange(2)))/2;
xlim(diodeGoodXLim);
drawnow();
legend( 'Diode 1 (mean amplitude)',...
        'Diode 1 (linear fit)',...
        'Diode 2 (mean amplitude)',...
        'Diode 2 (linear fit)',...
        'Diode 3 (mean amplitude)',...
        'Diode 3 (linear fit)',...
        'Location','best');
if savePlots; savePlot(savePlotDir,'LinFitDiodeVsVolts_zoom'); end;


figure;
for m=1:nMixers
    plot(powerInVolts,sqrtMeanDiodes(m,:),'o-','LineWidth',1,'Color',mixerColours(m,:));
    hold all;
    plot(powerInVolts,linFitSqrtDiodePoints(m,:),'--','LineWidth',2,'Color',mixerColours(m,:));
end
drawnow();
legend( 'Diode 1 (mean amplitude)',...
        'Diode 1 (linear fit)',...
        'Diode 2 (mean amplitude)',...
        'Diode 2 (linear fit)',...
        'Diode 3 (mean amplitude)',...
        'Diode 3 (linear fit)',...
        'Location','best');
plot([powerInVolts(diodePowerLevelsToFit(1)) powerInVolts(diodePowerLevelsToFit(1))],get(gca,'YLim'),'k');
plot([powerInVolts(diodePowerLevelsToFit(end)) powerInVolts(diodePowerLevelsToFit(end))],get(gca,'YLim'),'k');
title('sqrt(Diode) vs. Input Voltage')
xlabel('Input [V]');
ylabel('sqrt(Diode) [sqrt(V)]');
grid on;
ylim([0 2]);
if savePlots; savePlot(savePlotDir,'LinFitSqrtDiodeVsVolts'); end;

xlim(diodeGoodXLim);
drawnow();
legend( 'Diode 1 (mean amplitude)',...
        'Diode 1 (linear fit)',...
        'Diode 2 (mean amplitude)',...
        'Diode 2 (linear fit)',...
        'Diode 3 (mean amplitude)',...
        'Diode 3 (linear fit)',...
        'Location','best');
if savePlots; savePlot(savePlotDir,'LinFitSqrtDiodeVsVolts_zoom'); end;


figure;
for m=1:nMixers
    plot(meanDiodes(m,:),meanAmpMixers(m,:),'o-','LineWidth',1,'Color',mixerColours(m,:));
    hold all;
    plot(meanDiodes(m,:),linFitMixerVsDiodePoints(m,:),'--','LineWidth',2,'Color',mixerColours(m,:));
end
drawnow();
legend( 'Mixer 1 (mean amplitude)',...
        'Mixer 1 (linear fit)',...
        'Mixer 2 (mean amplitude)',...
        'Mixer 2 (linear fit)',...
        'Mixer 3 (mean amplitude)',...
        'Mixer 3 (linear fit)',...
        'Location','best');
% plot([powerInVolts(diodePowerLevelsToFit(1)) powerInVolts(diodePowerLevelsToFit(1))],get(gca,'YLim'),'k');
% plot([powerInVolts(diodePowerLevelsToFit(end)) powerInVolts(diodePowerLevelsToFit(end))],get(gca,'YLim'),'k');
title('Max Mixer Output vs. Diode Output')
xlabel('Diode [V]');
ylabel('Mixer [V]');
grid on;
if savePlots; savePlot(savePlotDir,'LinFitMixerVsDiodeVs'); end;

tmpGoodRange = diodeGoodPlotRange;
if (tmpGoodRange(1)>1)
    tmpGoodRange(1) = tmpGoodRange(1)-1;
end
if (tmpGoodRange(2)<nPowerLevels)
    tmpGoodRange(2) = tmpGoodRange(2)+1;
end
tmpXLim = meanDiodes(:,tmpGoodRange);
diodeGoodXLim = [min(tmpXLim(:)) max(tmpXLim(:))];
xlim(diodeGoodXLim);
drawnow();
if savePlots; savePlot(savePlotDir,'LinFitMixerVsDiode_zoom'); end;

% figure;
% for m=1:nMixers
%     semilogy(powerLevels,meanAmpMixers(m,:),'-','LineWidth',2,'Color',mixerColours(m,:));
%     hold all;
%     semilogy(powerLevels,linFitMixerAmp(m,:),'--','Color',mixerColours(m,:));
% end
% legend( 'Mixer 1 (mean amplitude)',...
%         'Mixer 1 (linear fit)',...
%         'Mixer 2 (mean amplitude)',...
%         'Mixer 2 (linear fit)',...
%         'Mixer 3 (mean amplitude)',...
%         'Mixer 3 (linear fit)',...
%         'Location','best');
% title('Mixer Output vs. Input Power (Log Axis)')
% xlabel('Input Power [dBm]');
% ylabel('Output [V]');
% grid on;
% ylim([0 1.3]);
% if savePlots; savePlot(savePlotDir,'LinFitMixerVsPower_Log'); end;

%% attempt to reconstruct phase from fits

% % fitted log(mixer amplitude) vs. power coefficients
% a = linFitLogMixer(:,1);
% b = linFitLogMixer(:,2);
% 
% % fitted offset of sin fit on mixer
% %mixerPowerLevelsToUse = 1:6;
% %c = nanmean(mixerCalConsts(:,mixerPowerLevelsToUse,4),2);
% c = mixerCalConsts(:,:,4);
% 
% % fitted sqrt(Diode) amplitude vs. power coefficients
% d = linFitSqrtDiode(:,1);
% e = linFitSqrtDiode(:,2);
% 
% % just use mean diodes as crude way to ignore cross-talk
% phases_MeanDiode = NaN(nMixers,nPowerLevels,nSamples);
% for m=1:nMixers
%     for p=1:nPowerLevels
%         tmpMixer = mixers(m,p,:);
%         tmpDiode = meanDiodes(m,p);
%         tmpA = a(m);
%         tmpB = b(m);
%         tmpC = c(m,p);
%         tmpD = d(m);
%         tmpE = e(m);
%         phases_MeanDiode(m,p,:) = getPhaseNewNormalisation(tmpMixer,tmpDiode,tmpA,tmpB,tmpC,tmpD,tmpE);
%     end
% end
% 
% powerLevelsToUse = [1:2 7:nPowerLevels];
% %powerLevelsToUse = 3:6;
% for m=1:nMixers
%     figure;
%     legLabels = cell(1,nPowerLevels);
%     for p=powerLevelsToUse
%         plot(squeeze(phases_MeanDiode(m,p,:)));
%         hold all;
%         legLabels{p} = sprintf('%d dBm',powerLevels(p));
%     end
%     title(sprintf('Reconstructed Mixer %d Phase',m))
%     xlabel('Sample No.');
%     ylabel('Phase [degrees]');
%     grid on;
%     legend(legLabels(powerLevelsToUse));
%     if savePlots; savePlot([saveDir '/ReconstructedPhase/InSaturation'],sprintf('PhaseMixer%d',m)); end;
% end
% 
% % using complete diode (with cross-talk) but no attempt to remove
% % the cross-talk -> see effect on phase.
% phases_FullDiode = NaN(nMixers,nPowerLevels,nSamples);
% for m=1:nMixers
%     for p=1:nPowerLevels
%         tmpMixer = mixers(m,p,:);
%         tmpDiode = diodes(m,p,:);
%         tmpA = a(m);
%         tmpB = b(m);
%         tmpC = c(m,p);
%         tmpD = d(m);
%         tmpE = e(m);
%         phases_FullDiode(m,p,:) = getPhaseNewNormalisation(tmpMixer,tmpDiode,tmpA,tmpB,tmpC,tmpD,tmpE);
%     end
% end

%% save results to a structure;
if savePlots; save([savePlotDir '/processedData.mat']); end;