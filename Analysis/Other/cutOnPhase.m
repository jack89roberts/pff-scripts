% Makes some plots with pulses split in to bins based on mean upstream
% phase.

%dataSetName= '20150707_2037_FFGain_63_Interleaved_Odd';
%dataSetName= '20150707_2138_FFGain_43_Interleaved_Odd';
%dataSetName= '20150707_2052_FFGain_53_Interleaved_Odd';
%dataSetName= '20150707_2102_FFGain_43_Interleaved_Even';
dataSetName= '20150715_1414_FF_Gain63_Gate185_350_Int_Even';

dataDir = '/home/jack/PhaseFeedforward/Analysis/201507';
saveDir = '/home/jack/PhaseFeedforward/Analysis/201507/CutOnPhase';
savePlots = true;

cutStds = 0.25:0.25:3;
cutOffset = 0;

plotXRange = [550 700];%[290 335];%[550 700];%[sampleRange(1)-20 sampleRange(end)+20]

%% load data
addpath('../');
addpath('../../');

savePlotDir = [saveDir '/' dataSetName];

load([dataDir '/' dataSetName '.mat']);
% sampleRange = 619:628;

format_plots;
close all;
%%
acceptRangeLow =  cutOffset-(stdMeanPulsePhase(2,:).*cutStds);
acceptRangeHigh = cutOffset+(stdMeanPulsePhase(2,:).*cutStds);

nWindows = length(acceptRangeLow);
binColours = varycolor(nWindows);

for mon=1:nMons
    phases(mon,:,:) = phases(mon,:,:)-subtractPhase(mon);
end
[~, phaseFlatness, ~, phaseFlatnessErr ] = nanMeanStdErr(phases(:,:,sampleRange),3);
[meanFlatness,~,meanFlatnessErr,~] = nanMeanStdErr(phaseFlatness,2);

downHist = figure;
hist(meanPulsePhase(3,:,:),20);
title('Distribution Downstream Phases')
xlabel('Phase [degrees]')
ylabel('No. Pulses');
format_plots;
downXLim = get(gca,'XLim');

upHist = figure;
hist(meanPulsePhase(2,:,:),20);
title('Distribution Upstream Phases')
xlabel('Phase [degrees]')
ylabel('No. Pulses');
format_plots;
upXLim = get(gca,'XLim');

histXLim = [min([downXLim(1) upXLim(1)]) max([downXLim(2) upXLim(2)])];

figure(downHist);
xlim(histXLim);
if (savePlots); savePlot(savePlotDir,'histDownstreamPhase'); end;

figure(upHist);
xlim(histXLim);
if (savePlots); savePlot(savePlotDir,'histUpstreamPhase'); end;

nPulsesPerBin = NaN(1,nWindows);

binnedMeanPhase = cell(1,nWindows);
binnedMeanPhaseAlong = cell(1,nWindows);
binnedFlatness = cell(1,nWindows);
binnedStdPhaseAlong = cell(1,nWindows);
binnedMeanPhaseErr = cell(1,nWindows);
binnedMeanPhaseAlongErr = cell(1,nWindows);
binnedFlatnessErr = cell(1,nWindows);
binnedStdPhaseAlongErr = cell(1,nWindows);

meanBinnedFlatness = NaN(nWindows,nMons);
meanBinnedMeanPhase = NaN(nWindows,nMons);
stdBinnedMeanPhase = NaN(nWindows,nMons);
meanBinnedFlatnessErr = NaN(nWindows,nMons);
meanBinnedMeanPhaseErr = NaN(nWindows,nMons);
stdBinnedMeanPhaseErr = NaN(nWindows,nMons);

% figUpPhas = figure;
% figDownPhas = figure;
figStdUpPhas = figure;
figStdDownPhas = figure;
legendLabels = cell(1,nWindows);

for i=1:nWindows

    acceptMon2 = meanPulsePhase(2,:)>=acceptRangeLow(i) & meanPulsePhase(2,:)<acceptRangeHigh(i);
    nPulsesPerBin(i) = sum(acceptMon2);
    binnedPhases = phases(:,acceptMon2,:);

    [binnedMeanPhase{i},~,binnedMeanPhaseErr{i},~] = nanMeanStdErr(binnedPhases(:,:,sampleRange),3);
    [meanBinnedMeanPhase(i,:),stdBinnedMeanPhase(i,:),meanBinnedMeanPhaseErr(i,:),stdBinnedMeanPhaseErr(i,:)] = nanMeanStdErr(binnedMeanPhase{i},2);
    [binnedMeanPhaseAlong{i},binnedStdPhaseAlong{i},binnedMeanPhaseAlongErr{i},binnedStdPhaseAlongErr{i}] = nanMeanStdErr(binnedPhases,2);
    [~,binnedFlatness{i},~,binnedFlatnessErr{i}] = nanMeanStdErr(binnedPhases(:,:,sampleRange),3);    
    [meanBinnedFlatness(i,:),~,meanBinnedFlatnessErr(i,:),~] = nanMeanStdErr(binnedFlatness{i},2);
    
    legendLabels{i} = sprintf('%.2f',cutStds(i));
    
%     figure(figUpPhas)
%     plot(binnedMeanPhaseAlong{i}(2,:),'LineWidth',lineWidthBig,'Color',binColours(i,:));
%     hold all;
%     
%     figure(figDownPhas)
%     plot(binnedMeanPhaseAlong{i}(3,:),'LineWidth',lineWidthBig,'Color',binColours(i,:));
%     hold all;

    figure(figStdUpPhas)
    plot(binnedStdPhaseAlong{i}(2,:),'LineWidth',lineWidthBig,'Color',binColours(i,:));
    hold all;

    figure(figStdDownPhas)
    plot(binnedStdPhaseAlong{i}(3,:),'LineWidth',lineWidthBig,'Color',binColours(i,:));
    hold all;

end

% figure(figDownPhas)
% title('DOWNSTREAM Phase Along Pulse')
% xlabel('Sample No.')
% ylabel('Phase [degrees]')
% xlim(plotXRange);
% format_colorbarLegend;
% format_plots; 
% plotDownPhaseRange = get(gca,'YLim');
% plot([sampleRange(1) sampleRange(1)],plotDownPhaseRange,'k');
% plot([sampleRange(end) sampleRange(end)],plotDownPhaseRange,'k');
% ylim(plotDownPhaseRange);
% if (savePlots); savePlot(savePlotDir, '/downstreamPhase'); end;
% 
% figure(figUpPhas)
% title('UPSTREAM Phase Along Pulse')
% xlabel('Sample No.')
% ylabel('Phase [degrees]')
% xlim(plotXRange);
% format_colorbarLegend;
% format_plots; 
% plot([sampleRange(1) sampleRange(1)],plotDownPhaseRange,'k');
% plot([sampleRange(end) sampleRange(end)],plotDownPhaseRange,'k');
% ylim(plotDownPhaseRange);
% if (savePlots); savePlot(savePlotDir, '/upstreamPhase'); end;

figure(figStdDownPhas)
title('DOWNSTREAM Std Phase Along Pulse')
xlabel('Sample No.')
ylabel('Phase [degrees]')
xlim(plotXRange);
format_colorbarLegend;
format_plots; 
plotDownStdRange = get(gca,'YLim');
plotDownStdRange(1) = 0;
plot([sampleRange(1) sampleRange(1)],plotDownStdRange,'k');
plot([sampleRange(end) sampleRange(end)],plotDownStdRange,'k');
ylim(plotDownStdRange);
if (savePlots); savePlot(savePlotDir, '/downstreamStd'); end;

figure(figStdUpPhas)
title('UPSTREAM Std Phase Along Pulse')
xlabel('Sample No.')
ylabel('Phase [degrees]')
xlim(plotXRange);
format_colorbarLegend;
format_plots;
plot([sampleRange(1) sampleRange(1)],plotDownStdRange,'k');
plot([sampleRange(end) sampleRange(end)],plotDownStdRange,'k');
ylim(plotDownStdRange);
if (savePlots); savePlot(savePlotDir, '/upstreamStd'); end;

% figure;
% herrorbar(meanBinnedMeanPhase(:,2),nPulsesPerBin,meanBinnedMeanPhaseErr(:,2),'o');
% title('No. Pulses Each Bin')
% xlabel('Mean Upstream Bin Phase [degrees]')
% ylabel('No. Pulses')
% format_plots;
% if (savePlots); savePlot(savePlotDir, '/nPulsesPerBin'); end;

% figure;
% herrorbar(meanBinnedMeanPhase(:,2),meanBinnedMeanPhase(:,3),meanBinnedMeanPhaseErr(:,2),'o');
% hold on;
% errorbar(meanBinnedMeanPhase(:,2),meanBinnedMeanPhase(:,3),meanBinnedMeanPhaseErr(:,3),'o');
% title('Mean Phase Downstream vs. Upstream')
% xlabel('Upstream Phase [degrees]')
% ylabel('Downstream Phase [degrees]')
% format_plots;
% if (savePlots); savePlot(savePlotDir, '/downstreamVsUpstreamPhase'); end;

figure;
a = errorbar(cutStds,meanBinnedFlatness(:,2),meanBinnedFlatnessErr(:,2),'bo');
hold on
b = errorbar(cutStds,meanBinnedFlatness(:,3),meanBinnedFlatnessErr(:,3),'ro');
title('Phase Flatness vs. Mean Upstream Phase');
xlabel('Cut on upstream phase [no. sigma]')
ylabel('Flatness [degrees]')
legend([a,b],'Upstream','Downstream')
format_plots;
if (savePlots); savePlot(savePlotDir, '/flatnessVsUpstreamPhase'); end;
