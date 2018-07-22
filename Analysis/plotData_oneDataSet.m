%% frascatiResolutionPlots.m - May 2015
% Uses output from frascatiResolution.m to make and save plots.
% Plots from this script are for results from within one dataset. 
% For plots formatted to compare datasets, use
% frascatiResolutionPlots_compareDatasets.m`
%%
close all;
%%
loadData = true; % if the otuput from frascati resolution that you want to plot is already
                  % in the workspace this can be set to false to plot
                  % without attempting to load any data from files.

if (loadData)
    dataSetName = '20161017_1511_Factor4TBL_-2ACRBend';
end

dataDescr = '';
mixerLabels = {...
    'Mon1',...
    'Mon2',...
    'Mon3'...
};

dataDir = '/home/jack/PhaseFeedforward/Analysis/201610'; %'/home/jack/PhaseFeedforward/Analysis/201504/Feedforward';
saveBaseDir = '/home/jack/PhaseFeedforward/Analysis/201610/Plots';%'/home/jack/PhaseFeedforward/Analysis/201504/Feedforward/Plots';
savePlots = false;
displayAll = true; % if true will make one figure window for every plot. if false will only make one (faster for saving)

plotStartSamp = 1;%500;
plotEndSamp = 768;%900;%768;

%% load data
addpath('../');

if (loadData)
    load([dataDir '/' dataSetName '.mat']);
end
if (savePlots)
    saveDir = [saveBaseDir '/' dataSetName];
    if (~exist(saveDir,'dir'))
        mkdir(saveDir);
    end
end

figure;

%% Frascati phase plots
plotData_oneDataSet_Frascati;

%% PETS
plotData_oneDataSet_PETS;

%% bpm correlation plots
% bpmName = strrep(bpmName,'_','.');
% 
% % correlations with mix1
% if (displayAll)
%     a = figure;
% end
% try
%     plot(meanPulsePhase(1,:),meanBPMH,'o','MarkerSize',7);
%     title(sprintf('%s vs. %sH (corr=%.2f, grad=%.2f)',mixerLabels{1},bpmName,corrMeanMix1BPMH,fitMeanMix1BPMH(1)));
%     hold all;
%     xVals = min(meanPulsePhase(1,:)):0.1:max(meanPulsePhase(1,:));
%     yVals = fitMeanMix1BPMH(1).*xVals + fitMeanMix1BPMH(2);
%     plot(xVals,yVals,'k','LineWidth',2);
%     xlabel('Mix 1 Phase [degrees]')
%     ylabel([bpmName 'H Position [mm]'])
%     format_plots;
%     if (savePlots)
%         saveName = [saveDir '/mix1vs' bpmName 'H'];
%         print([saveName '.png'],'-dpng');
%         savefig([saveName '.fig']);
%     end
% catch
%     if (displayAll)
%         close(a);
%     end  
% end
% hold off;
% 
% if (displayAll)
%     a = figure;
% end
% try
%     plot(meanPulsePhase(1,:),meanBPMS,'o','MarkerSize',7);
%     title(sprintf('%s vs. %sS (corr=%.2f, grad=%.2f)',mixerLabels{1},bpmName,corrMeanMix1BPMS,fitMeanMix1BPMS(1)));
%     hold all;
%     xVals = min(meanPulsePhase(1,:)):0.1:max(meanPulsePhase(1,:));
%     yVals = fitMeanMix1BPMS(1).*xVals + fitMeanMix1BPMS(2);
%     plot(xVals,yVals,'k','LineWidth',2);
%     xlabel('Mix 1 Phase [degrees]')
%     ylabel([bpmName ' Transmission [A]'])
%     format_plots;
%     if (savePlots)
%         saveName = [saveDir '/mix1vs' bpmName 'S'];
%         print([saveName '.png'],'-dpng');
%         savefig([saveName '.fig']);
%     end
% catch
%     if (displayAll)
%         close(a);
%     end  
% end
% hold off;
% 
% if (displayAll)
%     a = figure;
% end
% try
%     plot(meanPulsePhase(1,:),meanBPMV,'o','MarkerSize',7);
%     title(sprintf('%s vs. %sS (corr=%.2f, grad=%.2f)',mixerLabels{1},bpmName,corrMeanMix1BPMV,fitMeanMix1BPMV(1)));
%     hold all;
%     xVals = min(meanPulsePhase(1,:)):0.1:max(meanPulsePhase(1,:));
%     yVals = fitMeanMix1BPMV(1).*xVals + fitMeanMix1BPMV(2);
%     plot(xVals,yVals,'k','LineWidth',2);
%     xlabel('Mix 1 Phase [degrees]')
%     ylabel([bpmName 'V Position [mm]'])
%     format_plots;
%     if (savePlots)
%         saveName = [saveDir '/mix1vs' bpmName 'V'];
%         print([saveName '.png'],'-dpng');
%         savefig([saveName '.fig']);
%     end
% catch
%     if (displayAll)
%         close(a);
%     end  
% end
% hold off;
% 
% % correlations with mix2
% if (displayAll)
%     a = figure;
% end
% try
%     plot(meanPulsePhase(2,:),meanBPMH,'o','MarkerSize',7);
%     title(sprintf('%s vs. %sH (corr=%.2f, grad=%.2f)',mixerLabels{2},bpmName,corrMeanMix2BPMH,fitMeanMix2BPMH(1)));
%     hold all;
%     xVals = min(meanPulsePhase(2,:)):0.1:max(meanPulsePhase(2,:));
%     yVals = fitMeanMix2BPMH(1).*xVals + fitMeanMix2BPMH(2);
%     plot(xVals,yVals,'k','LineWidth',2);
%     xlabel('Mix 2 Phase [degrees]')
%     ylabel([bpmName 'H Position [mm]'])
%     format_plots;
%     if (savePlots)
%         saveName = [saveDir '/mix2vs' bpmName 'H'];
%         print([saveName '.png'],'-dpng');
%         savefig([saveName '.fig']);
%     end
% catch
%     if (displayAll)
%         close(a);
%     end  
% end
% hold off;
% 
% if (displayAll)
%     a = figure;
% end
% try
%     plot(meanPulsePhase(2,:),meanBPMS,'o','MarkerSize',7);
%     title(sprintf('%s vs. %sS (corr=%.2f, grad=%.2f)',mixerLabels{2},bpmName,corrMeanMix2BPMS,fitMeanMix2BPMS(1)));
%     hold all;
%     xVals = min(meanPulsePhase(2,:)):0.1:max(meanPulsePhase(2,:));
%     yVals = fitMeanMix2BPMS(1).*xVals + fitMeanMix2BPMS(2);
%     plot(xVals,yVals,'k','LineWidth',2);
%     xlabel('Mix 2 Phase [degrees]')
%     ylabel([bpmName ' Transmission [A]'])
%     format_plots;
%     if (savePlots)
%         saveName = [saveDir '/mix2vs' bpmName 'S'];
%         print([saveName '.png'],'-dpng');
%         savefig([saveName '.fig']);
%     end
% catch
%     if (displayAll)
%         close(a);
%     end  
% end
% hold off;
% 
% if (displayAll)
%     a = figure;
% end
% try
%     plot(meanPulsePhase(2,:),meanBPMV,'o','MarkerSize',7);
%     title(sprintf('%s vs. %sS (corr=%.2f, grad=%.2f)',mixerLabels{2},bpmName,corrMeanMix2BPMV,fitMeanMix2BPMV(1)));
%     hold all;
%     xVals = min(meanPulsePhase(2,:)):0.1:max(meanPulsePhase(2,:));
%     yVals = fitMeanMix2BPMV(1).*xVals + fitMeanMix2BPMV(2);
%     plot(xVals,yVals,'k','LineWidth',2);
%     xlabel('Mix 2 Phase [degrees]')
%     ylabel([bpmName 'V Position [mm]'])
%     format_plots;
%     if (savePlots)
%         saveName = [saveDir '/mix2vs' bpmName 'V'];
%         print([saveName '.png'],'-dpng');
%         savefig([saveName '.fig']);
%     end
% catch
%     if (displayAll)
%         close(a);
%     end  
% end
% hold off;
% 
% % correlations with mix3
% if (displayAll)
%     a = figure;
% end
% try
%     plot(meanPulsePhase(3,:),meanBPMH,'o','MarkerSize',7);
%     title(sprintf('%s vs. %sH (corr=%.2f, grad=%.2f)',mixerLabels{3},bpmName,corrMeanMix3BPMH,fitMeanMix3BPMH(1)));
%     hold all;
%     xVals = min(meanPulsePhase(3,:)):0.1:max(meanPulsePhase(3,:));
%     yVals = fitMeanMix3BPMH(1).*xVals + fitMeanMix3BPMH(2);
%     plot(xVals,yVals,'k','LineWidth',2);
%     xlabel('Mix 3 Phase [degrees]')
%     ylabel([bpmName 'H Position [mm]'])
%     format_plots;
%     if (savePlots)
%         saveName = [saveDir '/mix3vs' bpmName 'H'];
%         print([saveName '.png'],'-dpng');
%         savefig([saveName '.fig']);
%     end
% catch
%     if (displayAll)
%         close(a);
%     end  
% end
% hold off;
% 
% if (displayAll)
%     a = figure;
% end
% try
%     plot(meanPulsePhase(3,:),meanBPMS,'o','MarkerSize',7);
%     title(sprintf('%s vs. %sS (corr=%.2f, grad=%.2f)',mixerLabels{3},bpmName,corrMeanMix3BPMS,fitMeanMix3BPMS(1)));
%     hold all;
%     xVals = min(meanPulsePhase(3,:)):0.1:max(meanPulsePhase(3,:));
%     yVals = fitMeanMix3BPMS(1).*xVals + fitMeanMix3BPMS(2);
%     plot(xVals,yVals,'k','LineWidth',2);
%     xlabel('Mix 3 Phase [degrees]')
%     ylabel([bpmName ' Transmission [A]'])
%     format_plots;
%     if (savePlots)
%         saveName = [saveDir '/mix3vs' bpmName 'S'];
%         print([saveName '.png'],'-dpng');
%         savefig([saveName '.fig']);
%     end
% catch
%     if (displayAll)
%         close(a);
%     end  
% end
% hold off;
% 
% if (displayAll)
%     a = figure;
% end
% try
%     plot(meanPulsePhase(3,:),meanBPMV,'o','MarkerSize',7);
%     title(sprintf('%s vs. %sS (corr=%.2f, grad=%.2f)',mixerLabels{3},bpmName,corrMeanMix3BPMV,fitMeanMix3BPMV(1)));
%     hold all;
%     xVals = min(meanPulsePhase(3,:)):0.1:max(meanPulsePhase(3,:));
%     yVals = fitMeanMix3BPMV(1).*xVals + fitMeanMix3BPMV(2);
%     plot(xVals,yVals,'k','LineWidth',2);
%     xlabel('Mix 3 Phase [degrees]')
%     ylabel([bpmName 'V Position [mm]'])
%     format_plots;
%     if (savePlots)
%         saveName = [saveDir '/mix3vs' bpmName 'V'];
%         print([saveName '.png'],'-dpng');
%         savefig([saveName '.fig']);
%     end
% catch
%     if (displayAll)
%         close(a);
%     end  
% end
% hold off;
