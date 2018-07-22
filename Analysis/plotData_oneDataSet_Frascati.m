%% phase plots
if (displayAll)
    figure;
end
plot(timeAxis,meanPhaseAlongPulse(1,:),'LineWidth',2);
hold all;
plot(timeAxis,meanPhaseAlongPulse(2,:),'LineWidth',2);
plot(timeAxis,meanPhaseAlongPulse(3,:),'LineWidth',2);
title([dataDescr 'Mean Phase Along Pulse']);
yLim = get(gca,'YLim');
plot([timeAxis(min(sampleRange)) timeAxis(min(sampleRange))],yLim,'k');
plot([timeAxis(max(sampleRange)) timeAxis(max(sampleRange))],yLim,'k');
xlim([timeAxis(plotStartSamp) timeAxis(plotEndSamp)])
xlabel('Time [ns]');
ylabel('Phase [degrees]')
legend(mixerLabels,'Location','best');
format_plots;
if (savePlots)
    saveName = [saveDir '/meanPhaseAlongPulse'];
    print([saveName '.png'],'-dpng');
    savefig([saveName '.fig']);
end
hold off;

if (displayAll)
    figure;
end
plot(timeAxis,stdPhaseAlongPulse(1,:),'LineWidth',2);
hold all;
plot(timeAxis,stdPhaseAlongPulse(2,:),'LineWidth',2);
plot(timeAxis,stdPhaseAlongPulse(3,:),'LineWidth',2);
title([dataDescr 'Std Phase Along Pulse']);
xlim([timeAxis(plotStartSamp) timeAxis(plotEndSamp)])
legLabels = mixerLabels;
for mon=1:nMons
    legLabels{mon} = [mixerLabels{mon} sprintf(' (%.2f^o)',meanStdPhaseAlongPulse(mon))];
end
legend(legLabels,'Location','best');
yLim = get(gca,'YLim');
plot([timeAxis(min(sampleRange)) timeAxis(min(sampleRange))],yLim,'k');
plot([timeAxis(max(sampleRange)) timeAxis(max(sampleRange))],yLim,'k');
xlabel('Time [ns]');
ylabel('Phase [degrees]')
format_plots;
if (savePlots)
    saveName = [saveDir '/stdPhaseAlongPulse'];
    print([saveName '.png'],'-dpng');
    savefig([saveName '.fig']);
end
hold off;

if (displayAll)
    figure;
end
plot(meanPulsePhase(1,:),'LineWidth',2);
hold all;
plot(meanPulsePhase(2,:),'LineWidth',2);
plot(meanPulsePhase(3,:),'LineWidth',2);
title([dataDescr 'Mean Pulse Phase']);
legLabels = mixerLabels;
for mon=1:nMons
    legLabels{mon} = [mixerLabels{mon} sprintf(' (std = %.2f^o)',stdMeanPulsePhase(mon))];
end
format_plots;
legend(legLabels,'Location','best');
if (savePlots)
    saveName = [saveDir '/meanPulsePhase'];
    print([saveName '.png'],'-dpng');
    savefig([saveName '.fig']);
end
hold off;

%% frascati phase - correlation mean phase

if (displayAll)
    a = figure;
end
try
    plot(meanPulsePhase(1,:),meanPulsePhase(2,:),'o','MarkerSize',7);
    title(sprintf('%s%s vs. %s (corr=%.2f, grad=%.2f)',dataDescr,mixerLabels{1},mixerLabels{2},corrMeanMix1Mix2,fitMeanMix1Mix2(1)));
    hold all;
    xVals = min(meanPulsePhase(1,:)):0.1:max(meanPulsePhase(1,:));
    yVals = fitMeanMix1Mix2(1).*xVals + fitMeanMix1Mix2(2);
    plot(xVals,yVals,'k','LineWidth',2);
    xlabel('Mix 1 Phase [degrees]')
    ylabel('Mix 2 Phase [degrees]')
    format_plots;
    if (savePlots)
        saveName = [saveDir '/mix1Vsmix2'];
        print([saveName '.png'],'-dpng');
        savefig([saveName '.fig']);
    end
catch
    if (displayAll)
        close(a);
    end    
end
hold off;

if (displayAll)
    a = figure;
end
try
    plot(meanPulsePhase(1,:),meanPulsePhase(3,:),'o','MarkerSize',7);
    title(sprintf('%s vs. %s (corr=%.2f, grad=%.2f)',mixerLabels{1},mixerLabels{3},corrMeanMix1Mix3,fitMeanMix1Mix3(1)));
    hold all;
    xVals = min(meanPulsePhase(1,:)):0.1:max(meanPulsePhase(1,:));
    yVals = fitMeanMix1Mix3(1).*xVals + fitMeanMix1Mix3(2);
    plot(xVals,yVals,'k','LineWidth',2);
    xlabel('Mix 1 Phase [degrees]')
    ylabel('Mix 3 Phase [degrees]')
    format_plots;
    if (savePlots)
        saveName = [saveDir '/mix1vsmix3'];
        print([saveName '.png'],'-dpng');
        savefig([saveName '.fig']);
    end
catch
    if (displayAll)
        close(a);
    end  
end
hold off;

if (displayAll)
    a = figure;
end
try
    plot(meanPulsePhase(2,:),meanPulsePhase(3,:),'o','MarkerSize',7);
    title(sprintf('%s vs. %s (corr=%.2f, grad=%.2f)',mixerLabels{2},mixerLabels{3},corrMeanMix2Mix3,fitMeanMix2Mix3(1)));
    hold all;
    xVals = min(meanPulsePhase(2,:)):0.1:max(meanPulsePhase(2,:));
    yVals = fitMeanMix2Mix3(1).*xVals + fitMeanMix2Mix3(2);
    plot(xVals,yVals,'k','LineWidth',2);
    xlabel('Mix 2 Phase [degrees]')
    ylabel('Mix 3 Phase [degrees]')
    format_plots;
    if (savePlots)
        saveName = [saveDir '/mix2vsmix3'];
        print([saveName '.png'],'-dpng');
        savefig([saveName '.fig']);
    end
catch
    if (displayAll)
        close(a);
    end  
end
hold off;

if (displayAll)
    figure;
end
plot(timeAxis,resolution12,'LineWidth',2);
hold all;
plot(timeAxis,resolution13,'LineWidth',2);
plot(timeAxis,resolution23,'LineWidth',2);
xlim([timeAxis(plotStartSamp) timeAxis(plotEndSamp)])
yLim = get(gca,'YLim');
plot([timeAxis(min(sampleRange)) timeAxis(min(sampleRange))],yLim,'k');
plot([timeAxis(max(sampleRange)) timeAxis(max(sampleRange))],yLim,'k');
legLabels = mixerLabels;
legLabels{1} = sprintf('%s, %s (%.2f^o)',mixerLabels{1},mixerLabels{2},meanResolution12);
legLabels{2} = sprintf('%s, %s (%.2f^o)',mixerLabels{1},mixerLabels{3},meanResolution13);
legLabels{3} = sprintf('%s, %s (%.2f^o)',mixerLabels{2},mixerLabels{3},meanResolution23);
format_plots;
legend(legLabels,'Location','best');
title([dataDescr 'Resolution']);
xlabel('Time [ns]');
ylabel('Resolution [degrees]');
if (savePlots)
    saveName = [saveDir '/resolution'];
    print([saveName '.png'],'-dpng');
    savefig([saveName '.fig']);
end
hold off;

%% frascati phase - correlation vs. sample no.
if (displayAll)
    a = figure;
end
try
    plot(timeAxis,corrSamplesMix1Mix2,'LineWidth',2);
    hold all;
    plot(timeAxis,corrSamplesMix1Mix3,'LineWidth',2);
    plot(timeAxis,corrSamplesMix2Mix3,'LineWidth',2);
    
    title(sprintf('%sCorrelation vs. Sample No.',dataDescr));
    xlabel('Time along pulse [ns]')
    ylabel('Correlation')
    
    xlim([timeAxis(plotStartSamp) timeAxis(plotEndSamp)])

    format_plots;
    
    legLabels = mixerLabels;
    legLabels{1} = sprintf('%s, %s (%.2f)',mixerLabels{1},mixerLabels{2},meanCorrSamplesMix1Mix2);
    legLabels{2} = sprintf('%s, %s (%.2f)',mixerLabels{1},mixerLabels{3},meanCorrSamplesMix1Mix3);
    legLabels{3} = sprintf('%s, %s (%.2f)',mixerLabels{2},mixerLabels{3},meanCorrSamplesMix2Mix3);
    legend(legLabels,'Location','best');
    
    yLim = get(gca,'YLim');
    if (yLim(1) < -1)
        yLim(1) = -1;
    end
    if (yLim(2) > 1)
        yLim(2) = 1;
    end
    ylim(yLim);
    
    plot([timeAxis(min(sampleRange)) timeAxis(min(sampleRange))],yLim,'k');
    plot([timeAxis(max(sampleRange)) timeAxis(max(sampleRange))],yLim,'k');

    if (savePlots)
        saveName = [saveDir '/corrSamples'];
        print([saveName '.png'],'-dpng');
        savefig([saveName '.fig']);
    end
catch
    if (displayAll)
        close(a);
    end    
end
hold off;

%% frascati phase - correlation pulse shape vs. time
if (displayAll)
    a = figure;
end
try
    plot(corrShapeMix1Mix2,'LineWidth',2);
    hold all;
    plot(corrShapeMix1Mix3,'LineWidth',2);
    plot(corrShapeMix2Mix3,'LineWidth',2);
    
    title(sprintf('%sCorrelation Pulse Shape vs. Time',dataDescr));
    xlabel('Time [Pulse no.]')
    ylabel('Correlation')
    
    format_plots;
    
    legLabels = mixerLabels;
    legLabels{1} = sprintf('%s, %s (%.2f)',mixerLabels{1},mixerLabels{2},meanCorrShapeMix1Mix2);
    legLabels{2} = sprintf('%s, %s (%.2f)',mixerLabels{1},mixerLabels{3},meanCorrShapeMix1Mix3);
    legLabels{3} = sprintf('%s, %s (%.2f)',mixerLabels{2},mixerLabels{3},meanCorrShapeMix2Mix3);
    legend(legLabels,'Location','best');
    
    if (savePlots)
        saveName = [saveDir '/corrShape'];
        print([saveName '.png'],'-dpng');
        savefig([saveName '.fig']);
    end
catch
    if (displayAll)
        close(a);
    end    
end
hold off;

%% baseline noise plots
if (displayAll)
    figure;
end
if (isFONTData)
    plot(timeAxis(baselineStartSamples),meanSampleMixBaseSt','LineWidth',2)
    ylabel('Output [counts]');
else
    plot(timeAxis(baselineStartSamples),1000*meanSampleMixBaseSt','LineWidth',2)
    ylabel('Output [mV]');
end
legend(mixerLabels,'Location','best');
title([dataDescr 'Mean Mixer Before Pulse']);
xlabel('Time [ns]')
format_plots;
if (savePlots)
    saveName = [saveDir '/meanSampleMixBaseSt'];
    print([saveName '.png'],'-dpng');
    savefig([saveName '.fig']);
end
hold off;

if (displayAll)
    figure;
end
if (isFONTData)
    plot(timeAxis(baselineStartSamples),stdSampleMixBaseSt','LineWidth',2)
    ylabel('Jitter [counts]');
else
    plot(timeAxis(baselineStartSamples),1000*stdSampleMixBaseSt','LineWidth',2)
    ylabel('Jitter [mV]');
end
legend(mixerLabels,'Location','best');
title([dataDescr 'Std Mixer Before Pulse']);
xlabel('Time [ns]')
format_plots;
if (savePlots)
    saveName = [saveDir '/stdSampleMixBaseSt'];
    print([saveName '.png'],'-dpng');
    savefig([saveName '.fig']);
end
hold off;

if (displayAll)
    figure;
end
if (isFONTData)
    plot(timeAxis(baselineEndSamples),meanSampleMixBaseEn','LineWidth',2)
    ylabel('Output [counts]');    
else
    plot(timeAxis(baselineEndSamples),1000*meanSampleMixBaseEn','LineWidth',2)
    ylabel('Output [mV]');    
end
legend(mixerLabels,'Location','best');
title([dataDescr 'Mean Mixer After Pulse']);
xlabel('Time [ns]')
format_plots;
if (savePlots)
    saveName = [saveDir '/meanSampleMixBaseEn'];
    print([saveName '.png'],'-dpng');
    savefig([saveName '.fig']);
end
hold off;

if (displayAll)
    figure;
end
if (isFONTData)
    plot(timeAxis(baselineEndSamples),stdSampleMixBaseEn','LineWidth',2)
    ylabel('Jitter [counts]');    
else
    plot(timeAxis(baselineEndSamples),1000*stdSampleMixBaseEn','LineWidth',2)
    ylabel('Jitter [mV]');
end
legend(mixerLabels,'Location','best');
title([dataDescr 'Std Mixer After Pulse']);
xlabel('Time [ns]')
format_plots;
if (savePlots)
    saveName = [saveDir '/stdSampleMixBaseEn'];
    print([saveName '.png'],'-dpng');
    savefig([saveName '.fig']);
end
hold off;

if (displayAll)
    figure;
end
plot(timeAxis(baselineStartSamples),stdSamplePhBaseSt','LineWidth',2)
legend(mixerLabels,'Location','best');
title([dataDescr 'Std "Phase" Before Pulse']);
xlabel('Time [ns]')
ylabel('Jitter [degrees]');
format_plots;
if (savePlots)
    saveName = [saveDir '/stdSamplePhBaseSt'];
    print([saveName '.png'],'-dpng');
    savefig([saveName '.fig']);
end
hold off;

if (displayAll)
    figure;
end
plot(timeAxis(baselineEndSamples),stdSamplePhBaseEn','LineWidth',2)
legend(mixerLabels,'Location','best');
title([dataDescr 'Std "Phase" After Pulse']);
xlabel('Time [ns]')
ylabel('Jitter [degrees]');
format_plots;
if (savePlots)
    saveName = [saveDir '/stdSamplePhBaseEn'];
    print([saveName '.png'],'-dpng');
    savefig([saveName '.fig']);
end
hold off;