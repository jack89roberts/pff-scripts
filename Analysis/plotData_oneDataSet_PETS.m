%% phase plots
if (displayAll)
    figure;
end
plot(timeAxis,meanPhaseAlongPulse(1,:),'LineWidth',2);
hold all;
plot(timeAxis,meanPhaseAlongPulse(2,:),'LineWidth',2);
plot(timeAxis,meanPhaseAlongPulse(3,:),'LineWidth',2);
plot(petsTimeAxis,petsMeanPhaseAlongPulse,'LineWidth',2);
title([dataDescr 'Mean Phase Along Pulse']);
yLim = get(gca,'YLim');
plot([timeAxis(min(sampleRange)) timeAxis(min(sampleRange))],yLim,'k');
plot([timeAxis(max(sampleRange)) timeAxis(max(sampleRange))],yLim,'k');
xlim([timeAxis(plotStartSamp) timeAxis(plotEndSamp)]);
xlabel('Time [ns]');
ylabel('Phase [degrees]')
legend([mixerLabels 'PETS'],'Location','best');
format_plots;
if (savePlots)
    saveName = [saveDir '/petsMeanPhaseAlongPulse'];
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
plot(petsTimeAxis,petsStdPhaseAlongPulse,'LineWidth',2);
title([dataDescr 'Std Phase Along Pulse']);
xlim([timeAxis(plotStartSamp) timeAxis(plotEndSamp)])
legLabels = mixerLabels;
for mon=1:nMons
    legLabels{mon} = [mixerLabels{mon} sprintf(' (%.2f^o)',meanStdPhaseAlongPulse(mon))];
end
legLabels{4} = ['PETS' sprintf(' (%.2f^o)',petsMeanStdPhaseAlongPulse)];
legend(legLabels,'Location','best');
yLim = get(gca,'YLim');
plot([timeAxis(min(sampleRange)) timeAxis(min(sampleRange))],yLim,'k');
plot([timeAxis(max(sampleRange)) timeAxis(max(sampleRange))],yLim,'k');
xlabel('Time [ns]');
ylabel('Phase [degrees]')
format_plots;
if (savePlots)
    saveName = [saveDir '/petsStdPhaseAlongPulse'];
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
plot(petsMeanPulsePhase,'LineWidth',2);
title([dataDescr 'Mean Pulse Phase']);
legLabels = mixerLabels;
for mon=1:nMons
    legLabels{mon} = [mixerLabels{mon} sprintf(' (std = %.2f^o)',stdMeanPulsePhase(mon))];
end
legLabels{4} = ['PETS' sprintf(' (std = %.2f^o)',petsStdMeanPulsePhase)];
format_plots;
legend(legLabels,'Location','best');
if (savePlots)
    saveName = [saveDir '/petsMeanPulsePhase'];
    print([saveName '.png'],'-dpng');
    savefig([saveName '.fig']);
end
hold off;

%% pets phase correlations
if (displayAll)
    a = figure;
end
try
    plot(meanPulsePhase(3,:),petsMeanPulsePhase,'o','MarkerSize',7);
    title(sprintf('%sPETS vs. %s (corr=%.2f, grad=%.2f)',dataDescr,mixerLabels{3},corrMeanMix3PETS,fitMeanMix3PETS(1)));
    hold all;
    xVals = min(meanPulsePhase(3,:)):0.1:max(meanPulsePhase(3,:));
    yVals = fitMeanMix3PETS(1).*xVals + fitMeanMix3PETS(2);
    plot(xVals,yVals,'k','LineWidth',2);
    xlabel('Mon 3 Phase [degrees]')
    ylabel('PETS Phase [degrees]')
    format_plots;
    if (savePlots)
        saveName = [saveDir '/mix3VsPETS'];
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
    plot(meanPulsePhase(2,:),petsMeanPulsePhase,'o','MarkerSize',7);
    title(sprintf('%sPETS vs. %s (corr=%.2f, grad=%.2f)',dataDescr,mixerLabels{2},corrMeanMix2PETS,fitMeanMix2PETS(1)));
    hold all;
    xVals = min(meanPulsePhase(2,:)):0.1:max(meanPulsePhase(2,:));
    yVals = fitMeanMix2PETS(1).*xVals + fitMeanMix2PETS(2);
    plot(xVals,yVals,'k','LineWidth',2);
    xlabel('Mon 2 Phase [degrees]')
    ylabel('PETS Phase [degrees]')
    format_plots;
    if (savePlots)
        saveName = [saveDir '/mix2VsPETS'];
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
    plot(meanPulsePhase(1,:),petsMeanPulsePhase,'o','MarkerSize',7);
    title(sprintf('%sPETS vs. %s (corr=%.2f, grad=%.2f)',dataDescr,mixerLabels{1},corrMeanMix1PETS,fitMeanMix1PETS(1)));
    hold all;
    xVals = min(meanPulsePhase(1,:)):0.1:max(meanPulsePhase(1,:));
    yVals = fitMeanMix1PETS(1).*xVals + fitMeanMix1PETS(2);
    plot(xVals,yVals,'k','LineWidth',2);
    xlabel('Mon 1 Phase [degrees]')
    ylabel('PETS Phase [degrees]')
    format_plots;
    if (savePlots)
        saveName = [saveDir '/mix1VsPETS'];
        print([saveName '.png'],'-dpng');
        savefig([saveName '.fig']);
    end
catch
    if (displayAll)
        close(a);
    end    
end
hold off;