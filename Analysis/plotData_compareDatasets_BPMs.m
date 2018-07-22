return;
%%
if (displayAll)
    a = figure;
end
for ds=1:nDataSets
    %plot(data{ds}.meanBPMH{data{ds}.refBPM2Index},data{ds}.meanPulsePhase(2,:),markerType,'MarkerSize',markerSizeBig,'MarkerEdgeColor',datSetColors(ds,:),'MarkerFaceColor',datSetColors(ds,:))
    hold all;
    plot(data{ds}.meanBPMH{data{ds}.refBPM2Index},data{ds}.meanPulsePhase(3,:),markerType,'MarkerSize',markerSizeBig,'MarkerEdgeColor',datSetColors(ds,:),'MarkerFaceColor',datSetColors(ds,:))
    title('');
%     xlim([-1 9])
%     ylim([-10 25])
    xlabel('CT.608H [mm]');
    ylabel('Mon3 Phase [degrees]');
end
format_compare_legend;
format_plots;
mon3XLim = get(gca,'XLim');
mon3YLim = get(gca,'YLim');
if (savePlots)
    savePlot(saveDir,'Mon3VsCT608');
end
% for ds=1:nDataSets
%     fitCoeffs = nanpolyfit(data{ds}.meanBPMH{data{ds}.refBPM2Index},data{ds}.meanPulsePhase(3,:),2);
%     xVals = mon3XLim(1):0.1:mon3XLim(2);
%     plot(xVals,polyval(fitCoeffs,xVals),'Color',datSetColors(ds,:));
%     %plot(data{ds}.meanBPMH{data{ds}.refBPM2Index},data{ds}.meanPulsePhase(3,:),'.','MarkerEdgeColor',datSetColors(ds,:),'MarkerFaceColor',datSetColors(ds,:))
%     hold all;
% end
% xlim(mon3XLim);
% ylim(mon3YLim);

if (displayAll)
    a = figure;
end
for ds=1:nDataSets
    %plot(data{ds}.meanBPMH{data{ds}.refBPM2Index},data{ds}.meanPulsePhase(2,:),markerType,'MarkerSize',markerSizeBig,'MarkerEdgeColor',datSetColors(ds,:),'MarkerFaceColor',datSetColors(ds,:))
    hold all;
    plot(data{ds}.meanBPMH{data{ds}.refBPM2Index},data{ds}.meanPulsePhase(2,:),markerType,'MarkerSize',markerSizeBig,'MarkerEdgeColor',datSetColors(ds,:),'MarkerFaceColor',datSetColors(ds,:))
    title(sprintf('R56 = %.3f m',dataSetValues(ds)));
    xlim(mon3XLim)
    ylim(mon3YLim)
    xlabel('CT.608H [mm]');
    ylabel('Mon2 Phase [degrees]');
end
format_compare_legend;
format_plots;
if (savePlots)
    savePlot(saveDir,'Mon2VsCT608');
end

for ds=1:nDataSets
    if (displayAll)
        a = figure;
    end
    scatter(data{ds}.meanBPMH{data{ds}.refBPM2Index},data{ds}.meanPulsePhase(2,:),'b')
    hold all;
    scatter(data{ds}.meanBPMH{data{ds}.refBPM12ndex},data{ds}.meanPulsePhase(3,:),'r')
    title(sprintf('R56 = %.3f m',dataSetValues(ds)));
%     xlim([-1 9])
    ylim([-10 25])
    xlabel('CT.608H [mm]');
    ylabel('Phase [degrees]');
    format_plots;
    
%     if (savePlots)
%         savePlot(saveDir,sprintf('PhaseVsCT608_R56_%.3f',dataSetValues(ds)));
%     end

end

for ds=1:nDataSets
    if (displayAll)
        a = figure;
    end
    scatter(data{ds}.meanBPMH{data{ds}.refBPM1Index},data{ds}.meanPulsePhase(2,:),'b')
    hold all;
    scatter(data{ds}.meanBPMH{data{ds}.refBPM1Index},data{ds}.meanPulsePhase(3,:),'r')
    title(sprintf('R56 = %.3f m',dataSetValues(ds)));
%     xlim([-1 9])
    ylim([-10 25])
    xlabel('CL.502H [mm]');
    ylabel('Phase [degrees]');
    format_plots;
    
%     if (savePlots)
%         savePlot(saveDir,sprintf('PhaseVsCL502_R56_%.3f',dataSetValues(ds)));
%     end

end

for ds=1:nDataSets
    if (displayAll)
        a = figure;
    end
    scatter(data{ds}.meanBPMH{data{ds}.refBPM1Index},data{ds}.meanBPMH{data{ds}.refBPM2Index},'b')
    title(sprintf('R56 = %.3f m',dataSetValues(ds)));
%     xlim([-1 9])
%     ylim([-10 25])
    xlabel('CL.502H [mm]');
    ylabel('CT.608H [mm]');
    format_plots;
    
    if (savePlots)
        savePlot(saveDir,sprintf('PhaseVsCL502_R56_%.3f',dataSetValues(ds)));
    end

end

for ds=1:nDataSets
    if (displayAll)
        a = figure;
    end
    scatter(data{ds}.meanBPMH{1},data{ds}.meanBPMH{data{ds}.refBPM2Index},'b')
    title(sprintf('R56 = %.3f m',dataSetValues(ds)));
%     xlim([-1 9])
%     ylim([-10 25])
    xlabel('CL.402H [mm]');
    ylabel('CT.608H [mm]');
    format_plots;
    
    if (savePlots)
        savePlot(saveDir,sprintf('CL402VsCT608_R56_%.3f',dataSetValues(ds)));
    end

end
%%
% corrMon1_BPMH = NaN(nDataSets,data{1}.nBPMs);
% corrMon2_BPMH = NaN(nDataSets,data{1}.nBPMs);
% corrMon3_BPMH = NaN(nDataSets,data{1}.nBPMs);
% corrMon1_BPMS = NaN(nDataSets,data{1}.nBPMs);
% corrMon2_BPMS = NaN(nDataSets,data{1}.nBPMs);
% corrMon3_BPMS = NaN(nDataSets,data{1}.nBPMs);
% for ds=1:nDataSets
%     corrMon1_BPMH(ds,:) = data{ds}.corrMon1_BPMH;
%     corrMon2_BPMH(ds,:) = data{ds}.corrMon2_BPMH;
%     corrMon3_BPMH(ds,:) = data{ds}.corrMon3_BPMH;
%     corrMon1_BPMS(ds,:) = data{ds}.corrMon1_BPMS;
%     corrMon2_BPMS(ds,:) = data{ds}.corrMon2_BPMS;
%     corrMon3_BPMS(ds,:) = data{ds}.corrMon3_BPMS;
% end
% 
% % temprorarily added to look at position scan data
% ct430H = NaN(nDataSets,1);
% ct430S = NaN(nDataSets,1);
% ct430V = NaN(nDataSets,1);
% for ds=1:nDataSets
%     ct430H(ds) = nanmean(data{ds}.meanBPMH{10});
%     ct430S(ds) = nanmean(data{ds}.meanBPMS{10});    
%     ct430V(ds) = nanmean(data{ds}.meanBPMV{10});
% end
% 
% % mean sample jitter vs. data set
% if (displayAll)
%     a = figure;
% end
% plot(dataSetValues,ct430H,markerType,'MarkerSize',markerSizeBig);
% hold all;
% xlabel(dataSetValueLabel);
% ylabel('Position [mm]')
% title([dataDescr 'CT.430 H']);
% format_plots;
% if (savePlots)
%     saveName = [saveDir '/ct430H'];
%     print([saveName '.png'],'-dpng');
%     savefig([saveName '.fig']);
% end
% hold off;
% 
% % mean sample jitter vs. data set
% if (displayAll)
%     a = figure;
% end
% plot(dataSetValues,ct430S,markerType,'MarkerSize',markerSizeBig);
% hold all;
% xlabel(dataSetValueLabel);
% ylabel('Current [A]')
% title([dataDescr 'CT.430 S']);
% format_plots;
% if (savePlots)
%     saveName = [saveDir '/ct430S'];
%     print([saveName '.png'],'-dpng');
%     savefig([saveName '.fig']);
% end
% hold off;
% 
% % mean sample jitter vs. data set
% if (displayAll)
%     a = figure;
% end
% plot(dataSetValues,ct430V,markerType,'MarkerSize',markerSizeBig);
% hold all;
% xlabel(dataSetValueLabel);
% ylabel('Position [mm]')
% title([dataDescr 'CT.430 V']);
% format_plots;
% if (savePlots)
%     saveName = [saveDir '/ct430V'];
%     print([saveName '.png'],'-dpng');
%     savefig([saveName '.fig']);
% end
% hold off;

%%
% bpmName = strrep(data{1}.bpmName,'_','.');
% 
% % correlations with BPMH
% corrMeanMix1BPMH = NaN(1,nDataSets);
% corrMeanMix2BPMH = NaN(1,nDataSets);
% corrMeanMix3BPMH = NaN(1,nDataSets);
% for ds=1:nDataSets
%     meanStdPhaseAlongPulse(ds,:) = data{ds}.meanStdPhaseAlongPulse;
%     corrMeanMix1BPMH(ds) = data{ds}.corrMeanMix1BPMH;
%     corrMeanMix2BPMH(ds) = data{ds}.corrMeanMix2BPMH;
%     corrMeanMix3BPMH(ds) = data{ds}.corrMeanMix3BPMH;
% end
% 
% if (displayAll)
%     a = figure;
% end
% plot(dataSetValues, corrMeanMix1BPMH,markerType,'MarkerSize',markerSizeBig);
% hold all;
% plot(dataSetValues, corrMeanMix2BPMH,markerType,'MarkerSize',markerSizeBig);
% plot(dataSetValues, corrMeanMix3BPMH,markerType,'MarkerSize',markerSizeBig);
% 
% xlabel(dataSetValueLabel);
% ylabel('Correlation with BPMH')
% title([dataDescr 'Phase Correlation with ' bpmName 'H']);
% legend(mixerLabels);
% format_plots;
% if (savePlots)
%     saveName = [saveDir '/corrBPMH'];
%     print([saveName '.png'],'-dpng');
%     savefig([saveName '.fig']);
% end
% hold off;
% 
% % correlations with BPMS
% corrMeanMix1BPMS = NaN(1,nDataSets);
% corrMeanMix2BPMS = NaN(1,nDataSets);
% corrMeanMix3BPMS = NaN(1,nDataSets);
% for ds=1:nDataSets
%     meanStdPhaseAlongPulse(ds,:) = data{ds}.meanStdPhaseAlongPulse;
%     corrMeanMix1BPMS(ds) = data{ds}.corrMeanMix1BPMS;
%     corrMeanMix2BPMS(ds) = data{ds}.corrMeanMix2BPMS;
%     corrMeanMix3BPMS(ds) = data{ds}.corrMeanMix3BPMS;
% end
% 
% if (displayAll)
%     a = figure;
% end
% plot(dataSetValues, corrMeanMix1BPMS,markerType,'MarkerSize',markerSizeBig);
% hold all;
% plot(dataSetValues, corrMeanMix2BPMS,markerType,'MarkerSize',markerSizeBig);
% plot(dataSetValues, corrMeanMix3BPMS,markerType,'MarkerSize',markerSizeBig);
% 
% xlabel(dataSetValueLabel);
% ylabel('Correlation with BPMS')
% title([dataDescr 'Phase Correlation with ' bpmName 'S']);
% legend(mixerLabels);
% format_plots;
% if (savePlots)
%     saveName = [saveDir '/corrBPMS'];
%     print([saveName '.png'],'-dpng');
%     savefig([saveName '.fig']);
% end
% hold off;
% 
% % correlations with BPMV
% corrMeanMix1BPMV = NaN(1,nDataSets);
% corrMeanMix2BPMV = NaN(1,nDataSets);
% corrMeanMix3BPMV = NaN(1,nDataSets);
% for ds=1:nDataSets
%     meanStdPhaseAlongPulse(ds,:) = data{ds}.meanStdPhaseAlongPulse;
%     corrMeanMix1BPMV(ds) = data{ds}.corrMeanMix1BPMV;
%     corrMeanMix2BPMV(ds) = data{ds}.corrMeanMix2BPMV;
%     corrMeanMix3BPMV(ds) = data{ds}.corrMeanMix3BPMV;
% end
% 
% if (displayAll)
%     a = figure;
% end
% plot(dataSetValues, corrMeanMix1BPMV,markerType,'MarkerSize',markerSizeBig);
% hold all;
% plot(dataSetValues, corrMeanMix2BPMV,markerType,'MarkerSize',markerSizeBig);
% plot(dataSetValues, corrMeanMix3BPMV,markerType,'MarkerSize',markerSizeBig);
% 
% xlabel(dataSetValueLabel);
% ylabel('Correlation with BPMV')
% title([dataDescr 'Phase Correlation with ' bpmName 'V']);
% legend(mixerLabels);
% format_plots;
% if (savePlots)
%     saveName = [saveDir '/corrBPMV'];
%     print([saveName '.png'],'-dpng');
%     savefig([saveName '.fig']);
% end
% hold off;
