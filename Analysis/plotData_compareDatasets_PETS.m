% mean phase along pulse
% figure;
% for ds=1:nDataSets
%     plot(data{ds}.timeAxis,data{ds}.meanPhaseAlongPulse(mon,:),'LineWidth',lineWidthBig,'Color',datSetColors(ds,:));
%     hold all;
%     plot(data{ds}.timeAxis,data{ds}.petsMeanPhaseAlongPulse(mon,:),'LineWidth',lineWidthBig,'Color',datSetColors(ds,:));
% end
% title([dataDescr 'Mean Phase Along Pulse, ' mixerLabels{mon}]);
% xlim([data{1}.timeAxis(plotStartSamp) data{1}.timeAxis(plotEndSamp)])
% xlabel('Time [ns]');
% ylabel('Phase [degrees]')
% legLabels = dataSetLabels;
% yLim = get(gca,'YLim');
% plot([data{1}.timeAxis(min(data{1}.sampleRange)) data{1}.timeAxis(min(data{1}.sampleRange))],yLim,'k','LineWidth',lineWidthSmall);
% plot([data{1}.timeAxis(max(data{1}.sampleRange)) data{1}.timeAxis(max(data{1}.sampleRange))],yLim,'k','LineWidth',lineWidthSmall);
% format_compare_legend;
% format_plots;
% %yLim(yLim);
% 
% if (savePlots)
%     saveName = [saveDir '/meanPhaseAlongPulse_' mixerLabels{mon}];
%     print([saveName '.png'],'-dpng');
%     savefig([saveName '.fig']);
% end
% hold off;
% 
% % phase jitter along pulse
% for mon=1:nMons
%     if (plotMon(mon))
%         if (displayAll)
%             figure;
%         end
%         for ds=1:nDataSets
%             plot(data{ds}.timeAxis,data{ds}.stdPhaseAlongPulse(mon,:),'LineWidth',lineWidthBig,'Color',datSetColors(ds,:));
%             hold all;
%         end
%         title([dataDescr 'Std Phase Along Pulse, ' mixerLabels{mon}]);
%         xlim([data{1}.timeAxis(plotStartSamp) data{1}.timeAxis(plotEndSamp)])
%         legLabels = dataSetLabels;
%         for ds=1:nDataSets
%             legLabels{ds} = [dataSetLabels{ds} sprintf(' (%.2f^o)',data{ds}.meanStdPhaseAlongPulse(mon))];
%         end
%         format_compare_legend;
%         yLim = get(gca,'YLim');
%         plot([data{1}.timeAxis(min(data{1}.sampleRange)) data{1}.timeAxis(min(data{1}.sampleRange))],yLim,'k','LineWidth',lineWidthSmall);
%         plot([data{1}.timeAxis(max(data{1}.sampleRange)) data{1}.timeAxis(max(data{1}.sampleRange))],yLim,'k','LineWidth',lineWidthSmall);
%         xlabel('Time [ns]');
%         ylabel('Phase Jitter [degrees]')
%         format_plots;
%         if (savePlots)
%             saveName = [saveDir '/stdPhaseAlongPulse_' mixerLabels{mon}];
%             print([saveName '.png'],'-dpng');
%             savefig([saveName '.fig']);
%         end
%         hold off;
%     end
% end
% 
% % mean sample jitter vs. data set
% meanStdPhaseAlongPulse = NaN(nDataSets,nMons);
% for ds=1:nDataSets
%     meanStdPhaseAlongPulse(ds,:) = data{ds}.meanStdPhaseAlongPulse;
% end
% 
% if (displayAll)
%     a = figure;
% end
% for mon=1:nMons
%     if (plotMon(mon))
%         plot(dataSetValues, meanStdPhaseAlongPulse(:,mon),markerType,'MarkerSize',markerSizeBig);
%         hold all;
%     end
% end
% xlabel(dataSetValueLabel);
% ylabel('Phase Jitter [degrees]')
% title([dataDescr 'Mean Sample Jitter']);
% legend(mixerLabels(boolean(plotMon)));
% format_plots;
% if (savePlots)
%     saveName = [saveDir '/meanSampleStd'];
%     print([saveName '.png'],'-dpng');
%     savefig([saveName '.fig']);
% end
% hold off;
% 
% %%
% % mean diode along pulse
% for mon=1:nMons
%     if (plotMon(mon))
%         if (displayAll)
%             figure;
%         end
%         for ds=1:nDataSets
%             plot(data{ds}.timeAxis,data{ds}.meanDiodeAlongPulse(mon,:),'LineWidth',lineWidthBig,'Color',datSetColors(ds,:));
%             hold all;
%         end
%         title([dataDescr 'Mean Diode Along Pulse, ' mixerLabels{mon}]);
%         xlim([data{1}.timeAxis(plotStartSamp) data{1}.timeAxis(plotEndSamp)])
%         xlabel('Time [ns]');
%         ylabel('Diode [V]')
%         legLabels = dataSetLabels;
%         format_compare_legend;
%         format_plots;
%         yLim = get(gca,'YLim');
%         plot([data{1}.timeAxis(min(data{1}.sampleRange)) data{1}.timeAxis(min(data{1}.sampleRange))],yLim,'k','LineWidth',lineWidthSmall);
%         plot([data{1}.timeAxis(max(data{1}.sampleRange)) data{1}.timeAxis(max(data{1}.sampleRange))],yLim,'k','LineWidth',lineWidthSmall);
% 
%         if (savePlots)
%             saveName = [saveDir '/meanDiodeAlongPulse_' mixerLabels{mon}];
%             print([saveName '.png'],'-dpng');
%             savefig([saveName '.fig']);
%         end
%         hold off;
%     end
% end
% 
% % Diode jitter along pulse
% for mon=1:nMons
%     if (plotMon(mon))
%         if (displayAll)
%             figure;
%         end
%         for ds=1:nDataSets
%             plot(data{ds}.timeAxis,data{ds}.stdDiodeAlongPulse(mon,:),'LineWidth',lineWidthBig,'Color',datSetColors(ds,:));
%             hold all;
%         end
%         title([dataDescr 'Std Diode Along Pulse, ' mixerLabels{mon}]);
%         xlim([data{1}.timeAxis(plotStartSamp) data{1}.timeAxis(plotEndSamp)])
%         legLabels = dataSetLabels;
%         for ds=1:nDataSets
%             legLabels{ds} = [dataSetLabels{ds} sprintf(' (%.2f mV)',data{ds}.meanStdDiodeAlongPulse(mon)*1000)];
%         end
%         format_compare_legend;
%         yLim = get(gca,'YLim');
%         plot([data{1}.timeAxis(min(data{1}.sampleRange)) data{1}.timeAxis(min(data{1}.sampleRange))],yLim,'k','LineWidth',lineWidthSmall);
%         plot([data{1}.timeAxis(max(data{1}.sampleRange)) data{1}.timeAxis(max(data{1}.sampleRange))],yLim,'k','LineWidth',lineWidthSmall);
%         xlabel('Time [ns]');
%         ylabel('Diode Jitter [V]')
%         format_plots;
%         if (savePlots)
%             saveName = [saveDir '/stdDiodeAlongPulse_' mixerLabels{mon}];
%             print([saveName '.png'],'-dpng');
%             savefig([saveName '.fig']);
%         end
%         hold off;
%     end
% end
% 
% % mean diode vs. data set
% meanDiode = NaN(nDataSets,nMons);
% for ds=1:nDataSets
%     meanDiode(ds,:) = abs(nanmean(data{ds}.meanDiode,2));
% end
% 
% if (displayAll)
%     a = figure;
% end
% for mon=1:nMons
%     if (plotMon(mon))
%         plot(dataSetValues, meanDiode(:,mon),markerType,'MarkerSize',markerSizeBig);
%         hold all;
%     end
% end
% xlabel(dataSetValueLabel);
% ylabel('Diode [V]')
% title([dataDescr 'Mean Diode']);
% legend(mixerLabels(boolean(plotMon)));
% format_plots;
% if (savePlots)
%     saveName = [saveDir '/diodeMean'];
%     print([saveName '.png'],'-dpng');
%     savefig([saveName '.fig']);
% end
% hold off;
% 
% %%
% 
% % mean phase vs. time
% for mon=1:nMons
%     if (plotMon(mon))
%         if (displayAll)
%             figure;
%         end
%         for ds=1:nDataSets
%             plot(data{ds}.meanPulsePhase(mon,:),'LineWidth',lineWidthBig,'Color',datSetColors(ds,:));
%             hold all;
%         end
%         title([dataDescr 'Mean Pulse Phase, ' mixerLabels{mon}]);
%         legLabels = dataSetLabels;
%         for ds=1:nDataSets
%             legLabels{ds} = [dataSetLabels{ds} sprintf(' (std = %.2f^o)',data{ds}.stdMeanPulsePhase(mon))];
%         end
%         xlabel('Time [Pulse no.]');
%         ylabel('Phase [degrees]');
%         format_compare_legend;
%         format_plots;
%         if (savePlots)
%             saveName = [saveDir '/meanPulsePhase_' mixerLabels{mon}];
%             print([saveName '.png'],'-dpng');
%             savefig([saveName '.fig']);
%         end
%         hold off;
%     end
% end
% 
% % mean flatness vs. time
% for mon=1:nMons
%     if (plotMon(mon))
%         if (displayAll)
%             figure;
%         end
%         for ds=1:nDataSets
%             plot(data{ds}.pulsePhaseFlatness(mon,:),'LineWidth',lineWidthBig,'Color',datSetColors(ds,:));
%             hold all;
%         end
%         title([dataDescr 'Pulse Flatness, ' mixerLabels{mon}]);
%         legLabels = dataSetLabels;
%         for ds=1:nDataSets
%             legLabels{ds} = [dataSetLabels{ds} sprintf(' (std = %.2f^o)',data{ds}.meanPulsePhaseFlatness(mon))];
%         end
%         xlabel('Time [Pulse no.]');
%         ylabel('Phase [degrees]');
%         format_compare_legend;
%         format_plots;
%         if (savePlots)
%             saveName = [saveDir '/pulsePhaseFlatness_' mixerLabels{mon}];
%             print([saveName '.png'],'-dpng');
%             savefig([saveName '.fig']);
%         end
%         hold off;
%     end
% end
% 
% % mean phase vs. data set
% meanPulsePhase = NaN(nDataSets,nMons);
% for ds=1:nDataSets
%     meanPulsePhase(ds,:) = nanmean(data{ds}.meanPulsePhase,2);
% end
% if (displayAll)
%     a = figure;
% end
% for mon=1:nMons
%     if (plotMon(mon))
%         plot(dataSetValues, meanPulsePhase(:,mon),markerType,'MarkerSize',markerSizeBig);
%         hold all;
%     end
% end
% xlabel(dataSetValueLabel);
% ylabel('Phase [degrees]')
% title([dataDescr 'Mean Phase']);
% legend(mixerLabels(boolean(plotMon)));
% format_plots;
% if (savePlots)
%     saveName = [saveDir '/meanMeanPhase'];
%     print([saveName '.png'],'-dpng');
%     savefig([saveName '.fig']);
% end
% hold off;
% 
% 
% % jitter mean phase vs. data set
% stdMeanPulsePhase = NaN(nDataSets,nMons);
% stdMeanPulsePhase_err = NaN(nDataSets,nMons);
% for ds=1:nDataSets
%     stdMeanPulsePhase(ds,:) = data{ds}.stdMeanPulsePhase;
%     stdMeanPulsePhase_err(ds,:) = data{ds}.stdMeanPulsePhase_err;
% end
% 
% if (displayAll)
%     a = figure;
% end
% for mon=1:nMons
%     if (plotMon(mon))
%         plot(dataSetValues, stdMeanPulsePhase(:,mon),markerType,'MarkerSize',markerSizeBig);
%         hold all;
%     end
% end
% xlabel(dataSetValueLabel);
% ylabel('Phase Jitter [degrees]')
% title([dataDescr 'Jitter Mean Phase']);
% legend(mixerLabels(boolean(plotMon)));
% format_plots;
% if (savePlots)
%     saveName = [saveDir '/stdMeanPhase'];
%     print([saveName '.png'],'-dpng');
%     savefig([saveName '.fig']);
% end
% hold off;
% 
% % jitter ratio plots
% % if (displayAll)
% %     a = figure;
% % end
% % for ds=1:nDataSets
% %     plot(data{ds}.timeAxis,data{ds}.ratio12StdPhaseAlongPulse,'LineWidth',lineWidthBig);
% %     hold all;
% % end
% % xlabel('Time [ns]');
% % ylabel('Jitter Ratio')
% % title([dataDescr 'Jitter Along Pulse Ratio: ' mixerLabels{1} '/' mixerLabels{2}]);
% % legend(dataSetLabels);
% % format_plots;
% % if (savePlots)
% %     saveName = [saveDir '/stdMeanPhase'];
% %     print([saveName '.png'],'-dpng');
% %     savefig([saveName '.fig']);
% % end
% % 
% % hold off;
% 
% % ratio13StdPhaseAlongPulse
% % ratio23StdPhaseAlongPulse
% % 
% % ratio12StdMeanPulsePhase
% % ratio13StdMeanPulsePhase
% % ratio23StdMeanPulsePhase
% 
% 
% %% frascati phase - correlation mean phase
% 
% % extract correlation for each data set
% corrMeanMix1Mix2 = NaN(1,nDataSets);
% corrMeanMix1Mix3 = NaN(1,nDataSets);
% corrMeanMix2Mix3 = NaN(1,nDataSets);
% corrMeanMix1Mix2_err = NaN(1,nDataSets);
% corrMeanMix1Mix3_err = NaN(1,nDataSets);
% corrMeanMix2Mix3_err = NaN(1,nDataSets);
% fitMeanMix1Mix2 = NaN(nDataSets,2);
% fitMeanMix1Mix3 = NaN(nDataSets,2);
% fitMeanMix2Mix3 = NaN(nDataSets,2);
% for ds=1:nDataSets
%     corrMeanMix1Mix2(ds) = data{ds}.corrMeanMix1Mix2;
%     corrMeanMix1Mix3(ds) = data{ds}.corrMeanMix1Mix3;
%     corrMeanMix2Mix3(ds) = data{ds}.corrMeanMix2Mix3;
%     corrMeanMix1Mix2_err(ds) = data{ds}.corrMeanMix1Mix2_err;
%     corrMeanMix1Mix3_err(ds) = data{ds}.corrMeanMix1Mix3_err;
%     corrMeanMix2Mix3_err(ds) = data{ds}.corrMeanMix2Mix3_err;
%     fitMeanMix1Mix2(ds,:) = data{ds}.fitMeanMix1Mix2;
%     fitMeanMix1Mix3(ds,:) = data{ds}.fitMeanMix1Mix3;
%     fitMeanMix2Mix3(ds,:) = data{ds}.fitMeanMix2Mix3;
% 
% end
% 
% if (plotMon(1) && plotMon(2))
%     if (displayAll)
%         a = figure;
%     end
%     plot(dataSetValues, corrMeanMix1Mix2,markerType,'MarkerSize',markerSizeBig);
%     xlabel(dataSetValueLabel);
%     ylabel('Correlation');
%     title([dataDescr 'Correlation Mean Phase: ' mixerLabels{1} ' vs. ' mixerLabels{2}]);
%     format_plots;
%     if (savePlots)
%         saveName = [saveDir '/meanCorr_' mixerLabels{1} '_' mixerLabels{2}];
%         print([saveName '.png'],'-dpng');
%         savefig([saveName '.fig']);
%     end
% end
% 
% if (plotMon(1) && plotMon(3))
%     if (displayAll)
%         a = figure;
%     end
%     plot(dataSetValues, corrMeanMix1Mix3,markerType,'MarkerSize',markerSizeBig);
%     xlabel(dataSetValueLabel);
%     ylabel('Correlation');
%     title([dataDescr 'Correlation Mean Phase: ' mixerLabels{1} ' vs. ' mixerLabels{3}]);
%     format_plots;
%     if (savePlots)
%         saveName = [saveDir '/meanCorr_' mixerLabels{1} '_' mixerLabels{3}];
%         print([saveName '.png'],'-dpng');
%         savefig([saveName '.fig']);
%     end
% end
% 
% if (plotMon(2) && plotMon(3))
%     if (displayAll)
%         a = figure;
%     end
%     plot(dataSetValues, corrMeanMix2Mix3,markerType,'MarkerSize',markerSizeBig);
%     xlabel(dataSetValueLabel);
%     ylabel('Correlation');
%     title([dataDescr 'Correlation Mean Phase: ' mixerLabels{2} ' vs. ' mixerLabels{3}]);
%     format_plots;
%     if (savePlots)
%         saveName = [saveDir '/meanCorr_' mixerLabels{2} '_' mixerLabels{3}];
%         print([saveName '.png'],'-dpng');
%         savefig([saveName '.fig']);
%     end
% end
% 
% %% Mean sample correlation
% 
% % extract mean sample correlation for each data set
% meanCorrSamplesMix1Mix2 = NaN(1,nDataSets);
% meanCorrSamplesMix1Mix3 = NaN(1,nDataSets);
% meanCorrSamplesMix2Mix3 = NaN(1,nDataSets);
% for ds=1:nDataSets
%     meanCorrSamplesMix1Mix2(ds) = data{ds}.meanCorrSamplesMix1Mix2;
%     meanCorrSamplesMix1Mix3(ds) = data{ds}.meanCorrSamplesMix1Mix3;
%     meanCorrSamplesMix2Mix3(ds) = data{ds}.meanCorrSamplesMix2Mix3;
% end
% 
% % Mix1 Mix2
% if (plotMon(1) && plotMon(2))
%     if (displayAll)
%         a = figure;
%     end
%     plot(dataSetValues, meanCorrSamplesMix1Mix2,markerType,'MarkerSize',markerSizeBig);
%     xlabel(dataSetValueLabel);
%     ylabel('Correlation');
%     title([dataDescr 'Mean Sample Correlation: ' mixerLabels{1} ' vs. ' mixerLabels{2}]);
%     format_plots;
%     if (savePlots)
%         saveName = [saveDir '/meanSampleCorr_' mixerLabels{1} '_' mixerLabels{2}];
%         print([saveName '.png'],'-dpng');
%         savefig([saveName '.fig']);
%     end
% end
% 
% % Mix1 Mix3
% if (plotMon(1) && plotMon(3))
%     if (displayAll)
%         a = figure;
%     end
%     plot(dataSetValues, meanCorrSamplesMix1Mix3,markerType,'MarkerSize',markerSizeBig);
%     xlabel(dataSetValueLabel);
%     ylabel('Correlation');
%     title([dataDescr 'Mean Sample Correlation: ' mixerLabels{1} ' vs. ' mixerLabels{3}]);
%     format_plots;
%     if (savePlots)
%         saveName = [saveDir '/meanSampleCorr_' mixerLabels{1} '_' mixerLabels{3}];
%         print([saveName '.png'],'-dpng');
%         savefig([saveName '.fig']);
%     end
% end
% 
% % Mix2 Mix3
% if (plotMon(2) && plotMon(3))
%     if (displayAll)
%         a = figure;
%     end
%     plot(dataSetValues, meanCorrSamplesMix2Mix3,markerType,'MarkerSize',markerSizeBig);
%     xlabel(dataSetValueLabel);
%     ylabel('Correlation');
%     title([dataDescr 'Mean Sample Correlation: ' mixerLabels{2} ' vs. ' mixerLabels{3}]);
%     format_plots;
%     if (savePlots)
%         saveName = [saveDir '/meanSampleCorr_' mixerLabels{2} '_' mixerLabels{3}];
%         print([saveName '.png'],'-dpng');
%         savefig([saveName '.fig']);
%     end
% end
% 
% %% Correlation vs. sample no.
% 
% % Mix1 Mix2
% if (plotMon(1) && plotMon(2))
%     if (displayAll)
%         a = figure;
%     end
%     for ds=1:nDataSets
%         plot(data{ds}.timeAxis,data{ds}.corrSamplesMix1Mix2,'LineWidth',lineWidthBig,'Color',datSetColors(ds,:));
%         hold all;
%     end
%     title(sprintf('%sCorrelation vs. Sample No, %s vs %s',dataDescr,mixerLabels{1},mixerLabels{2}));
%     xlabel('Time along pulse [ns]')
%     ylabel('Correlation')
%     xlim([data{1}.timeAxis(plotStartSamp) data{1}.timeAxis(plotEndSamp)])
% 
%     legLabels = dataSetLabels;
%     for ds=1:nDataSets
%         legLabels{ds} = [dataSetLabels{ds} sprintf(' (%.2f)',data{ds}.meanCorrSamplesMix1Mix2)];
%     end
%     format_compare_legend;
% 
%     yLim = get(gca,'YLim');
%     if (yLim(1) < -1)
%         yLim(1) = -1;
%     end
%     if (yLim(2) > 1)
%         yLim(2) = 1;
%     end
%     ylim(yLim);
% 
%     plot([data{1}.timeAxis(min(data{1}.sampleRange)) data{1}.timeAxis(min(data{1}.sampleRange))],yLim,'k','LineWidth',lineWidthSmall);
%     plot([data{1}.timeAxis(max(data{1}.sampleRange)) data{1}.timeAxis(max(data{1}.sampleRange))],yLim,'k','LineWidth',lineWidthSmall);
% 
%     format_plots;
% 
%     if (savePlots)
%         saveName = [saveDir '/corrVsSample_' mixerLabels{1} '_' mixerLabels{2}];
%         print([saveName '.png'],'-dpng');
%         savefig([saveName '.fig']);
%     end
%     hold off;
% end
% 
% % Mix1 Mix3
% if (plotMon(1) && plotMon(3))
%     if (displayAll)
%         a = figure;
%     end
%     for ds=1:nDataSets
%         plot(data{ds}.timeAxis,data{ds}.corrSamplesMix1Mix3,'LineWidth',lineWidthBig,'Color',datSetColors(ds,:));
%         hold all;
%     end
% 
%     title(sprintf('%sCorrelation vs. Sample No, %s vs %s',dataDescr,mixerLabels{1},mixerLabels{3}));
%     xlabel('Time along pulse [ns]')
%     ylabel('Correlation')
%     xlim([data{1}.timeAxis(plotStartSamp) data{1}.timeAxis(plotEndSamp)])
% 
%     legLabels = dataSetLabels;
%     for ds=1:nDataSets
%         legLabels{ds} = [dataSetLabels{ds} sprintf(' (%.2f)',data{ds}.meanCorrSamplesMix1Mix3)];
%     end
%     format_compare_legend;
% 
%     yLim = get(gca,'YLim');
%     if (yLim(1) < -1)
%         yLim(1) = -1;
%     end
%     if (yLim(2) > 1)
%         yLim(2) = 1;
%     end
%     ylim(yLim);
% 
%     plot([data{1}.timeAxis(min(data{1}.sampleRange)) data{1}.timeAxis(min(data{1}.sampleRange))],yLim,'k','LineWidth',lineWidthSmall);
%     plot([data{1}.timeAxis(max(data{1}.sampleRange)) data{1}.timeAxis(max(data{1}.sampleRange))],yLim,'k','LineWidth',lineWidthSmall);
% 
%     format_plots;
% 
%     if (savePlots)
%         saveName = [saveDir '/corrVsSample_' mixerLabels{1} '_' mixerLabels{3}];
%         print([saveName '.png'],'-dpng');
%         savefig([saveName '.fig']);
%     end
%     hold off;
%     end
% 
% % Mix2 Mix3
% if (plotMon(2) && plotMon(3))
%     if (displayAll)
%         a = figure;
%     end
%     for ds=1:nDataSets
%         plot(data{ds}.timeAxis,data{ds}.corrSamplesMix2Mix3,'LineWidth',lineWidthBig,'Color',datSetColors(ds,:));
%         hold all;
%     end
% 
%     title(sprintf('%sCorrelation vs. Sample No, %s vs %s',dataDescr,mixerLabels{2},mixerLabels{3}));
%     xlabel('Time along pulse [ns]')
%     ylabel('Correlation')
%     xlim([data{1}.timeAxis(plotStartSamp) data{1}.timeAxis(plotEndSamp)])
% 
%     legLabels = dataSetLabels;
%     for ds=1:nDataSets
%         legLabels{ds} = [dataSetLabels{ds} sprintf(' (%.2f)',data{ds}.meanCorrSamplesMix2Mix3)];
%     end
%     format_compare_legend;
% 
%     yLim = get(gca,'YLim');
%     if (yLim(1) < -1)
%         yLim(1) = -1;
%     end
%     if (yLim(2) > 1)
%         yLim(2) = 1;
%     end
%     ylim(yLim);
% 
%     plot([data{1}.timeAxis(min(data{1}.sampleRange)) data{1}.timeAxis(min(data{1}.sampleRange))],yLim,'k','LineWidth',lineWidthSmall);
%     plot([data{1}.timeAxis(max(data{1}.sampleRange)) data{1}.timeAxis(max(data{1}.sampleRange))],yLim,'k','LineWidth',lineWidthSmall);
% 
%     format_plots;
% 
%     if (savePlots)
%         saveName = [saveDir '/corrVsSample_' mixerLabels{2} '_' mixerLabels{3}];
%         print([saveName '.png'],'-dpng');
%         savefig([saveName '.fig']);
%     end
%     hold off;
% end
% 
% %% frascati phase - mean correlation of pulse shape.
% 
% % extract mean sample correlation for each data set
% meanCorrShapeMix1Mix2 = NaN(1,nDataSets);
% meanCorrShapeMix1Mix3 = NaN(1,nDataSets);
% meanCorrShapeMix2Mix3 = NaN(1,nDataSets);
% for ds=1:nDataSets
%     meanCorrShapeMix1Mix2(ds) = data{ds}.meanCorrShapeMix1Mix2;
%     meanCorrShapeMix1Mix3(ds) = data{ds}.meanCorrShapeMix1Mix3;
%     meanCorrShapeMix2Mix3(ds) = data{ds}.meanCorrShapeMix2Mix3;
% end
% 
% % Mix1 Mix2
% if (plotMon(1) && plotMon(2))
%     if (displayAll)
%         a = figure;
%     end
%     plot(dataSetValues, meanCorrShapeMix1Mix2,markerType,'MarkerSize',markerSizeBig);
%     xlabel(dataSetValueLabel);
%     ylabel('Correlation');
%     title([dataDescr 'Mean Shape Correlation: ' mixerLabels{1} ' vs. ' mixerLabels{2}]);
%     format_plots;
%     if (savePlots)
%         saveName = [saveDir '/meanShapeCorr_' mixerLabels{1} '_' mixerLabels{2}];
%         print([saveName '.png'],'-dpng');
%         savefig([saveName '.fig']);
%     end
% end
% 
% % Mix1 Mix3
% if (plotMon(1) && plotMon(3))
%     if (displayAll)
%         a = figure;
%     end
%     plot(dataSetValues, meanCorrShapeMix1Mix3,markerType,'MarkerSize',markerSizeBig);
%     xlabel(dataSetValueLabel);
%     ylabel('Correlation');
%     title([dataDescr 'Mean Shape Correlation: ' mixerLabels{1} ' vs. ' mixerLabels{3}]);
%     format_plots;
%     if (savePlots)
%         saveName = [saveDir '/meanShapeCorr_' mixerLabels{1} '_' mixerLabels{3}];
%         print([saveName '.png'],'-dpng');
%         savefig([saveName '.fig']);
%     end
%     end
% 
% % Mix2 Mix3
% if (plotMon(2) && plotMon(3))
%     if (displayAll)
%         a = figure;
%     end
%     plot(dataSetValues, meanCorrShapeMix2Mix3,markerType,'MarkerSize',markerSizeBig);
%     xlabel(dataSetValueLabel);
%     ylabel('Correlation');
%     title([dataDescr 'Mean Shape Correlation: ' mixerLabels{2} ' vs. ' mixerLabels{3}]);
%     format_plots;
%     if (savePlots)
%         saveName = [saveDir '/meanShapeCorr_' mixerLabels{2} '_' mixerLabels{3}];
%         print([saveName '.png'],'-dpng');
%         savefig([saveName '.fig']);
%     end
% end
% 
% %% frascati phase - correlation pulse shape vs. time
% 
% % Mix1 Mix2
% if (plotMon(1) && plotMon(2))
%     if (displayAll)
%         a = figure;
%     end
%     for ds=1:nDataSets
%         plot(data{ds}.corrShapeMix1Mix2,'LineWidth',lineWidthBig,'Color',datSetColors(ds,:));
%         hold all;
%     end
% 
%     title(sprintf('%sCorrelation Pulse Shape vs. Time, %s vs. %s',dataDescr,mixerLabels{1},mixerLabels{2}));
%     xlabel('Time [Pulse no.]')
%     ylabel('Correlation')
% 
%     legLabels = dataSetLabels;
%     for ds=1:nDataSets
%         legLabels{ds} = [dataSetLabels{ds} sprintf(' (%.2f)',data{ds}.meanCorrShapeMix1Mix2)];
%     end
%     format_compare_legend;
% 
%     format_plots;
% 
%     if (savePlots)
%         saveName = [saveDir '/corrShape_' mixerLabels{1} '_' mixerLabels{2}];
%         print([saveName '.png'],'-dpng');
%         savefig([saveName '.fig']);
%     end
%     hold off;
% end
% 
% % Mix1 Mix3
% if (plotMon(1) && plotMon(3))
%     if (displayAll)
%         a = figure;
%     end
%     for ds=1:nDataSets
%         plot(data{ds}.corrShapeMix1Mix3,'LineWidth',lineWidthBig,'Color',datSetColors(ds,:));
%         hold all;
%     end
% 
%     title(sprintf('%sCorrelation Pulse Shape vs. Time, %s vs. %s',dataDescr,mixerLabels{1},mixerLabels{3}));
%     xlabel('Time [Pulse no.]')
%     ylabel('Correlation')
% 
%     legLabels = dataSetLabels;
%     for ds=1:nDataSets
%         legLabels{ds} = [dataSetLabels{ds} sprintf(' (%.2f)',data{ds}.meanCorrShapeMix1Mix3)];
%     end
%     format_compare_legend;
% 
%     format_plots;
% 
%     if (savePlots)
%         saveName = [saveDir '/corrShape_' mixerLabels{1} '_' mixerLabels{3}];
%         print([saveName '.png'],'-dpng');
%         savefig([saveName '.fig']);
%     end
%     hold off;
% end
% 
% % Mix2 Mix3
% if (plotMon(2) && plotMon(3))
%     if (displayAll)
%         a = figure;
%     end
%     for ds=1:nDataSets
%         plot(data{ds}.corrShapeMix2Mix3,'LineWidth',lineWidthBig,'Color',datSetColors(ds,:));
%         hold all;
%     end
% 
%     title(sprintf('%sCorrelation Pulse Shape vs. Time, %s vs. %s',dataDescr,mixerLabels{2},mixerLabels{3}));
%     xlabel('Time [Pulse no.]')
%     ylabel('Correlation')
% 
%     legLabels = dataSetLabels;
%     for ds=1:nDataSets
%         legLabels{ds} = [dataSetLabels{ds} sprintf(' (%.2f)',data{ds}.meanCorrShapeMix2Mix3)];
%     end
%     format_compare_legend;
% 
%     format_plots;
% 
%     if (savePlots)
%         saveName = [saveDir '/corrShape_' mixerLabels{2} '_' mixerLabels{3}];
%         print([saveName '.png'],'-dpng');
%         savefig([saveName '.fig']);
%     end
%     hold off;
% end
% 
% %% resolution
% 
% % Mon1 Mon2
% if (plotMon(1) && plotMon(2))
%     if (displayAll)
%         figure;
%     end
%     for ds=1:nDataSets
%         plot(data{ds}.timeAxis,data{ds}.resolution12,'LineWidth',2);
%         hold all;
%     end
%     xlim([data{1}.timeAxis(plotStartSamp) data{1}.timeAxis(plotEndSamp)])
%     yLim = get(gca,'YLim');
%     plot([data{1}.timeAxis(min(data{1}.sampleRange)) data{1}.timeAxis(min(data{1}.sampleRange))],yLim,'k');
%     plot([data{1}.timeAxis(max(data{1}.sampleRange)) data{1}.timeAxis(max(data{1}.sampleRange))],yLim,'k');
%     legLabels = dataSetLabels;
%     for ds=1:nDataSets
%         legLabels{ds} = [dataSetLabels{ds} sprintf(' (%.2f^o)',data{ds}.meanResolution12)];
%     end
%     title([dataDescr 'Resolution Mon1Mon2']);
%     xlabel('Time [ns]');
%     ylabel('Resolution [degrees]');
%     format_compare_legend;
%     format_plots;
%     if (savePlots)
%         saveName = [saveDir '/resolutionMon1Mon2'];
%         print([saveName '.png'],'-dpng');
%         savefig([saveName '.fig']);
%     end
%     hold off;
% end
% 
% % Mon1 Mon3
% if (plotMon(1) && plotMon(3))
%     if (displayAll)
%         figure;
%     end
%     for ds=1:nDataSets
%         plot(data{ds}.timeAxis,data{ds}.resolution13,'LineWidth',2);
%         hold all;
%     end
%     xlim([data{1}.timeAxis(plotStartSamp) data{1}.timeAxis(plotEndSamp)])
%     yLim = get(gca,'YLim');
%     plot([data{1}.timeAxis(min(data{1}.sampleRange)) data{1}.timeAxis(min(data{1}.sampleRange))],yLim,'k');
%     plot([data{1}.timeAxis(max(data{1}.sampleRange)) data{1}.timeAxis(max(data{1}.sampleRange))],yLim,'k');
%     legLabels = dataSetLabels;
%     for ds=1:nDataSets
%         legLabels{ds} = [dataSetLabels{ds} sprintf(' (%.2f^o)',data{ds}.meanResolution13)];
%     end
%     format_compare_legend;
%     format_plots;
%     title([dataDescr 'Resolution Mon1Mon3']);
%     xlabel('Time [ns]');
%     ylabel('Resolution [degrees]');
%     if (savePlots)
%         saveName = [saveDir '/resolutionMon1Mon3'];
%         print([saveName '.png'],'-dpng');
%         savefig([saveName '.fig']);
%     end
%     hold off;
% end
% 
% % Mon2 Mon3
% if (plotMon(2) && plotMon(3))
%     if (displayAll)
%         figure;
%     end
%     for ds=1:nDataSets
%         plot(data{ds}.timeAxis,data{ds}.resolution23,'LineWidth',2);
%         hold all;
%     end
%     xlim([data{1}.timeAxis(plotStartSamp) data{1}.timeAxis(plotEndSamp)])
%     yLim = get(gca,'YLim');
%     plot([data{1}.timeAxis(min(data{1}.sampleRange)) data{1}.timeAxis(min(data{1}.sampleRange))],yLim,'k');
%     plot([data{1}.timeAxis(max(data{1}.sampleRange)) data{1}.timeAxis(max(data{1}.sampleRange))],yLim,'k');
%     legLabels = dataSetLabels;
%     for ds=1:nDataSets
%         legLabels{ds} = [dataSetLabels{ds} sprintf(' (%.2f^o)',data{ds}.meanResolution23)];
%     end
%     format_compare_legend;
%     format_plots;
%     title([dataDescr 'Resolution Mon2Mon3']);
%     xlabel('Time [ns]');
%     ylabel('Resolution [degrees]');
%     if (savePlots)
%         saveName = [saveDir '/resolutionMon2Mon3'];
%         print([saveName '.png'],'-dpng');
%         savefig([saveName '.fig']);
%     end
%     hold off;
% end
% 
% 
% % extract mean resolution for each data set
% resolutionMix1Mix2 = NaN(1,nDataSets);
% resolutionMix1Mix3 = NaN(1,nDataSets);
% resolutionMix2Mix3 = NaN(1,nDataSets);
% for ds=1:nDataSets
%     resolutionMix1Mix2(ds) = data{ds}.meanResolution12;
%     resolutionMix1Mix3(ds) = data{ds}.meanResolution13;
%     resolutionMix2Mix3(ds) = data{ds}.meanResolution23;
% end
% 
% % Mix1 Mix2
% if (displayAll)
%     a = figure;
% end
% plot(dataSetValues, resolutionMix1Mix2,markerType,'MarkerSize',markerSizeBig);
% xlabel(dataSetValueLabel);
% ylabel('Resolution [degrees]');
% title([dataDescr 'Mean Resolution Mon1Mon2']);
% format_plots;
% if (savePlots)
%     saveName = [saveDir '/meanResolution12VsDataSet'];
%     print([saveName '.png'],'-dpng');
%     savefig([saveName '.fig']);
% end
% 
