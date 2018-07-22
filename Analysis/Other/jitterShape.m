load('/home/jack/PhaseFeedforward/Analysis/201507/20150715_1414_FF_Gain63_Gate185_350_Int_Odd.mat');
saveDir = '/home/jack/PhaseFeedforward/Analysis/201507/ShapeJitter/20150715_1414_FFOff';
longSampleRange = 545:705;
%%
close all;
addpath('../../');

diffPhases = NaN(size(phases));
for mon=1:nMons
    phases(mon,:,:) = phases(mon,:,:)-subtractPhase(mon);
    
    for p=1:nPulses
        diffPhases(mon,p,:) = squeeze(phases(mon,p,:))' - meanPhaseAlongPulse(mon,:);
    end   
end
meanDiffPhases = nanmean(diffPhases(:,:,sampleRange),3);
for mon=1:nMons
    for p=1:nPulses
        diffPhases(mon,p,:) = diffPhases(mon,p,:)-meanDiffPhases(mon,p);
    end
end

stdDiffPhases = nanstd(diffPhases(:,:,sampleRange),[],3);
meanStdDiffPhases = nanmean(stdDiffPhases,2);
longStdDiffPhases = nanstd(diffPhases(:,:,longSampleRange),[],3);
longMeanStdDiffPhases = nanmean(longStdDiffPhases,2);

meanLongPhase = nanmean(phases(:,:,longSampleRange),3);
stdMeanLongPhase = nanstd(meanLongPhase,[],2);

figure;
plot(squeeze(phases(2,:,:))');
hold all;
plot(meanPhaseAlongPulse(2,:),'k','LineWidth',2.5);
xlabel('Sample No.');
ylabel('Phase [degrees]')
title('All Upstream Phases (mean in black)');
xlim([pulseSampleRange{3}(1) pulseSampleRange{3}(end)]);
plot([sampleRange(1) sampleRange(1)],get(gca,'YLim'),'k');
plot([sampleRange(end) sampleRange(end)],get(gca,'YLim'),'k');
plot([longSampleRange(1) longSampleRange(1)],get(gca,'YLim'),'k--');
plot([longSampleRange(end) longSampleRange(end)],get(gca,'YLim'),'k--');
savePlot(saveDir,'allUpstream');

figure;
plot(squeeze(phases(3,:,:))');
hold all;
plot(meanPhaseAlongPulse(3,:),'k','LineWidth',2.5);
xlabel('Sample No.');
ylabel('Phase [degrees]')
title('All Downstream Phases (mean in black)');
xlim([pulseSampleRange{3}(1) pulseSampleRange{3}(end)]);
plot([sampleRange(1) sampleRange(1)],get(gca,'YLim'),'k');
plot([sampleRange(end) sampleRange(end)],get(gca,'YLim'),'k');
plot([longSampleRange(1) longSampleRange(1)],get(gca,'YLim'),'k--');
plot([longSampleRange(end) longSampleRange(end)],get(gca,'YLim'),'k--');
savePlot(saveDir,'allDownstream');

figure;
plot(squeeze(diffPhases(2,:,:))');
hold all;
xlabel('Sample No.');
ylabel('Diff Phase [degrees]')
title('Difference to Mean Phase Upstream');
ylim([-15 15]);
xlim([pulseSampleRange{3}(1) pulseSampleRange{3}(end)]);
plot([sampleRange(1) sampleRange(1)],get(gca,'YLim'),'k');
plot([sampleRange(end) sampleRange(end)],get(gca,'YLim'),'k');
plot([longSampleRange(1) longSampleRange(1)],get(gca,'YLim'),'k--');
plot([longSampleRange(end) longSampleRange(end)],get(gca,'YLim'),'k--');
savePlot(saveDir,'diffUpstream');

figure;
plot(squeeze(diffPhases(3,:,:))');
hold all;
xlabel('Sample No.');
ylabel('Diff Phase [degrees]')
title('Difference to Mean Phase Downstream');
ylim([-15 15]);
xlim([pulseSampleRange{3}(1) pulseSampleRange{3}(end)]);
plot([sampleRange(1) sampleRange(1)],get(gca,'YLim'),'k');
plot([sampleRange(end) sampleRange(end)],get(gca,'YLim'),'k');
plot([longSampleRange(1) longSampleRange(1)],get(gca,'YLim'),'k--');
plot([longSampleRange(end) longSampleRange(end)],get(gca,'YLim'),'k--');
savePlot(saveDir,'diffDownstream');


figure;
plot(longStdDiffPhases(2,:),'b');
hold all;
plot(longStdDiffPhases(3,:),'r');
plot([1 nPulses], [stdMeanLongPhase(2) stdMeanLongPhase(2)],'b--');
plot([1 nPulses], [stdMeanLongPhase(3) stdMeanLongPhase(3)],'r--');
xlabel('Pulse No.');
ylabel('Phase Jitter [degrees]');
legLabels = cell(1,4);
legLabels{1} = sprintf('Upstream Shape (%.2f^o)',longMeanStdDiffPhases(2));
legLabels{2} = sprintf('Downstream Shape (%.2f^o)',longMeanStdDiffPhases(3));
legLabels{3} = sprintf('Upstream Mean (%.2f^o)',stdMeanLongPhase(2));
legLabels{4} = sprintf('Downstream Mean (%.2f^o)',stdMeanLongPhase(3));
legend(legLabels,'Location','Best');
title({'Comparison Jitter of Pulse Shape and Jitter of Pulse Mean','Long Range'}) 
savePlot(saveDir,'jitterLong');

figure;
plot(stdDiffPhases(2,:),'b');
hold all;
plot(stdDiffPhases(3,:),'r');
plot([1 nPulses], [stdMeanPulsePhase(2) stdMeanPulsePhase(2)],'b--');
plot([1 nPulses], [stdMeanPulsePhase(3) stdMeanPulsePhase(3)],'r--');
xlabel('Pulse No.');
ylabel('Phase Jitter [degrees]');
legLabels = cell(1,4);
legLabels{1} = sprintf('Upstream Shape (%.2f^o)',meanStdDiffPhases(2));
legLabels{2} = sprintf('Downstream Shape (%.2f^o)',meanStdDiffPhases(3));
legLabels{3} = sprintf('Upstream Mean (%.2f^o)',stdMeanPulsePhase(2));
legLabels{4} = sprintf('Downstream Mean (%.2f^o)',stdMeanPulsePhase(3));
legend(legLabels,'Location','Best');
title({'Comparison Jitter of Pulse Shape and Jitter of Pulse Mean','Short Range'}) 
savePlot(saveDir,'jitterShort');
