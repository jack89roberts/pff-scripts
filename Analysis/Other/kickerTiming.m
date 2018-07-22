close all;
%%
% K1 and K2
% kickOnDataSets = {...
% '20151113_1629_Gain_2000_Delay0_Interleaved_Even',...
% '20151113_1632_Gain_2000_Delay3_Interleaved_Even',...
% '20151113_1634_Gain_2000_Delay6_Interleaved_Even',...
% '20151113_1637_Gain_2000_Delay9_Interleaved_Even',...
% '20151113_1639_Gain_2000_Delay12_Interleaved_Even',...
% '20151113_1642_Gain_2000_Delay15_Interleaved_Even',...
% '20151113_1644_Gain_2000_Delay18_Interleaved_Even'...
% };
% 
% kickOffDataSets = {...
% '20151113_1629_Gain_2000_Delay0_Interleaved_Odd',...
% '20151113_1632_Gain_2000_Delay3_Interleaved_Odd',...
% '20151113_1634_Gain_2000_Delay6_Interleaved_Odd',...
% '20151113_1637_Gain_2000_Delay9_Interleaved_Odd',...
% '20151113_1639_Gain_2000_Delay12_Interleaved_Odd',...
% '20151113_1642_Gain_2000_Delay15_Interleaved_Odd',...
% '20151113_1644_Gain_2000_Delay18_Interleaved_Odd'...
% };

kickOnDataSets = {...
    '20151012_2046_Gain_K1_4000_K2_4000_Delay_K1_0_K2_1_Odd',...
    '20151012_2050_Gain_K1_4000_K2_4000_Delay_K1_10_K2_11_Even',...
    '20151012_2053_Gain_K1_4000_K2_4000_Delay_K1_20_K2_21_Odd',...
    '20151012_2056_Gain_K1_4000_K2_4000_Delay_K1_30_K2_31_Even'...
};
kickOffDataSets = {...
    '20151012_2046_Gain_K1_4000_K2_4000_Delay_K1_0_K2_1_Even',...
    '20151012_2050_Gain_K1_4000_K2_4000_Delay_K1_10_K2_11_Odd',...
    '20151012_2053_Gain_K1_4000_K2_4000_Delay_K1_20_K2_21_Even',...
    '20151012_2056_Gain_K1_4000_K2_4000_Delay_K1_30_K2_31_Odd'...
};

% K1 only
% kickOnDataSets = {...
%     '20151012_2103_Gain_K1_4000_K2_0_Delay_K1_0_K2_0_Even',...
%     '20151012_2107_Gain_K1_4000_K2_0_Delay_K1_10_K2_0_Odd',...
%     '20151012_2108_Gain_K1_4000_K2_0_Delay_K1_20_K2_0_Even',...
%     '20151012_2110_Gain_K1_4000_K2_0_Delay_K1_30_K2_0_Odd'...
% };
% kickOffDataSets = {...
%     '20151012_2103_Gain_K1_4000_K2_0_Delay_K1_0_K2_0_Odd',...
%     '20151012_2107_Gain_K1_4000_K2_0_Delay_K1_10_K2_0_Even',...
%     '20151012_2108_Gain_K1_4000_K2_0_Delay_K1_20_K2_0_Odd',...
%     '20151012_2110_Gain_K1_4000_K2_0_Delay_K1_30_K2_0_Even'...
% };
% 
% K2 only
% kickOnDataSets = {...
%     '20151012_2118_Gain_K1_0_K2_4000_Delay_K1_0_K2_1_Even',...
%     '20151012_2117_Gain_K1_0_K2_4000_Delay_K1_0_K2_11_Even',...
%     '20151012_2113_Gain_K1_0_K2_4000_Delay_K1_0_K2_21_Odd',...
%     '20151012_2112_Gain_K1_0_K2_4000_Delay_K1_0_K2_31_Even'...
% };
% kickOffDataSets = {...
%     '20151012_2118_Gain_K1_0_K2_4000_Delay_K1_0_K2_1_Odd',...
%     '20151012_2117_Gain_K1_0_K2_4000_Delay_K1_0_K2_11_Odd',...
%     '20151012_2113_Gain_K1_0_K2_4000_Delay_K1_0_K2_21_Even',...
%     '20151012_2112_Gain_K1_0_K2_4000_Delay_K1_0_K2_31_Odd'...
% };

% dataBaseDir = '/home/jack/PhaseFeedforward/Analysis/201510';
savePlotDir = '/home/jack/PhaseFeedforward/Analysis/201510/KickerTimingScan/K1K2';
% dataBaseDir = '/user/ctf3op/PhaseFeedforward/Processed';
dataBaseDir = '/home/jack/PhaseFeedforward/Analysis/201510';

savePlots = 0;

kickDelaysFONT = [0 10 20 30];
% kickDelaysFONT = [0 3 6 9 12 15 18];
% kickDelaysFONT = [1 11 21 31];

fontSampInt = 1000/357;

bpmIndex = 43;
sampleRangeToUse = 620:700;
plotXRange = [600 720];
plotYRange = [-5 2];
%%
addpath('../');
addpath('../../');

nDataSets = length(kickOnDataSets);
kickDelaysTime = kickDelaysFONT.*fontSampInt;

for datInd=1:nDataSets
    % load kick on data
    load([dataBaseDir '/' kickOnDataSets{datInd} '.mat']);
    
    % initialise arrays
    if (datInd==1)
        nBPMSamples = length(meanBPMSAlongPulse{bpmIndex});
        nPhaseSamples = length(meanPhaseAlongPulse(2,:));
        
        kickOnBPM = NaN(nDataSets,nBPMSamples);
        kickOffBPM = NaN(nDataSets,nBPMSamples);
        kickOnDownstreamPhase = NaN(nDataSets,nPhaseSamples);
        
        kickOnBPMS = NaN(nDataSets,nBPMSamples);
        kickOffBPMS = NaN(nDataSets,nBPMSamples);
        ffKickBPM = NaN(nDataSets,nBPMSamples);
        downstreamPhase = NaN(nDataSets,nPhaseSamples);
        upstreamPhase = NaN(nDataSets,nPhaseSamples);
        downstreamDiode = NaN(nDataSets,nPhaseSamples);
        upstreamDiode = NaN(nDataSets,nPhaseSamples);
        
        normKickOnBPMS = NaN(nDataSets,nBPMSamples);
        normKickOffBPMS = NaN(nDataSets,nBPMSamples);
        normFFKickBPM = NaN(nDataSets,nBPMSamples);
        normDownstreamPhase = NaN(nDataSets,nPhaseSamples);
        normUpstreamPhase = NaN(nDataSets,nPhaseSamples);
        normDownstreamDiode = NaN(nDataSets,nPhaseSamples);
        normUpstreamDiode = NaN(nDataSets,nPhaseSamples);      
    end
    
    kickOnBPM(datInd,:) = meanBPMHAlongPulse{bpmIndex};
    kickOnBPMS(datInd,:) = meanBPMSAlongPulse{bpmIndex};
    kickOnDownstreamPhase(datInd,:) = meanPhaseAlongPulse(3,:);
    
    % load kick off data
    load([dataBaseDir '/' kickOffDataSets{datInd} '.mat']);
    kickOffBPM(datInd,:) = meanBPMHAlongPulse{bpmIndex};
    kickOffBPMS(datInd,:) = meanBPMSAlongPulse{bpmIndex};
    downstreamPhase(datInd,:) = meanPhaseAlongPulse(3,:);
    upstreamPhase(datInd,:) = meanPhaseAlongPulse(2,:);
    meanDiodes = squeeze(nanmean(diodes,2));
    downstreamDiode(datInd,:) = meanDiodes(3,:);
    upstreamDiode(datInd,:) = meanDiodes(2,:);
    
    ffKickBPM(datInd,:) = kickOnBPM(datInd,:)-kickOffBPM(datInd,:);

    % make all phase sags same sign
    fitUpstreamPhase = polyfit(sampleRangeToUse,upstreamPhase(datInd,sampleRangeToUse),2);
    fitDownstreamPhase = polyfit(sampleRangeToUse,downstreamPhase(datInd,sampleRangeToUse),2);
    fitFFKickBPM = polyfit(sampleRangeToUse,ffKickBPM(datInd,sampleRangeToUse),2);    
    if (sign(fitUpstreamPhase(1))==1); upstreamPhase(datInd,:)= -upstreamPhase(datInd,:); end;
    if (sign(fitDownstreamPhase(1))==1); downstreamPhase(datInd,:)= -downstreamPhase(datInd,:); end;
    if (sign(fitFFKickBPM(1))==1); ffKickBPM(datInd,:) = -ffKickBPM(datInd,:); end;
    
    % normalise (to max, not absolute max, in sample range)
    normUpstreamPhase(datInd,:) = upstreamPhase(datInd,:)-nanmean(upstreamPhase(datInd,sampleRangeToUse));
    normUpstreamPhase(datInd,:) = normUpstreamPhase(datInd,:)./(max(normUpstreamPhase(datInd,sampleRangeToUse)));
    
    normDownstreamPhase(datInd,:) = downstreamPhase(datInd,:)-nanmean(downstreamPhase(datInd,sampleRangeToUse));
    normDownstreamPhase(datInd,:) = normDownstreamPhase(datInd,:)./(max(normDownstreamPhase(datInd,sampleRangeToUse)));
    
    normFFKickBPM(datInd,:) = ffKickBPM(datInd,:)-nanmean(ffKickBPM(datInd,sampleRangeToUse));
    normFFKickBPM(datInd,:) = normFFKickBPM(datInd,:)./(max(normFFKickBPM(datInd,sampleRangeToUse)));

    
    normKickOffBPMS(datInd,:) = kickOffBPMS(datInd,:);%-nanmean(kickOffBPMS(sampleRangeToUse));
    normKickOffBPMS(datInd,:) = normKickOffBPMS(datInd,:)./max(abs(normKickOffBPMS(datInd,sampleRangeToUse)));
    
    normKickOnBPMS(datInd,:) = kickOnBPMS(datInd,:);%-nanmean(kickOnBPMS(sampleRangeToUse));
    normKickOnBPMS(datInd,:) = normKickOnBPMS(datInd,:)./max(abs(normKickOnBPMS(datInd,sampleRangeToUse)));

    normUpstreamDiode(datInd,:) = upstreamDiode(datInd,:);%-nanmean(upstreamDiode(sampleRangeToUse));
    normUpstreamDiode(datInd,:) = normUpstreamDiode(datInd,:)./max(abs(normUpstreamDiode(datInd,sampleRangeToUse)));

    normDownstreamDiode(datInd,:) = downstreamDiode(datInd,:);%-nanmean(downstreamDiode(sampleRangeToUse));
    normDownstreamDiode(datInd,:) = normDownstreamDiode(datInd,:)./max(abs(normDownstreamDiode(datInd,sampleRangeToUse)));
    
end

diffDownstreamPhase = kickOnDownstreamPhase-downstreamPhase;
figure;
for datInd=1:nDataSets
    plot(diffDownstreamPhase);
    hold all;
end
% return;

% Find delay based on location of max value
% This works only for data where a wiggle along the pulse gives the max
% output in the sample range!!!
[~,upstreamPeakLocs] = max(normUpstreamPhase(:,sampleRangeToUse),[],2);
[~,downstreamPeakLocs] = max(normDownstreamPhase(:,sampleRangeToUse),[],2);
[~,ffKickBPMPeakLocs] = max(normFFKickBPM(:,sampleRangeToUse),[],2);

ffKickBPMDelays = (ffKickBPMPeakLocs-upstreamPeakLocs).*sampInterval; % only valid if phase mon sampling freq and bpm sampling freq is the same
[fitFFKickBPMDelays,fitrsq,fitconf] = nanpolyfit(kickDelaysTime,ffKickBPMDelays',1);
fitFFKickBPMDelays_err = (fitFFKickBPMDelays-fitconf(1,:))/2;
fitOptimalDelay = -fitFFKickBPMDelays(2)./fitFFKickBPMDelays(1);
fitOptimalDelayFONT = round(fitOptimalDelay./fontSampInt);

figure;
plot(kickDelaysTime,ffKickBPMDelays,'o','MarkerFaceColor','b');
hold all;
plot(kickDelaysTime,fitFFKickBPMDelays(1).*kickDelaysTime + fitFFKickBPMDelays(2),'k','LineWidth',2);
xlabel('Applied Output Delay [ns]');
ylabel('Kick - Beam Offset [ns]');
legend('Data',sprintf('Fit: Gradient: %.1f%c%.1f, Offset: %.0f%c%.0f ns',fitFFKickBPMDelays(1),char(177),fitFFKickBPMDelays_err(1),fitFFKickBPMDelays(2),char(177),fitFFKickBPMDelays_err(2)),'Location','North');
title('Fitted Offest Between Kick and Beam')
ylim([-40 65])
xlim([0 85])
format_plots;
if savePlots; savePlot(savePlotDir,'fitResults'); end;

for datInd = 1:nDataSets
    figure;
    plot(normFFKickBPM(datInd,:),'b','LineWidth',2);
    hold all;
    plot(normUpstreamPhase(datInd,:),'g','LineWidth',2);
    plot(normDownstreamPhase(datInd,:),'r','LineWidth',2);
    legend('BPM', 'Upstream Phase','Downstream Phase','Location','South');
    xlim(plotXRange);
    xlabel('Sample No.');
    ylabel('Output [a.u.]');
    ylim(plotYRange);
    plot([sampleRangeToUse(upstreamPeakLocs(datInd)) sampleRangeToUse(upstreamPeakLocs(datInd))],plotYRange,'k','LineWidth',2);
    plot([sampleRangeToUse(ffKickBPMPeakLocs(datInd)) sampleRangeToUse(ffKickBPMPeakLocs(datInd))],plotYRange,'k','LineWidth',2);
    title(sprintf('Applied Delay: %.0f ns, Calculated Offset: %.0f ns',kickDelaysTime(datInd),ffKickBPMDelays(datInd)));
    ylim([-3 1.5])
    format_plots;
    if savePlots; savePlot(savePlotDir,sprintf('delay%.0fPhase',kickDelaysFONT(datInd))); end;


    figure;
    plot(normKickOnBPMS(datInd,:),'b','LineWidth',2);
    hold all;
%     plot(normKickOffBPMS(datInd,:));
    plot(normUpstreamDiode(datInd,:),'g','LineWidth',2);
    plot(normDownstreamDiode(datInd,:),'r','LineWidth',2);
%     legend('Kick on BPM Transmission','Kick off BPM transmission', 'Upstream Diode','Downstream Diode');
    legend('BPM', 'Upstream','Downstream','Location','North');
    xlim([530 730]);
    xlabel('Sample No.');
    ylabel('Output [a.u.]')
    yLim = get(gca,'YLim');
    plot([sampleRangeToUse(upstreamPeakLocs(datInd)) sampleRangeToUse(upstreamPeakLocs(datInd))],yLim,'k','LineWidth',2);
    plot([sampleRangeToUse(ffKickBPMPeakLocs(datInd)) sampleRangeToUse(ffKickBPMPeakLocs(datInd))],yLim,'k','LineWidth',2);
    title(sprintf('Transmission, Applied Delay: %.0f ns',kickDelaysTime(datInd)));
    format_plots;
    if savePlots; savePlot(savePlotDir,sprintf('delay%.0fDiode',kickDelaysFONT(datInd))); end;
end

nSamplesToDelay = -ffKickBPMDelays(1)./sampInterval;
% nSamplesToDelay = round(fitOptimalDelay./sampInterval);
delayedFFKickBPM = delaySignal(normFFKickBPM(1,:),nSamplesToDelay);
delayedDownstream = delaySignal(normDownstreamPhase(1,:),-1);
optGain = corrMeanMix2Mix3.*(stdMeanPulsePhase(3,:)./stdMeanPulsePhase(2,:));
delayedNoNormDownstream = delaySignal(downstreamPhase(1,:),-1);
simFFDownstream = delayedNoNormDownstream-optGain.*upstreamPhase(1,:);
normSimFF = simFFDownstream-nanmean(simFFDownstream(sampleRangeToUse));
normSimFF = normSimFF./(max(normSimFF(sampleRangeToUse)));


figure;
plot(delayedFFKickBPM,'b','LineWidth',2);
hold all;
plot(normUpstreamPhase(1,:),'g','LineWidth',2);
plot(delayedDownstream,'r','LineWidth',2);
% plot(delayedDownstream-delayedFFKickBPM);
% plot(normSimFF);
plot([sampleRangeToUse(upstreamPeakLocs(datInd)) sampleRangeToUse(upstreamPeakLocs(datInd))],get(gca,'YLim'),'k','LineWidth',2);
legend('BPM', 'Upstream Phase','Downstream Phase','Location','NorthWest');
xlim(plotXRange);
xlabel('Sample No.');
ylabel('Output [a.u.]');
ylim([-3 1.5])
title('Optimal Delay')
format_plots;
if savePlots; savePlot(savePlotDir,sprintf('optimalDelay_%.0f',fitOptimalDelay)); end;