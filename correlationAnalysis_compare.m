r56_03 = load('/home/jack/PhaseFeedforward/Analysis/CorrelationAnalysis/R56/20141215_1411_allForStability_R56_+0.3/correlationAnalysis.mat');%load('/home/jack/PhaseFeedforward/Analysis/CorrelationAnalysis/R56/20141215_1411_allForStability_R56_+0.3/correlationAnalysis.mat'); %load('/home/jack/PhaseFeedforward/Analysis/CorrelationAnalysis/R56/20141205_1330_R56_+0.3/correlationAnalysis.mat');
r56_05 = load('/home/jack/PhaseFeedforward/Analysis/CorrelationAnalysis/R56/20141205_1125_R56_+0.5/correlationAnalysis.mat');
r56_neg02_posGain = load('/home/jack/PhaseFeedforward/Analysis/CorrelationAnalysis/R56/20141205_1125_R56_-0.2_Gain26/correlationAnalysis.mat');
r56_neg02_negGain = load('/home/jack/PhaseFeedforward/Analysis/CorrelationAnalysis/R56/20141205_1125_R56_-0.2_Gain-26/correlationAnalysis.mat');
close all


%% mon3 vs. ct.285h
tmpPhase = [r56_neg02_negGain.meanMon3Phase' r56_neg02_posGain.meanMon3Phase'];
phase_neg02 = tmpPhase-nanmean(tmpPhase);
bpmCT285Index = find(ismember(r56_neg02_posGain.bpmNames,'CT.SVBPM0285'));
tmpPos = [r56_neg02_negGain.meanBPMH(bpmCT285Index,:) r56_neg02_posGain.meanBPMH(bpmCT285Index,:)];
pos_neg02 = tmpPos-nanmean(tmpPos);

phase03 = r56_03.meanMon3Phase-nanmean(r56_03.meanMon3Phase);
bpmCT285Index = find(ismember(r56_03.bpmNames,'CT.SVBPM0285'));
pos03 = r56_03.meanBPMH(bpmCT285Index,:)-nanmean(r56_03.meanBPMH(bpmCT285Index,:));

phase05 = r56_05.meanMon3Phase-nanmean(r56_05.meanMon3Phase);
bpmCT285Index = find(ismember(r56_05.bpmNames,'CT.SVBPM0285'));
pos05 = r56_05.meanBPMH(bpmCT285Index,:)-nanmean(r56_05.meanBPMH(bpmCT285Index,:));

figure;
plot(phase_neg02,pos_neg02,'o','Color',[0,0.8,0]);
hold all;
plot(phase03,pos03,'bx');
plot(phase05,pos05,'r^');
legend('R56 = -0.2','R56 = +0.3','R56 = +0.5');
xlabel('Mon3 Phase [degrees]');
ylabel('CT.285 H [mm]');
axXLim = get(gca,'XLim');
axYLim = get(gca,'YLim');

goodPulse = ~(isnan(phase_neg02) | isnan(pos_neg02));
fit_neg02 = polyfit(phase_neg02(goodPulse),pos_neg02(goodPulse),1);

goodPulse = ~(isnan(phase03)' | isnan(pos03));
fit03 = polyfit(phase03(goodPulse)',pos03(goodPulse),1);

goodPulse = ~(isnan(phase05)' | isnan(pos05));
fit05 = polyfit(phase05(goodPulse)',pos05(goodPulse),1);

xDat = axXLim(1):axXLim(2);
yDat = fit_neg02(1).*xDat + fit_neg02(2);
plot(xDat,yDat,'Color',[0,0.8,0],'LineWidth',2);

xDat = axXLim(1):axXLim(2);
yDat = fit03(1).*xDat + fit03(2);
plot(xDat,yDat,'b','LineWidth',2);

xDat = axXLim(1):axXLim(2);
yDat = fit05(1).*xDat + fit05(2);
plot(xDat,yDat,'r','LineWidth',2);

ylim(axYLim);

%% mon3 correlation with bpmH

corr_neg02 = nanmean([r56_neg02_negGain.corrMon3_BPMH' r56_neg02_posGain.corrMon3_BPMH'],2);
corr03 = r56_03.corrMon3_BPMH;
corr05 = r56_05.corrMon3_BPMH;

xAxCorr03 = 1:length(r56_03.bpmNames);

xAxCorr_neg02 = NaN(1,length(r56_neg02_negGain.bpmNames));
for i=1:length(r56_neg02_negGain.bpmNames)
    xAxCorr_neg02(i) = find(ismember(r56_03.bpmNames,r56_neg02_negGain.bpmNames{i}));
end

xAxCorr05 = NaN(1,length(r56_05.bpmNames));
for i=1:length(r56_05.bpmNames)
    xAxCorr05(i) = find(ismember(r56_03.bpmNames,r56_05.bpmNames{i}));
end

figure;
plot(xAxCorr_neg02(1:18),corr_neg02(1:18),'Color',[0 0.8 0],'Linewidth',2);
hold all;
plot(xAxCorr03,corr03,'b','Linewidth',2);
plot(xAxCorr05(1:18),corr05(1:18),'r','Linewidth',2);
legend('R56 = -0.2','R56 = +0.3','R56 = +0.5')
plot(xAxCorr_neg02(19:end),corr_neg02(19:end),'Color',[0 0.8 0],'Linewidth',2);
plot(xAxCorr05(19:end),corr05(19:end),'r','Linewidth',2);
xlabel('BPM Index')
ylabel('Correlation')
title('Correlation with BPM H for Different R56 in TL1')
grid