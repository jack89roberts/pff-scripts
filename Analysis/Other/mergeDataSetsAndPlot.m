% % run plotData_compareDatasets on these files first
dataSetNames = {...
...% ff off
'20151120_1526_Gain-950_R56_0.1_Interleaaved_Odd',...
'20151120_1538_Gain-800_R56_0.1_Interleaaved_Odd',...
'20151120_1551_Gain-750_R56_0.1_Interleaved_Odd',...
'20151120_1600_Gain950_R56_0.1_Interleaaved_Odd',...
'20151120_1604_Gain-775_R56_0.1_Interleaved_Even',...
'20151120_1611_Gain-775_R56_0.1_Interleaved_Even',...
'20151120_1617_Gain-825_R56_0.1_Interleaved_Odd',...
'20151120_1622_Gain-825_R56_0.1_Interleaved_Even',...
'20151120_1639_Gain-790_R56_0.1_Interleaved_Odd',...
'20151120_1644_Gain-790_R56_0.1_Interleaved_Even',...
'20151120_1650_Gain-810_R56_0.1_Interleaved_Odd',...
'20151120_1654_Gain-810_R56_0.1_Interleaved_Odd',...
'20151120_1700_Gain-850_R56_0.1_Interleaved_Even',...
'20151120_1721_Gain-850_R56_0.1_Interleaved_Odd',...
'20151120_1800_Gain-850_R56_0.1_Interleaved_Odd'...
...% ff on
% '20151120_1526_Gain-950_R56_0.1_Interleaaved_Even',...
% '20151120_1538_Gain-800_R56_0.1_Interleaaved_Even',...
% '20151120_1551_Gain-750_R56_0.1_Interleaved_Even',...
% '20151120_1600_Gain950_R56_0.1_Interleaaved_Even',...
% '20151120_1604_Gain-775_R56_0.1_Interleaved_Odd',...
% '20151120_1611_Gain-775_R56_0.1_Interleaved_Odd',...
% '20151120_1617_Gain-825_R56_0.1_Interleaved_Even',...
% '20151120_1622_Gain-825_R56_0.1_Interleaved_Odd',...
% '20151120_1639_Gain-790_R56_0.1_Interleaved_Even',...
% '20151120_1644_Gain-790_R56_0.1_Interleaved_Odd',...
% '20151120_1650_Gain-810_R56_0.1_Interleaved_Even',...
% '20151120_1654_Gain-810_R56_0.1_Interleaved_Even',...
% '20151120_1700_Gain-850_R56_0.1_Interleaved_Odd',...
% '20151120_1721_Gain-850_R56_0.1_Interleaved_Even',...
% '20151120_1800_Gain-850_R56_0.1_Interleaved_Even',...
};
%% run when ff on files loaded

ct608HistoryFFOn = [];
mon2HistoryFFOn = []; 
mon3HistoryFFOn = [];
nPulsesFFOn = NaN(1,length(data));
stdMon2FFOnDatSet = NaN(1,length(data));
stdMon3FFOnDatSet = NaN(1,length(data));
stdMon2FFOnDatSet_err = NaN(1,length(data));
stdMon3FFOnDatSet_err = NaN(1,length(data));

for i=1:nDataSets
    ct608HistoryFFOn = [ct608HistoryFFOn; data{i}.meanBPMH{data{i}.refBPM2Index}];
    mon2HistoryFFOn = [mon2HistoryFFOn, data{i}.meanPulsePhase(2,:)]; 
    mon3HistoryFFOn = [mon3HistoryFFOn, data{i}.meanPulsePhase(3,:)];
    nPulsesFFOn(i) = data{i}.nPulses;
    stdMon2FFOnDatSet(i) = data{i}.stdMeanPulsePhase(2);
    stdMon3FFOnDatSet(i) = data{i}.stdMeanPulsePhase(3);
    stdMon2FFOnDatSet_err(i) = data{i}.stdMeanPulsePhase_err(2);
    stdMon3FFOnDatSet_err(i) = data{i}.stdMeanPulsePhase_err(3);

end

[~,stdMon2FFOn,~,stdMon2FFOn_err] = nanMeanStdErr(mon2HistoryFFOn);
[~,stdMon3FFOn,~,stdMon3FFOn_err] = nanMeanStdErr(mon3HistoryFFOn);
[corrFFOn,corrFFOn_err] = nancorrcoef(mon2HistoryFFOn, mon3HistoryFFOn);
fitFFOn = nanpolyfit(mon2HistoryFFOn, mon3HistoryFFOn,1);

ct608HistoryFFOn = ct608HistoryFFOn - nanmean(ct608HistoryFFOn);

figure;
plot(mon2HistoryFFOn);
hold all;
plot(mon3HistoryFFOn);
plot(ct608HistoryFFOn);

%% run when ff off files loaded
ct608HistoryFFOff = [];
mon2HistoryFFOff = []; 
mon3HistoryFFOff = [];
nPulsesFFOff = NaN(1,length(data));
stdMon2FFOffDatSet = NaN(1,length(data));
stdMon3FFOffDatSet = NaN(1,length(data));
stdMon2FFOffDatSet_err = NaN(1,length(data));
stdMon3FFOffDatSet_err = NaN(1,length(data));
corrFFOffDatSet = NaN(1,length(data));
corrFFOffDatSet_err = NaN(1,length(data));

for i=1:nDataSets
    ct608HistoryFFOff = [ct608HistoryFFOff; data{i}.meanBPMH{data{i}.refBPM2Index}];
    mon2HistoryFFOff = [mon2HistoryFFOff, data{i}.meanPulsePhase(2,:)]; 
    mon3HistoryFFOff = [mon3HistoryFFOff, data{i}.meanPulsePhase(3,:)];
    nPulsesFFOff(i) = data{i}.nPulses;
    stdMon2FFOffDatSet(i) = data{i}.stdMeanPulsePhase(2);
    stdMon3FFOffDatSet(i) = data{i}.stdMeanPulsePhase(3);
    stdMon2FFOffDatSet_err(i) = data{i}.stdMeanPulsePhase_err(2);
    stdMon3FFOffDatSet_err(i) = data{i}.stdMeanPulsePhase_err(3);
    corrFFOffDatSet(i) = data{i}.corrMeanMix2Mix3;
    corrFFOffDatSet_err(i) = data{i}.corrMeanMix2Mix3_err;
end

ct608HistoryFFOff = ct608HistoryFFOff - nanmean(ct608HistoryFFOff);

normMon2HistoryFFOff = mon2HistoryFFOff-nanmean(mon2HistoryFFOff);
normMon2HistoryFFOff = normMon2HistoryFFOff./max(abs(normMon2HistoryFFOff));
normMon3HistoryFFOff = mon3HistoryFFOff-nanmean(mon3HistoryFFOff);
normMon3HistoryFFOff = normMon3HistoryFFOff./max(abs(normMon3HistoryFFOff));

[~,stdMon2FFOff,~,stdMon2FFOff_err] = nanMeanStdErr(mon2HistoryFFOff);
[~,stdMon3FFOff,~,stdMon3FFOff_err] = nanMeanStdErr(mon3HistoryFFOff);
[corrFFOff,corrFFOff_err] = nancorrcoef(mon2HistoryFFOff, mon3HistoryFFOff);
fitFFOff = nanpolyfit(mon2HistoryFFOff, mon3HistoryFFOff,1);

figure;
plot(mon2HistoryFFOff);
hold all;
plot(mon3HistoryFFOff);
plot(ct608HistoryFFOff);

%% mean phase history
figure;
plot(mon2HistoryFFOff,'Color',[0 0.8 0]);
hold all;
plot(mon3HistoryFFOff,'b');
plot(mon3HistoryFFOn,'r');

datSetBoundaries = cumsum(nPulsesFFOn);
for i=1:length(data);
    plot([datSetBoundaries(i) datSetBoundaries(i)],[-30 30],'k','LineWidth',1);
end
% ylim([-25 0]);
ylim([-5 5]);

datSetMid = datSetBoundaries-([datSetBoundaries(1) diff(datSetBoundaries)]/2);
for i=1:length(data)
    tmpTime = strsplit(dataSetNames{i},'_');
    tmpTime = tmpTime{2};
    tmpTime = [tmpTime(1:2) ':' tmpTime(3:4)];
    %text(datSetMid(i),-2,tmpTime,'HorizontalAlignment','center','FontWeight','bold')
    text(datSetMid(i),-4,tmpTime,'HorizontalAlignment','center','FontWeight','bold')

end
xlabel('Pulse No.')
ylabel('Phase [degrees]')
legend('Upstream','Downstream FF Off','Downstream FF On')
title('20th November 2015: Mean Phase History')
format_plots;
grid off;
ax = gca;
set(ax,'YGrid','on');


%% ff sim (run when ff off files loaded)

simRange = 5.5;

% opt gain, unlimited correction range
mon3FFSimUnlim = [];
optGain = NaN(1,nDataSets);
simUnlimJitDatset = NaN(1,nDataSets);
simUnlimJitDatset_err = NaN(1,nDataSets);
for i=1:nDataSets
    simFFResults = getSimulatedFF(data{i});
%     simFFResults = getSimulatedFF(data{i},[],6,-2);
    optGain(i) = simFFResults.gain;
    mon3FFSimUnlim = [mon3FFSimUnlim; simFFResults.meanSimFFPhase];
    simUnlimJitDatset(i) = simFFResults.stdSimFFPhase;
    simUnlimJitDatset_err(i) = simFFResults.stdSimFFPhase_err;

end
[~,stdSimUnlim,~,stdSimUnlim_err] = nanMeanStdErr(mon3FFSimUnlim);
[corrSimUnlim,corrSimUnlim_err] = nancorrcoef(mon2HistoryFFOff, mon3FFSimUnlim);
fitSimUnlim = nanpolyfit(mon2HistoryFFOff, mon3FFSimUnlim,1);


% opt gain, limited 6 degree correction range
mon3FFSim = [];
optGain = NaN(1,nDataSets);
simOptJitDatset = NaN(1,nDataSets);
simOptJitDatset_err = NaN(1,nDataSets);
for i=1:nDataSets
    simFFResults = getSimulatedFF(data{i},[],simRange);
%     simFFResults = getSimulatedFF(data{i},[],6,-2);
    optGain(i) = simFFResults.gain;
    mon3FFSim = [mon3FFSim; simFFResults.meanSimFFPhase];
    simOptJitDatset(i) = simFFResults.stdSimFFPhase;
    simOptJitDatset_err(i) = simFFResults.stdSimFFPhase_err;

end
% for i=1:nDataSets
%     tmpGain = (data{i}.corrMeanMix2Mix3).*(data{i}.stdMeanPulsePhase(3)./data{i}.stdMeanPulsePhase(2));
%     tmpFFSim = data{i}.meanPulsePhase(3,:)-tmpGain.*data{i}.meanPulsePhase(2,:);
%     mon3FFSim = [mon3FFSim, tmpFFSim];
% end

[~,stdSimOpt,~,stdSimOpt_err] = nanMeanStdErr(mon3FFSim);
[corrSimOpt,corrSimOpt_err] = nancorrcoef(mon2HistoryFFOff, mon3FFSim);
fitSimOpt = nanpolyfit(mon2HistoryFFOff, mon3FFSim,1);

% real gain
mon3FFSimRealGain = [];
simRealJitDatset = NaN(1,nDataSets);
simRealJitDatset_err = NaN(1,nDataSets);
realGainFONT = [950 800 750 950 775 775 825 825 790 790 810 810 850 850 850];
realGain = realGainFONT./624;%710.64;
for i=1:nDataSets
    simFFResults = getSimulatedFF(data{i},realGain(i),simRange);
    mon3FFSimRealGain = [mon3FFSimRealGain; simFFResults.meanSimFFPhase];
    simRealJitDatset(i) = simFFResults.stdSimFFPhase;
    simRealJitDatset_err(i) = simFFResults.stdSimFFPhase_err;
end
% for i=1:nDataSets
%     tmpFFSim = data{i}.meanPulsePhase(3,:)-realGain(i).*data{i}.meanPulsePhase(2,:);
%     mon3FFSimRealGain = [mon3FFSimRealGain, tmpFFSim];
% end

[~,stdSimReal,~,stdSimReal_err] = nanMeanStdErr(mon3FFSimRealGain);
[corrSimReal,corrSimReal_err] = nancorrcoef(mon2HistoryFFOff, mon3FFSimRealGain);
fitSimReal = nanpolyfit(mon2HistoryFFOff, mon3FFSimRealGain,1);


% opt gain with offset
offsets = -[0.5658 -0.2350 -1.8099 -0.6079 -1.4615 -0.7200 -1.0374 -1.9211 -2.9281 -3.2253 -2.4033 -1.0277 -1.0294 -1.9405 -1.5774];

mon3FFSimOffset = [];
simOffsetJitDatset = NaN(1,nDataSets);
simOffsetJitDatset_err = NaN(1,nDataSets);
for i=1:nDataSets
    simFFResults = getSimulatedFF(data{i},[],simRange,offsets(i)); % 
    mon3FFSimOffset = [mon3FFSimOffset; simFFResults.meanSimFFPhase];
    simOffsetJitDatset(i) = simFFResults.stdSimFFPhase;
    simOffsetJitDatset_err(i) = simFFResults.stdSimFFPhase_err;
end
% for i=1:nDataSets
%     tmpFFSim = data{i}.meanPulsePhase(3,:)-realGain(i).*data{i}.meanPulsePhase(2,:);
%     mon3FFSimRealGain = [mon3FFSimRealGain, tmpFFSim];
% end

[~,stdSimRealOffset ,~,stdSimRealOffset_err] = nanMeanStdErr(mon3FFSimOffset);
[corrSimRealOffset,corrSimRealOffset_err] = nancorrcoef(mon2HistoryFFOff, mon3FFSimOffset);
fitSimRealOffset = nanpolyfit(mon2HistoryFFOff, mon3FFSimOffset,1);



% real gain with offset calculated from FFOn FFOff difference.
% NB: Offsets must be calculated with both ff on and ff off data loaded,
% and using data without mean subtraction.
if (length(mon3HistoryFFOff)==length(mon3HistoryFFOn))
    diffMon3FFOnOff = mon3HistoryFFOff-mon3HistoryFFOn;
elseif (length(mon3HistoryFFOff)>length(mon3HistoryFFOn))
    diffMon3FFOnOff = mon3HistoryFFOff(1:end-1)-mon3HistoryFFOn;
elseif (length(mon3HistoryFFOff)==length(mon3HistoryFFOn))
    diffMon3FFOnOff = mon3HistoryFFOff-mon3HistoryFFOn(1:end-1);
end

% offsets = NaN(1,nDataSets);
% offsets(1) = nanmean(diffMon3FFOnOff(1:datSetBoundaries(1)));
% for i=2:nDataSets
%     offsets(i) = nanmean(diffMon3FFOnOff((datSetBoundaries(i-1)+1):datSetBoundaries(i)));
% end

offsets = -[0.5658 -0.2350 -1.8099 -0.6079 -1.4615 -0.7200 -1.0374 -1.9211 -2.9281 -3.2253 -2.4033 -1.0277 -1.0294 -1.9405 -1.5774];

mon3FFSimALL = [];
simALLJitDatset = NaN(1,nDataSets);
simALLJitDatset_err = NaN(1,nDataSets);
for i=1:nDataSets
    simFFResults = getSimulatedFF(data{i},realGain(i),simRange,offsets(i)); % 
    mon3FFSimALL = [mon3FFSimALL; simFFResults.meanSimFFPhase];
    simALLJitDatset(i) = simFFResults.stdSimFFPhase;
    simALLJitDatset_err(i) = simFFResults.stdSimFFPhase_err;
end
% for i=1:nDataSets
%     tmpFFSim = data{i}.meanPulsePhase(3,:)-realGain(i).*data{i}.meanPulsePhase(2,:);
%     mon3FFSimRealGain = [mon3FFSimRealGain, tmpFFSim];
% end

[~,stdSimALL ,~,stdSimALL_err] = nanMeanStdErr(mon3FFSimALL);
[corrSimALL,corrSimALL_err] = nancorrcoef(mon2HistoryFFOff, mon3FFSimALL);
fitSimRealOffset = nanpolyfit(mon2HistoryFFOff, mon3FFSimALL,1);

%%
figure;
plot(simUnlimJitDatset);
hold all;
plot(simOptJitDatset);
plot(simRealJitDatset);
plot(simOffsetJitDatset);
plot(simALLJitDatset);
% plot(simOffsetJitDatset-simUnlimJitDatset,'k--')
% plot(stdMon3FFOnDatSet,'k','LineWidth',3);
legend('Unlimited','6deg Range','Real Gain','Real Offset','All effects')
ylabel('Downstream Phase Jitter [degrees]')
xlabel('Dataset')
title('Simulated Phase Jitter')
set(gca,'xtick',1:2:15)
set(gca,'xticklabel',{'15:26','15:51','16:04','16:17','16:39','16:50','17:00','18:00'})
rotateXLabels(gca,30)

figure;
plot(simOptJitDatset-simUnlimJitDatset);
hold all;
plot(simRealJitDatset-simUnlimJitDatset);
plot(simOffsetJitDatset-simUnlimJitDatset);
plot(simALLJitDatset-simUnlimJitDatset);
% plot(stdMon3FFOnDatSet,'k','LineWidth',3);
legend('sim opt','sim real','sim offset','sim all')

figure;
% a=shadedErrorBar(1:15,simOptJitDatset,simOptJitDatset_err,'b',0.5);
% hold all;
% b=shadedErrorBar(1:15,simRealJitDatset,simRealJitDatset_err,'r',0.5);
e=shadedErrorBar(1:15,stdMon3FFOffDatSet,stdMon3FFOffDatSet_err,'b',0.5);
hold all;
d=shadedErrorBar(1:15,stdMon3FFOnDatSet,stdMon3FFOnDatSet_err,'r',0.5);
c=shadedErrorBar(1:15,simALLJitDatset,simALLJitDatset_err,'g',0.5);
f=shadedErrorBar(1:15,stdMon2FFOnDatSet,stdMon2FFOnDatSet_err,'m',0.5);
% legend([a.mainLine b.mainLine c.mainLine d.mainLine],'SimOpt','SimReal','SimOffset','Real')
legend([e.mainLine d.mainLine c.mainLine f.mainLine],'FF Off','FF On','Sim All Effects','Upstream')
format_plots;
ylabel('Phase Jitter [degrees]')
xlabel('Dataset')
title('Achieved Phase Jitter')
set(gca,'xtick',1:2:15)
set(gca,'xticklabel',{'15:26','15:51','16:04','16:17','16:39','16:50','17:00','18:00'})
rotateXLabels(gca,30)
ylim([0.2 2.2])
set(gca,'ytick',0.2:0.4:2.2)

figure;
plot(stdMon2FFOffDatSet);
hold all;
plot(stdMon2FFOnDatSet);

%% scatter plots

figure;
plot(mon2HistoryFFOff,mon3HistoryFFOff,'o','MarkerFaceColor','b','MarkerEdgeColor',[0.5 1 0]);
xlim([-3 3])
ylim([-5 5])
xlabel('Upstream Phase [degrees]');
ylabel('Downstream Phase [degrees]');
title({'FF OFF', sprintf('(std = %.2f%c%.2f^o, corr = %.2f%c%.2f)',stdMon3FFOff,char(177),stdMon3FFOff_err,corrFFOff,char(177),corrFFOff_err)});
format_plots;
set(gca,'YTick',-5:5)

figure;
plot(mon2HistoryFFOn,mon3HistoryFFOn,'o','MarkerFaceColor','r','MarkerEdgeColor',[1 1 0]);
xlim([-3 3])
ylim([-5 5])
xlabel('Upstream Phase [degrees]');
ylabel('Downstream Phase [degrees]');
title({'FF ON',sprintf('(std = %.2f%c%.2f^o, corr = %.2f%c%.2f)',stdMon3FFOn,char(177),stdMon3FFOn_err,corrFFOn,char(177),corrFFOn_err)});
format_plots;
set(gca,'YTick',-5:5)

figure;
plot(mon2HistoryFFOff,mon3FFSimUnlim,'o','MarkerFaceColor','g','MarkerEdgeColor',[0 0 0]);
xlim([-3 3])
ylim([-5 5])
xlabel('Upstream Phase [degrees]');
ylabel('Downstream Phase [degrees]');
title({'FF Sim Unlimited',sprintf('(std = %.2f%c%.2f^o, corr = %.2f%c%.2f)',stdSimUnlim,char(177),stdSimUnlim_err,corrSimUnlim,char(177),corrSimUnlim_err)});
format_plots;
set(gca,'YTick',-5:5)

figure;
plot(mon2HistoryFFOff,mon3FFSim,'o','MarkerFaceColor','g','MarkerEdgeColor',[0 0 0]);
xlim([-3 3])
ylim([-5 5])
xlabel('Upstream Phase [degrees]');
ylabel('Downstream Phase [degrees]');
title({'FF Sim 6deg Range',sprintf('(std = %.2f%c%.2f^o, corr = %.2f%c%.2f)',stdSimOpt,char(177),stdSimOpt_err,corrSimOpt,char(177),corrSimOpt_err)});
format_plots;
set(gca,'YTick',-5:5)

figure;
plot(mon2HistoryFFOff,mon3FFSimRealGain,'o','MarkerFaceColor','k','MarkerEdgeColor',[0.7 0.7 0.7]);
xlim([-3 3])
ylim([-5 5])
xlabel('Upstream Phase [degrees]');
ylabel('Downstream Phase [degrees]');
title({'FF Sim Real Gain',sprintf('(std = %.2f%c%.2f^o, corr = %.2f%c%.2f)',stdSimReal,char(177),stdSimReal_err,corrSimReal,char(177),corrSimReal_err)});
format_plots;
set(gca,'YTick',-5:5)

figure;
plot(mon2HistoryFFOff,mon3FFSimOffset,'o','MarkerFaceColor','c','MarkerEdgeColor','m');
xlim([-3 3])
ylim([-5 5])
xlabel('Upstream Phase [degrees]');
ylabel('Downstream Phase [degrees]');
title({'FF Sim Used Gain With Offset',sprintf('(std = %.2f%c%.2f^o, corr = %.2f%c%.2f)',stdSimRealOffset,char(177),stdSimRealOffset_err,corrSimRealOffset,char(177),corrSimRealOffset_err)});
format_plots;
set(gca,'YTick',-5:5)

figure;
plot(mon2HistoryFFOff,mon3FFSimALL,'o','MarkerFaceColor','k','MarkerEdgeColor',[0.7 0.7 0.7]);
xlim([-3 3])
ylim([-5 5])
xlabel('Upstream Phase [degrees]');
ylabel('Downstream Phase [degrees]');
title({'FF Sim All Effects',sprintf('(std = %.2f%c%.2f^o, corr = %.2f%c%.2f)',stdSimALL,char(177),stdSimALL_err,corrSimALL,char(177),corrSimALL_err)});
format_plots;
set(gca,'YTick',-5:5)

%%
bpmNearMon2 = [];
for i=1:nDataSets
    bpmNearMon2 = [bpmNearMon2; data{i}.meanBPMH{9}];
end
figure;
plot(bpmNearMon2);

bpmNearMon3 = [];
for i=1:nDataSets
    bpmNearMon3 = [bpmNearMon3; data{i}.meanBPMH{43}];
end
figure;
plot(bpmNearMon3);



%% bpm loop and plot

colours = varycolor(60);
figure;
for bpm=1:60
    tmpBPMH = [];
    for i=1:nDataSets
        tmpBPMH = [tmpBPMH; data{i}.meanBPMH{bpm}];
    end
    tmpBPMH = tmpBPMH-nanmean(tmpBPMH);
    tmpBPMH = tmpBPMH./max(abs(tmpBPMH));
    %tmpBPMH = tmpBPMH+bpm/2;
    %plot(tmpBPMH,'Color',colours(bpm,:))
    subplot(2,2,1:2)
    plot(tmpBPMH,'b')
    hold all;
    plot(normMon2HistoryFFOff,'r');
    plot(normMon3HistoryFFOff,'Color',[0 0.8 0],'LineWidth',1);
    hold off;
    title(data{1}.bpmNames{bpm});
    legend('BPM','Upstream Phase','Downstream Phase');
    xlabel('Pulse No.')
    ylabel('Normalised Output [a.u.]')
    subplot(2,2,3)
    plot(tmpBPMH,'b')
    hold all;
    plot(normMon2HistoryFFOff,'r');
    plot(normMon3HistoryFFOff,'Color',[0 0.8 0],'LineWidth',1);
    xlabel('Pulse No.')
    ylabel('Normalised Output [a.u.]')
    hold off;
    title('Region of Interest');
    xlim([350 550])
    subplot(2,2,4)
    plot(tmpBPMH,'b')
    hold all;
    plot(-normMon2HistoryFFOff,'r');
    plot(-normMon3HistoryFFOff,'Color',[0 0.8 0],'LineWidth',1);
    xlabel('Pulse No.')
    ylabel('Normalised Output [a.u.]')
    hold off;
    title('Region of Interest (phase sign flipped)');
    xlim([350 550])
    input(sprintf('BPM %d...', bpm))
%     savePlot('/home/jack/PhaseFeedforward/Analysis/201511/Plots/20151120_AllAfternoon/BPMsNormalised',sprintf('BPM%d',bpm));
    %hold all;
end

% colours = varycolor(1542);
% c=1;
% figure;
% for i=1:nDataSets
%     for p=1:data{i}.nPulses
%         tmpOrbit = NaN(1,60);
%         for bpm=1:60
%             try
%                 tmpOrbit(bpm) = data{i}.meanBPMH{bpm}(p);
%             catch
%             end
%         end
%         plot(tmpOrbit,'Color',colours(c,:));
%         hold all;
%         c=c+1;
% 
%     end
% end


%% bpm SUM loop and plot

colours = varycolor(60);
figure;
for bpm=1:60
    tmpBPMS = [];
    for i=1:nDataSets
        tmpBPMS = [tmpBPMS; data{i}.meanBPMS{bpm}];
    end
    tmpBPMS = tmpBPMS-nanmean(tmpBPMS);
    tmpBPMS = tmpBPMS./max(abs(tmpBPMS));
    subplot(2,2,1:2)
    plot(tmpBPMS,'b')
    hold all;
    plot(normMon2HistoryFFOff,'r');
    plot(normMon3HistoryFFOff,'Color',[0 0.8 0],'LineWidth',1);
    hold off;
    title(['BEAM CURRENT ' data{1}.bpmNames{bpm}]);
    legend('BPM','Upstream Phase','Downstream Phase');
    xlabel('Pulse No.')
    ylabel('Normalised Output [a.u.]')
    subplot(2,2,3)
    plot(tmpBPMS,'b')
    hold all;
    plot(normMon2HistoryFFOff,'r');
    plot(normMon3HistoryFFOff,'Color',[0 0.8 0],'LineWidth',1);
    xlabel('Pulse No.')
    ylabel('Normalised Output [a.u.]')
    hold off;
    title('Region of Interest');
    xlim([350 550])
    subplot(2,2,4)
    plot(tmpBPMS,'b')
    hold all;
    plot(-normMon2HistoryFFOff,'r');
    plot(-normMon3HistoryFFOff,'Color',[0 0.8 0],'LineWidth',1);
    xlabel('Pulse No.')
    ylabel('Normalised Output [a.u.]')
    hold off;
    title('Region of Interest (phase sign flipped)');
    xlim([350 550])
    input(sprintf('BPM %d...', bpm))
    savePlot('/home/jack/PhaseFeedforward/Analysis/201511/Plots/20151120_AllAfternoon/BPMCurrentNormalised',sprintf('BPM%d',bpm));
    %hold all;
end

%%
timeStrings = cell(1,nDataSets);
for i=1:length(data)
    tmpTime = strsplit(dataSetNames{i},'_');
    tmpTime = tmpTime{2};
    tmpTime = [tmpTime(1:2) ':' tmpTime(3:4)];
    timeStrings{i} = tmpTime;
end

figure;
plot(offsets,'b','LineWidth',2);
set(gca,'XTick',1:2:nDataSets);
set(gca,'XTickLabel',timeStrings(1:2:end));
rotateXLabels(gca,30)
xlabel('Dataset')
ylabel('Offset [degrees]')
title('Upstream Phase Offset wrt PFF Window')
format_plots; 

figure;
plot(realGain,'b','LineWidth',2);
hold all;
plot(optGain,'r','LineWidth',2);
set(gca,'XTick',1:2:nDataSets);
set(gca,'XTickLabel',timeStrings(1:2:end));
rotateXLabels(gca,30)
xlabel('Dataset')
ylabel('Gain')
legend('Used Gain','Optimal Gain')
format_plots; 

%%
correlsFFOn = NaN(1,nDataSets);
correlsFFOn_err = NaN(1,nDataSets);
for i=1:nDataSets
    correlsFFOn(i) = data{i}.corrMeanMix2Mix3;
    correlsFFOn_err(i) = data{i}.corrMeanMix2Mix3_err;
end

correlsFFOff = NaN(1,nDataSets);
correlsFFOff_err = NaN(1,nDataSets);
for i=1:nDataSets
    correlsFFOff(i) = data{i}.corrMeanMix2Mix3;
    correlsFFOff_err(i) = data{i}.corrMeanMix2Mix3_err;
end

figure;
a=shadedErrorBar(1:15,correlsFFOff,correlsFFOff_err,'b');
hold all;
b=shadedErrorBar(1:15,correlsFFOn,correlsFFOn_err,'r');
set(gca,'XTick',1:2:nDataSets);
set(gca,'XTickLabel',timeStrings(1:2:end));
rotateXLabels(gca,30)
xlabel('Dataset')
ylabel('Correlation')
legend([a.mainLine,b.mainLine],{'FF Off','FF On'})
format_plots; 
title('Correlation')

%% correction saturation (no mean sub and fixed correction window)
subMon2HistoryFFOff = mon2HistoryFFOff-nanmean(mon2HistoryFFOff);
bestGain = nanmean(optGain);
corrRange = 6;
offset = 0;

% figure;
% plot(subMon2HistoryFFOff);
% hold all;
% plot([0 1600],[6 6],'k','LineWidth',2);
% plot([0 1600],[-6 -6],'k','LineWidth',2);

mon2PhasesFFOff = [];
mon3PhasesFFOff = [];
for i=1:nDataSets
    mon2PhasesFFOff = [mon2PhasesFFOff; squeeze(data{i}.phases(2,:,:))];
    mon3PhasesFFOff = [mon3PhasesFFOff; squeeze(data{i}.phases(3,:,:))];
end
mon2PhasesFFOff = mon2PhasesFFOff-nanmean(mon2HistoryFFOff);
mon3PhasesFFOff = mon3PhasesFFOff-nanmean(mon3HistoryFFOff);

figure;
plot(mon2PhasesFFOff'.*bestGain)
hold all
plot([200 800],[-6 -6],'k','LineWidth',2);
plot([200 800],[6 6],'k','LineWidth',2);
plot([602 602], [-60 80],'k','LineWidth',2)
plot([682 682], [-60 80],'k','LineWidth',2)
ylabel('Upstream Phase*Gain [degrees]')
xlabel('Sample No.')

figure;
plot(mon3PhasesFFOff')
hold all
plot([602 602], [-60 80],'k','LineWidth',2)
plot([682 682], [-60 80],'k','LineWidth',2)
ylabel('Downstream Phase [degrees]')
xlabel('Sample No.')

outRange = ((mon2PhasesFFOff+offset).*bestGain>corrRange)|((mon2PhasesFFOff+offset)*bestGain<-corrRange);
fracInRange = 1 - (sum(outRange)./length(mon2HistoryFFOff));

% figure;
plot(fracInRange);