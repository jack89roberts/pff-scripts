mon1On = squeeze(data{1}.phases(1,:,:));
mon3On = squeeze(data{1}.phases(3,:,:));

mon1Off = squeeze(data{2}.phases(1,:,:));
mon3Off = squeeze(data{2}.phases(3,:,:));

meanMon1On = nanmean(mon1On(:,463:843),2);
meanMon1Off = nanmean(mon1Off(:,463:843),2);

figure;
plot(meanMon1On,data{1}.meanPulsePhase(3,:),'bo');
hold all;
plot(meanMon1Off,data{2}.meanPulsePhase(3,:),'ro');

corrOn = nancorrcoef(meanMon1On,data{1}.meanPulsePhase(3,:))
corrOff = nancorrcoef(meanMon1Off,data{2}.meanPulsePhase(3,:))


meanMon1Off = meanMon1Off-nanmean(meanMon1Off);
meanMon3Off = data{2}.meanPulsePhase(3,:)-nanmean(data{2}.meanPulsePhase(3,:));
stdMon1Off= nanstd(meanMon1Off);
stdMon3Off = nanstd(meanMon3Off);

figure;
plot(meanMon1Off./stdMon1Off)
hold all
plot(meanMon3Off./stdMon3Off)
legend('Mon1','Mon3')
title('20161214.1410.x4')
xlabel('Pulse No.')
ylabel('Phase [degs]')


%% sub-pulse ranges

endUp = 848;
sub4 = (endUp-99):endUp;
sub3 = (endUp-199):(endUp-100);
sub2 = (endUp-299):(endUp-200);
sub1 = (endUp-399):(endUp-300);

plotRange = [395 900];

figure;
subplot(1,2,1)
plot(meanDiodeAlongPulse(1,:,:),'b','LineWidth',2)
hold all;
plot([endUp endUp],get(gca,'YLim'),'k--','LineWidth',2);
plot([endUp endUp]-99,get(gca,'YLim'),'k--','LineWidth',2);
plot([endUp endUp]-199,get(gca,'YLim'),'k--','LineWidth',2);
plot([endUp endUp]-299,get(gca,'YLim'),'k--','LineWidth',2);
plot([endUp endUp]-399,get(gca,'YLim'),'k--','LineWidth',2);
text(490,-1000,'1','FontSize',24)
text(590,-1000,'2','FontSize',24)
text(690,-1000,'3','FontSize',24)
text(790,-1000,'4','FontSize',24)
xlim(plotRange)
title('Upstream Diode')
xlabel('Sample No.')
ylabel('Output [counts]')

subplot(1,2,2)
plot(meanPhaseAlongPulse(1,:,:),'b','LineWidth',2)
hold all;
plot([endUp endUp],get(gca,'YLim'),'k--','LineWidth',2);
plot([endUp endUp]-99,get(gca,'YLim'),'k--','LineWidth',2);
plot([endUp endUp]-199,get(gca,'YLim'),'k--','LineWidth',2);
plot([endUp endUp]-299,get(gca,'YLim'),'k--','LineWidth',2);
plot([endUp endUp]-399,get(gca,'YLim'),'k--','LineWidth',2);
xlim(plotRange)
title('Upstream Phase')
xlabel('Sample No.')
ylabel('Phase [degrees]')
text(490,-30,'1','FontSize',24)
text(590,-30,'2','FontSize',24)
text(690,-30,'3','FontSize',24)
text(790,-30,'4','FontSize',24)

figure;
subplot(1,2,1)
plot(meanPhaseAlongPulse(1,sub1),'LineWidth',2);
hold all;
plot(meanPhaseAlongPulse(1,sub2),'LineWidth',2);
plot(meanPhaseAlongPulse(1,sub3),'LineWidth',2);
plot(meanPhaseAlongPulse(1,sub4),'LineWidth',2);
xlabel('Sample No.')
ylabel('Phase [degrees]')
xlim([1 100])
legend('Sub-pulse 1','Sub-pulse 2','Sub-pulse 3','Sub-pulse 4')
title('Upstream Sub-pulse Phases')

subplot(1,2,2)
sumSubs = meanPhaseAlongPulse(1,sub1);
sumSubs = sumSubs+meanPhaseAlongPulse(1,sub2);
sumSubs = sumSubs+meanPhaseAlongPulse(1,sub3);
sumSubs = sumSubs+meanPhaseAlongPulse(1,sub4);
sumSubs = sumSubs./4;
plot(sumSubs,'k','LineWidth',2);
title('Predicted Downstream Phase')
xlabel('Sample No.')
ylabel('Phase [degrees]')
xlim([1 100])