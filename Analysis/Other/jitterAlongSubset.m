pulseRange = 294:343;
plotRange = [135 1080];

pffOff = squeeze(data{1}.phases(3,pulseRange,:));
pffOn = squeeze(data{2}.phases(3,pulseRange,:));
tAx = data{2}.timeAxis;

[pffOffMeanAlong,pffOffStdAlong]=nanMeanStdErr(pffOff);
[pffOnMeanAlong,pffOnStdAlong]=nanMeanStdErr(pffOn);

figure;
plot(tAx,pffOffMeanAlong,'b','LineWidth',2);
hold all;
plot(tAx,pffOnMeanAlong,'r','LineWidth',2);
legend('PFF Off','PFF On')
xlabel('Time [ns]')
ylabel('Phase [degrees]')
xlim(plotRange)
title('Phase Along Pulse')

figure;
plot(tAx,pffOffStdAlong,'b','LineWidth',2);
hold all;
plot(tAx,pffOnStdAlong,'r','LineWidth',2);
plot(get(gca,'XLim'),[0.2 0.2],'k','LineWidth',2)
plot(get(gca,'XLim'),[0.3 0.3],'k--','LineWidth',1)
legend('PFF Off','PFF On')
xlabel('Time [ns]')
ylabel('Phase Jitter [degrees]')
xlim(plotRange)
title('Phase Jitter Along Pulse')