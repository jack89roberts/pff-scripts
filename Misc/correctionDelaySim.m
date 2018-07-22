load('/home/jack/PhaseFeedforward/Analysis/201511/20151120_1538_Gain-800_R56_0.1_Interleaaved_Odd.mat')
%%
upPhas = averageSamples(meanPhaseAlongPulse(2,:),5);

addSamps = 615:625;
upPhas(addSamps) = upPhas(addSamps) -0.2*(((-5:5).^2)-25);

optCorr = -upPhas;
optDown = upPhas+optCorr;

delay = 7;
slowDel = [NaN(1,delay) upPhas(1:(end-delay))];

myx = [515.1303 715.7401];
time = 5.2*(1:length(upPhas));
time = time-time(round(myx(1)));

close all;
figure;
subplot(6,1,1:3)
plot(time,upPhas,'b','LineWidth',2);
hold all;
plot(time,slowDel,'r','LineWidth',2)
plot(time([620 620]),get(gca,'YLim'),'b--');
plot(time([635 635]),get(gca,'YLim'),'r--');
xlim(time(round(myx)));
xlabel('Time [ns]')
ylabel('Phase [degrees]')
format_plots;
grid off;
ylim([-5 40])
subplot(6,1,5:6)
plot(time,upPhas-slowDel,'g','LineWidth',2);
hold all;
plot(time([620 620]),get(gca,'YLim'),'b--');
plot(time([635 635]),get(gca,'YLim'),'r--');
xlim(time(round(myx)));
xlabel('Time [ns]')
ylabel('Phase [degrees]')
ylim([-10 10])
format_plots;
grid off;
