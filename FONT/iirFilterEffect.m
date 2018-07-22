IIROff = load('/home/jack/PhaseFeedforward/Analysis/201511/20151118_1732_Resolution_NoIIR_181115.mat');
IIROn = load('/home/jack/PhaseFeedforward/Analysis/201511/20151118_1748_Resolution_80MHzLowPass_181115.mat');
saveDir = '/home/jack/Documents/FONT/Droop/201511';
%%
figure;
plot(IIROff.meanDiodeAlongPulse(1,:),'b','LineWidth',2);
hold all;
plot(IIROn.meanDiodeAlongPulse(1,:),'r','LineWidth',2);
xlabel('Time [Sample No.]')
ylabel('Output [counts]')
title('FONT5a Mon1 Diode');
xlim([300 900])
legend('Filter Off', 'Filter On','Location','South');
format_plots;
savePlot(saveDir,'diodeFiltOffOn');

figure;
plot(IIROff.meanMixerAlongPulse(1,:),'b','LineWidth',2);
hold all;
plot(IIROn.meanMixerAlongPulse(1,:),'r','LineWidth',2);
xlabel('Time [Sample No.]')
ylabel('Output [counts]')
title('FONT5a Mon1 Mixer');
xlim([300 900])
ylim([-2000 2000])
legend('Filter Off', 'Filter On');
format_plots;
savePlot(saveDir,'mixerFiltOffOn');

figure;
plot(IIROff.meanPhaseAlongPulse(1,:),'b','LineWidth',2);
hold all;
plot(IIROn.meanPhaseAlongPulse(1,:),'r','LineWidth',2);
xlabel('Time [Sample No.]')
ylabel('Phase [degrees]')
title('FONT5a Mon1 Phase');
xlim([300 900])
ylim([-50 20])
legend('Filter Off', 'Filter On');
format_plots;
savePlot(saveDir,'phaseFiltOffOn');

sampRange = 453:844; %453:844; % 472:835
y = IIROff.meanDiodeAlongPulse(1,sampRange);
n = length(y);
% x = 1:n;
fitresult = expPlusCFit(sampRange,y,ones(1,n));
yFit = fitresult.a*exp(-sampRange./fitresult.b)+fitresult.c;

figure;
plot(sampRange,y,'b','LineWidth',2);
hold all;
plot(sampRange,yFit,'r','LineWidth',2);
legend('DATA','FIT')
xlim([450 850])
xlabel('Time [Sample No.]')
ylabel('Output [counts]')
title('FONT5a Mon1 Diode Fit');
format_plots;
savePlot(saveDir,'diodeFit');

figure;
plot(sampRange,y-yFit,'r','LineWidth',2);
xlim([450 850])
xlabel('Time [Sample No.]')
ylabel('Residual [counts]')
title('FONT5a Mon1 Diode Fit Residuals');
format_plots;
savePlot(saveDir,'diodeFitResid');


