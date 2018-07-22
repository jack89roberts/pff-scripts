p= data{1}.meanBPMH{14};
p = p-nanmean(p);
p = p./(1000*0.605);

figure;
plot(p,data{1}.meanPulsePhase(2,:),'o','MarkerFaceColor','b')
xlabel('Relative Energy Offset')
ylabel('Upstream Phase [degrees]')
title('Upstream Phase vs. Energy')
xlim([-0.0075 0.0075])
legend(['correlation: 0.94' char(177) '0.04'])
ylim([-5 4])
format_plots

figure;
for i=1:nDataSets
plot(data{i}.stdBPMHAlongPulse{14})
hold all;
end

figure;
plot((data{1}.meanBPMHAlongPulse{14}-nanmean(data{1}.meanBPMHAlongPulse{14}(550:700)))/(1000*0.6),'b','LineWidth',2);
xlim([550 700])
xlabel('Sample No.');
ylabel('Relative Energy Offset')
title('Beam Energy Along Pulse')
format_plots;

figure;
plot((data{1}.stdBPMHAlongPulse{14})/(1000*0.6),'b','LineWidth',2);
xlim([550 700])
xlabel('Sample No.');
ylabel('Relative Energy Jitter')
title('Energy Jitter Along Pulse')
format_plots;

figure;
hold all;
plot(data{6}.stdPhaseAlongPulse(3,:))
% plot(data{6}.stdPhaseAlongPulse(2,:))

ds=6;
bpm = data{ds}.stdBPMHAlongPulse{14}.^2;
bpm = bpm-nanmean(bpm(data{ds}.sampleRange));
bpm = bpm./max(abs(bpm(data{ds}.sampleRange)));
phas = -data{1}.stdPhaseAlongPulse(3,:)-data{ds}.stdPhaseAlongPulse(3,:);
phas = phas-nanmean(phas(data{ds}.sampleRange));
phas = phas./max(abs(phas(data{ds}.sampleRange)));
figure;
plot(phas,'b','LineWidth',2);
hold all;
plot(bpm);
xlim([550 700])
ylim([-2 2])
xlabel('Sample No.')
ylabel('Output [a.u.]')
legend('Phase Jitter','Energy Jitter')
title('Phase Jitter and Energy Jitter Along Pulse')
format_plots;

%%

meanPhases = NaN(nDataSets,data{1}.nPulses);
meanEnergies = NaN(nDataSets,data{1}.nPulses);
fitCoeffs = NaN(nDataSets,3);
fitErrs = NaN(nDataSets,3);
for ds=1:nDataSets
    meanPhases(ds,:) = -data{ds}.meanPulsePhase(3,:)-nanmean(data{ds}.meanPulsePhase(3,:));
    meanEnergies(ds,:) = (data{ds}.meanBPMH{14}-nanmean(data{ds}.meanBPMH{14}))./(-0.605*1000);
    [fitCoeffs(ds,:),~,tmpConf] = nanpolyfit(meanEnergies(ds,:),meanPhases(ds,:),2);
    fitErrs(ds,:) = (fitCoeffs(ds,:)-tmpConf(1,:))/2;
end

figure;
plot(meanEnergies(1,:),meanPhases(1,:),'ro','MarkerFaceColor','r')
hold all;
plot(meanEnergies(6,:),meanPhases(6,:),'bo','MarkerFaceColor','b')
plot(meanEnergies(end,:),meanPhases(end,:),'go','MarkerFaceColor','g')
legend('R56 = -0.100 m','R56 = 0.075 m', 'R56 = 0.300 m')
plot(linspace(-8e-3,8e-3,100),polyval(fitCoeffs(1,:),linspace(-8e-3,8e-3,100)),'r','LineWidth',2);
plot(linspace(-8e-3,8e-3,100),polyval(fitCoeffs(6,:),linspace(-8e-3,8e-3,100)),'b','LineWidth',2);
plot(linspace(-8e-3,8e-3,100),polyval(fitCoeffs(end,:),linspace(-8e-3,8e-3,100)),'g','LineWidth',2);
xlabel('Relative Energy Offset')
ylabel('Phase [degrees]')
title('Downstream Phase vs. Beam Energy')
xlim([-8e-3 8e-3])
format_plots;

figure;
errorbar(dataSetValues,fitCoeffs(:,1)*(0.025/360),fitErrs(:,1)*(0.025/360),'bo','MarkerFaceColor','b');
title('T566 Fit (Upstream-Downstream)')
ylabel('T566 [m]')
xlabel('R56 in TL1 [m]')
format_plots;

figure;
errorbar(dataSetValues,fitCoeffs(:,2)*(0.025/360),fitErrs(:,2)*(0.025/360),'bo','MarkerFaceColor','b');
title('R56 Fit (Upstream-Downstream)')
ylabel('R56 [m]')
xlabel('R56 in TL1 [m]')
format_plots;

figure;
errorbar(dataSetValues,fitCoeffs(:,3)*(0.025/360),fitErrs(:,3)*(0.025/360),'bo','MarkerFaceColor','b');

myCols = varycolor(nDataSets);
figure;
for ds=1:nDataSets
    plot(meanEnergies(ds,:),meanPhases(ds,:),'o','Color',myCols(ds,:))
    hold all;
    plot(linspace(-8e-3,8e-3,100),polyval(fitCoeffs(ds,:),linspace(-8e-3,8e-3,100)),'Color',myCols(ds,:),'LineWidth',2);
end
