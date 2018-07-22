%%
% load R56 scan data first

%%
DX = 0.61; % dispersion at CT.608
R56TL2 = -0.17;%-0.1; %-0.17% r56 TL2
T566 = -13.5; % t566 between upstream-downstream
ENOffset = -2e-3; % mean energy offset for simulation
simShift = 0;%+0.055; % shift simulation results along x-axis (make min jitter at same location, probably better to do this by adjusting other parameters)
%%
residR56Mod = (R56TL2 + dataSetValues)*(2*pi/0.025);
t566Mod = T566*(2*pi/0.025);

rhoUP = NaN(1,nDataSets);
rhoUDMeas = NaN(1,nDataSets);
sigU = NaN(1,nDataSets);
sigP = NaN(1,nDataSets);
sigDMeas = NaN(1,nDataSets);
rhoUP_err = NaN(1,nDataSets);
rhoUDMeas_err = NaN(1,nDataSets);
sigU_err = NaN(1,nDataSets);
sigP_err = NaN(1,nDataSets);
sigDMeas_err = NaN(1,nDataSets);
sigCCBPM = NaN(1,nDataSets);
sigCCBPM_err = NaN(1,nDataSets);
meanCCBPM = NaN(1,nDataSets);
meanCCBPM_err = NaN(1,nDataSets);
sigD_t566 = NaN(1,nDataSets);

stdDAlong = NaN(nDataSets,data{1}.nSamples);

for ds=1:nDataSets
    rhoUP(ds) = data{ds}.corrMon2_BPMH(data{ds}.refBPM2Index);
    rhoUDMeas(ds) = data{ds}.corrMeanMix2Mix3;
    sigP(ds) = data{ds}.stdMeanBPMH{data{ds}.refBPM2Index};
    sigU(ds) = data{ds}.stdMeanPulsePhase(2);
    sigDMeas(ds) = data{ds}.stdMeanPulsePhase(3);
    sigCCBPM(ds) = data{ds}.stdMeanBPMS{42};
    meanCCBPM(ds) = nanmean(data{ds}.meanBPMS{42});
    
    rhoUP_err(ds) = data{ds}.corrMon2_BPMH_err(data{ds}.refBPM2Index);
    rhoUDMeas_err(ds) = data{ds}.corrMeanMix2Mix3_err;
    sigP_err(ds) = data{ds}.stdMeanBPMH_err{data{ds}.refBPM2Index};
    sigU_err(ds) = data{ds}.stdMeanPulsePhase_err(2);
    sigDMeas_err(ds) = data{ds}.stdMeanPulsePhase_err(3); 
    sigCCBPM_err(ds) = data{ds}.stdMeanBPMS_err{42};
    [~,~,meanCCBPM_err(ds)] = nanMeanStdErr(data{ds}.meanBPMS{42});
    
    p = data{ds}.meanBPMH{data{ds}.refBPM2Index}';
    p = ((p-nanmean(p))./(1000*DX))+ENOffset;
    u = (pi/180)*data{ds}.meanPulsePhase(2,:);
    t566Sim = u + p.*residR56Mod(ds) + (p.^2).*t566Mod;
    sigD_t566(ds) = nanstd(t566Sim);
    tmpCorr = nancorrcoef(u,t566Sim);
    rhoUD_t566(ds) = tmpCorr(1,2);
   
    

end

sigP = sigP./(1000*DX);
sigP_err = sigP_err./(1000*DX);
sigU = sigU.*(pi/180);

rhoUP = rhoUP;

sigDSim = sqrt(sigU.^2 + (residR56Mod.^2).*(sigP.^2) + 2.*residR56Mod.*rhoUP.*sigU.*sigP);
rhoUDSim = (sigU + residR56Mod.*rhoUP.*sigP)./sigDSim;


% sim error calculations
% err^2 = sum(partialDiffs^2*paramErrs^2)
sdErr_su = (0.5./sigDSim).*(2.*sigU + 2.*residR56Mod.*rhoUP.*sigP);
sdErr_su = (sdErr_su.^2).*(sigU_err.^2);
sdErr_sp = (0.5./sigDSim).*(2.*residR56Mod.*sigP + 2.*residR56Mod.*rhoUP.*sigU);
sdErr_sp = (sdErr_sp.^2).*(sigP_err.^2);
sdErr_rup = (0.5./sigDSim).*(2.*residR56Mod.*sigU.*sigP);
sdErr_rup = (sdErr_rup.^2).*(rhoUP_err.^2);
sigDSim_err = sqrt(sdErr_su + sdErr_sp + sdErr_rup);



sigDSim = sigDSim*(180/pi);
sigDSim_err = sigDSim_err*(180/pi);
sigU = sigU*(180/pi);
sigD_t566 = sigD_t566.*(180/pi);


%%

figure;
errorbar(dataSetValues,sigDMeas,sigDMeas_err,'o')
hold all
plot(dataSetValues-simShift,sigDSim)
plot(dataSetValues-simShift,sigD_t566);
legend('meas','r56','t566')
title('jitter down')

figure;
errorbar(dataSetValues,rhoUDMeas,rhoUDMeas_err,'o','MarkerFaceColor','b')
hold all
plot(dataSetValues-simShift,rhoUDSim)
plot(dataSetValues-simShift,rhoUD_t566);
legend('meas','r56','t566')
title('corr UD')

return;
%%
figure;
subplot(4,1,1)
errorbar(dataSetValues,rhoUP,rhoUP_err,'o','MarkerFaceColor','b');
title('Upstream Phase-Energy Correlation')
ylabel('Correlation')
xlabel('R56 in TL1 [m]')
xlim([-0.1 0.6])
ylim([-0.8 0.6])
format_plots;

subplot(4,1,2)
errorbar(dataSetValues,sigP,sigP_err,'o','MarkerFaceColor','r')
title('Relative Energy Jitter')
ylabel('Jitter')
xlabel('R56 in TL1 [m]')
xlim([-0.1 0.6])
ylim([0.4 1.2]*1e-3)
format_plots;

subplot(4,1,3)
errorbar(dataSetValues,sigU,sigU_err,'o','MarkerFaceColor',[0 0.8 0]);
title('Upstream Phase Jitter')
ylabel('Jitter [degrees]')
xlabel('R56 in TL1 [m]')
xlim([-0.1 0.6])
ylim([0.3 1.0])
format_plots;

% subplot(5,1,4)
% errorbar(dataSetValues,meanCCBPM,meanCCBPM_err,'o','MarkerFaceColor','m');
% title('Downstream Beam Transmission')
% ylabel('Beam Current [A]')
% xlabel('R56 in TL1 [m]')
% xlim([-0.1 0.6])
% % ylim([0.3 1.0])
% format_plots;

subplot(4,1,4)
errorbar(dataSetValues,sigCCBPM,sigCCBPM_err,'o','MarkerFaceColor','k');
title('Downstream Beam Current Stability')
ylabel('Jitter [A]')
xlabel('R56 in TL1 [m]')
xlim([-0.1 0.6])
ylim([0.01 0.17])
format_plots;


%%
%5,11,19
figure;
for ds=[5,12,19]
plot(data{ds}.meanPhaseAlongPulse(3,:)-nanmean(data{ds}.meanPulsePhase(3,:)))
hold all
end
plot(data{11}.meanPhaseAlongPulse(2,:)-nanmean(data{ds}.meanPulsePhase(2,:)),'k')

%%
figure;
for ds=1:nDataSets
    plot(data{ds}.timeAxis,data{ds}.meanPhaseAlongPulse(3,:)-data{6}.meanPhaseAlongPulse(3,:),'LineWidth',lineWidthBig,'Color',datSetColors(ds,:));
    hold all;
end
xlim([data{1}.timeAxis(plotStartSamp) data{1}.timeAxis(plotEndSamp)])
xlabel('Time [ns]');
ylabel('Phase [degrees]')
legLabels = dataSetLabels;
yLim = get(gca,'YLim');
format_compare_legend;
format_plots;

figure;
bpmAlong  = data{13}.meanBPMHAlongPulse{data{13}.refBPM2Index};
phaseAlong = data{end}.meanPhaseAlongPulse(3,:)-data{13}.meanPhaseAlongPulse(3,:);%data{6}.meanPhaseAlongPulse(3,:);
bpmAlong = bpmAlong-nanmean(bpmAlong(data{13}.bpmSampleRange{data{13}.refBPM2Index}));
phaseAlong = phaseAlong-nanmean(phaseAlong(data{13}.sampleRange));
bpmAlong = bpmAlong./max(abs(bpmAlong(data{13}.bpmSampleRange{data{13}.refBPM2Index})));
phaseAlong = phaseAlong./max(abs(phaseAlong(data{13}.sampleRange)));
bpmAlong = delaySignal(bpmAlong,2);
plot(phaseAlong,'b','LineWidth',2);
hold all;
plot(bpmAlong,'r','LineWidth',2);
xlim([540 710]);
xlabel('Sample No.');
ylabel('Output [a.u.]')
legend('Phase','Energy')
title({'Phase and Energy Along Pulse' '(phase is diff between R56=0.3 and R56=opt)'})
format_plots;
% savePlot(plotDir,'r56Scan_comparisonPhaseEnergy')

figure;
plot(data{6}.meanPhaseAlongPulse(2,:));
hold all;
plot(data{6}.meanPhaseAlongPulse(3,:));
plot(data{3}.meanPhaseAlongPulse(3,:));
plot(data{end}.meanPhaseAlongPulse(3,:));

hold all;
plot(data{6}.timeAxis,data{6}.meanPhaseAlongPulse(2,:),'k--','LineWidth',2);

%%
figure;
errorbar(scan1Vals,scan1,scan1Err,'o');
hold all;

figure;
errorbar(scan2Vals,scan2,scan2Err,'o');
figure;
errorbar(scan3Vals,scan3,scan3Err,'o');

