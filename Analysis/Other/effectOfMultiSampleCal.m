%%
sampRange = 500:710;
load('/home/jack/PhaseFeedforward/Analysis/201511/20151116_1750_InitConds.mat')
load('/home/jack/PhaseFeedforward/CTFData/201511/FrascatiCalibrations/frascatiCalibration_20151116_1651_calVsSample.mat');

mixers = mixers(:,:,sampRange);
diodes = diodes(:,:,sampRange);
meanPhaseAlongPulse = meanPhaseAlongPulse(:,sampRange);
stdPhaseAlongPulse = stdPhaseAlongPulse(:,sampRange);
timeAxis = timeAxis(sampRange);

A = calFactors(:,:,1);
b = calFactors(:,:,2);
c = calFactors(:,:,3);
d = calFactors(:,:,4);

[nM,nP,nS] = size(mixers);
newPhases = NaN(nM,nP,nS);
for m=1:nM
    for p=1:nP
        for s=1:nS
            tmp(m,p,s) = (mixers(m,p,s)-d(m,s))./A(m,s);
            newPhases(m,p,s) = real(asin((mixers(m,p,s)-d(m,s))./A(m,s)));
%             newPhases(m,p,s) = (newPhases(m,p,s)-c(m,s))./b(m,s);
        end
    end
end
newPhases = newPhases*(180/pi);

[newMeanAlong,newStdAlong] = nanMeanStdErr(newPhases,2);
newMeanAlong(1,:) = newMeanAlong(1,:)-newMeanAlong(1,150);
newMeanAlong(2,:) = newMeanAlong(2,:)-newMeanAlong(2,150);
newMeanAlong(3,:) = newMeanAlong(3,:)-newMeanAlong(3,150);

[~,newRes,~,newRes_err] = nanMeanStdErr(squeeze(newPhases(1,:,:))-squeeze(newPhases(2,:,:)));
newRes = newRes./sqrt(2);
newRes_err = newRes_err./sqrt(2);

figure;
subplot(3,1,1:2)
plot(timeAxis,meanPhaseAlongPulse(2,:)-meanPhaseAlongPulse(2,150),'b','LineWidth',2);
hold all;
plot(timeAxis,newMeanAlong(2,:),'r','LineWidth',2)
title('Upstream Phase (Mon 2)')
ylabel('Phase [degrees]')
% xlabel('Time [ns]')
xlim([40 1094])
ylim([-5 40])
legend('Single Sample','Multi-Sample')
format_plots
subplot(3,1,3)
plot(timeAxis,newMeanAlong(2,:)-meanPhaseAlongPulse(2,:)+meanPhaseAlongPulse(2,150),'k','LineWidth',2);
title('Difference')
ylabel('Diff [degs]')
xlabel('Time [ns]')
xlim([40 1094])
ylim([-1.5 0.5])
format_plots

figure;
subplot(3,1,1:2)
plot(timeAxis,meanPhaseAlongPulse(3,:)-meanPhaseAlongPulse(3,150),'b','LineWidth',2);
hold all;
plot(timeAxis,newMeanAlong(3,:),'r','LineWidth',2)
title('Downstream Phase (Mon 3)')
ylabel('Phase [degrees]')
% xlabel('Time [ns]')
xlim([40 1094])
ylim([-5 40])
legend('Single Sample','Multi-Sample')
format_plots
subplot(3,1,3)
plot(timeAxis,newMeanAlong(3,:)-meanPhaseAlongPulse(3,:)+meanPhaseAlongPulse(3,150),'k','LineWidth',2);
title('Difference')
ylabel('Diff [degs]')
xlabel('Time [ns]')
xlim([40 1094])
ylim([-1.5 0.5])
format_plots


% 
% figure;
% plot(stdPhaseAlongPulse(2,:),'b');
% hold all;
% plot(newStdAlong(2,:),'r')
% 
% figure;
% plot(stdPhaseAlongPulse(3,:),'b');
% hold all;
% plot(newStdAlong(3,:),'r')
% 
% figure;
% plot(newMeanAlong(1,:));
% hold all
% plot(newMeanAlong(2,:));
% plot(newMeanAlong(3,:));

