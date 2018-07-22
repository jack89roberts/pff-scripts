load('/home/jack/PhaseFeedforward/Analysis/201511/20151116_1400_Resolution_0dBPhMonInputs.mat')

samps = 500:715;
res = resolution12(samps);
phas = meanPhaseAlongPulse(2,samps)-meanPhaseAlongPulse(2,637);
nBins = 20;

binRes = NaN(1,nBins);
binPhas = NaN(1,nBins);
binResErr = NaN(1,nBins);
binPhasErr = NaN(1,nBins);
nBin = NaN(1,nBins);

maxPhas = max(phas);
minPhas = min(phas);
binInt = (maxPhas-minPhas)/nBins;

binStart = minPhas+((0:(nBins-1))*binInt);
binEnd=minPhas+((1:nBins)*binInt);

for i=1:nBins
    binSamps = (phas>=binStart(i)) & (phas<binEnd(i));
    nBin(i) = sum(binSamps);
    
    [binRes(i),~,binResErr(i),~] = nanMeanStdErr(res(binSamps));
    [binPhas(i),~,binPhasErr(i),~] = nanMeanStdErr(phas(binSamps));

end

figure;
errorbar(binPhas-binPhas(9),binRes,binResErr,'bo','LineWidth',2,'MarkerFaceColor','b');
% hold all;
% herrorbar(binPhas,binRes,binPhasErr,'ro');
xlabel('Phase [degrees]')
ylabel('Resolution [degrees]')
title('Resolution Dependence on Small Phase Offsets')
format_plots

%% with sample averaging
mon1=squeeze(phases(1,:,:));
mon2=squeeze(phases(2,:,:));

% f1=figure;
% subplot(1,2,1)
% plot(mon1')
% f2=figure;
% subplot(1,2,1)
% plot(mon2')

figure;
a=shadedErrorBar(timeAxis,resolution12,resolution12_err,'b');
hold all;

nAvg=5;
%[avMon1,midSamps]=averageSamples(mon1,nAvg);
[avMon1,midSamps]=averageSamples(mon1,nAvg,'box');
%avMon2=averageSamples(mon2,nAvg);
avMon2=averageSamples(mon2,nAvg,'box');
[~,avRes,~,avRes_err] = nanMeanStdErr(avMon1-avMon2);
avRes = avRes./sqrt(2);
avRes_err = avRes_err./sqrt(2);

%b=shadedErrorBar(timeAxis,avRes,avRes_err,'r');
b=shadedErrorBar(timeAxis(midSamps),avRes,avRes_err,'r');

% [meanAvRes,meanAvResErr]=combineMeas(avRes(521:713),avRes_err(521:713))
[meanAvRes,meanAvResErr]=combineMeas(avRes(104:143),avRes_err(104:143))

xlim([70 950]);
ylim([0 0.3]);

legend([a.mainLine b.mainLine c.mainLine],'Original','5 Samples Averaged','50 Samples Averaged')
xlabel('Time [ns]')
ylabel('Resolution [degrees]')
title('Resolution with Sample Averaging')
format_plots


% figure(f1)
% subplot(1,2,2)
% plot(mon1');
% figure(f2);
% subplot(1,2,2)
% plot(mon2')
% 

%% resolution on mean phase

figure;
plot(meanPulsePhase(1,:)-meanPulsePhase(2,:))

[~,meanRes,~,meanResErr]=nanMeanStdErr((meanPulsePhase(1,:)-meanPulsePhase(2,:))./sqrt(2))


%%
sampRange = 500:710;
load('/home/jack/PhaseFeedforward/Analysis/201511/20151116_1400_Resolution_0dBPhMonInputs.mat')
load('/home/jack/PhaseFeedforward/CTFData/201511/FrascatiCalibrations/frascatiCalibration_20151116_1334_calVsSample.mat');

mixers = mixers(:,:,sampRange);
diodes = diodes(:,:,sampRange);

A = calFactors(:,:,1);
b = calFactors(:,:,2);
c = calFactors(:,:,3);
d = calFactors(:,:,4);

[nM,nP,nS] = size(mixers);
newPhases = NaN(nM,nP,nS);

for m=1:nM
    for p=1:nP
        for s=1:nS
            newPhases(m,p,s) = asin((mixers(m,p,s)-d(m,s))./A(m,s));
%             newPhases(m,p,s) = (newPhases(m,p,s)-c(m,s))./b(m,s);
        end
    end
end
newPhases = newPhases*(180/pi);

newMeanAlong = squeeze(nanmean(newPhases,2));
newMeanAlong(1,:) = newMeanAlong(1,:)-newMeanAlong(1,150);
newMeanAlong(2,:) = newMeanAlong(2,:)-newMeanAlong(2,150);
newMeanAlong(3,:) = newMeanAlong(3,:)-newMeanAlong(3,150);

[~,newRes,~,newRes_err] = nanMeanStdErr(squeeze(newPhases(1,:,:))-squeeze(newPhases(2,:,:)));
newRes = newRes./sqrt(2);
newRes_err = newRes_err./sqrt(2);


figure;
plot(newMeanAlong(1,:));
hold all
plot(newMeanAlong(2,:));

% figure;
% plot(squeeze(newPhases(1,:,:))','b')
% hold all;
% plot(squeeze(newPhases(2,:,:))','r')
% 
% figure;
% plot(squeeze(phases(1,:,sampRange))','b')
% hold all;
% plot(squeeze(phases(2,:,sampRange))','r')

figure;
a=shadedErrorBar(timeAxis(sampRange),resolution12(sampRange),resolution12_err(sampRange),'b',0.5);
hold all;
b=shadedErrorBar(timeAxis(sampRange),newRes,newRes_err,'r',0.5);
xlim([200 1200])
legend([a.mainLine b.mainLine],'Single Sample','Multi-Sample')
xlabel('Time [ns]');
ylabel('Resolution [degrees]')
title('Resolution Dependence on Calibration Method')
format_plots

[origR,origR_err]=combineMeas(resolution12(sampRange),resolution12_err(sampRange))
[newR,newR_err]=combineMeas(newRes,newRes_err)