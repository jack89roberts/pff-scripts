clearvars; close all;
%%

datDir = '/home/jack/PhaseFeedforward/CTFData/201610/20161025_1358_AcqCheck';
saveDir = '/home/jack/PhaseFeedforward/Analysis/201610/Plots/20161025_1358_AcqCheck';
savePlots = 0;

figSize = [66  1  1535  821];

%%

myFiles = dir([datDir '/*.mat']);

nPulses = length(myFiles);

tPFFBridgeAcq = NaN(1,nPulses);
tPFFBridgeLV = NaN(1,nPulses);
tScope = NaN(1,nPulses);
tBPM = NaN(1,nPulses);

adc1 = NaN(nPulses,900);
scopeCH3 = NaN(nPulses,768);
errorADC = zeros(1,nPulses);
errorScope = zeros(1,nPulses);
errorBPM = zeros(1,nPulses);
for i=1:nPulses
    fprintf('File %d out of %d: %s\n',i,nPulses,myFiles(i).name);
    load([datDir '/' myFiles(i).name]);
    
    try
        tPFFBridgeAcq(i) = myDataStruct.CT_PFFBRIDGE1.Acquisition.timeStamp;
        tPFFBridgeLV(i) = myDataStruct.CT_PFFBRIDGE1.Acquisition.value.timeStampPFF;
        adc1(i,:) = myDataStruct.CT_PFFBRIDGE1.Acquisition.value.ADC1;
        adc1(i,:) = adc1(i,:)-nanmean(adc1(i,1:10));
    catch
        errorADC(i) = 1;
    end
    try  
        tScope(i) = myDataStruct.CT_SCOPE01_CH01.Acquisition.value.timeStamp;

        scopeCH3(i,:) = myDataStruct.CT_SCOPE01_CH03.Acquisition.value.value;
        scopeCH3(i,:) = scopeCH3(i,:)-nanmean(scopeCH3(i,1:10));
    catch
        errorScope(i) = 1;
    end
    try
        tBPM(i) = myDataStruct.CC_SVBPM0845S.Samples.samples.timeStamp;
    catch
        errorBPM(i) = 1;
    end
    
end
%%
adcErrors = find(errorADC)
[~,dupADC] = removeDuplicatePulses(adc1);
adcDuplicates = find(dupADC)

scopeErrors = find(errorScope)
[~,dupScope] = removeDuplicatePulses(scopeCH3);
scopeDuplicates = find(dupScope)

bpmErrors = find(errorBPM)
[~,dupBPM] = removeDuplicatePulses(tBPM);
bpmDuplicates = find(dupBPM)

%%
figure;
plot(adc1');
figure;
plot(scopeCH3');

histADC1 = adc1(:,200);
histScope = scopeCH3(:,300);
histADC1 = abs(histADC1);
histScope=  abs(histScope);
histADC1 = histADC1./max(histADC1);
histScope = histScope./max(histScope);

figure;
plot(histADC1,'b','LineWidth',1);
hold all;
plot(histScope,'r','LineWidth',1);
xlabel('Pulse No.')
ylabel('Output [a.u.]')
title('Synchronisation of FONT and SiS Data')
legend('FONT','SiS')
set(gcf,'Position',figSize);
if savePlots; savePlot(saveDir,'syncFONTSiS'); end;
%%
figure;
plot(diff(tPFFBridgeAcq))
title('PFF Bridge Diff Acq Time Stamp')
ylabel('Time Diff [ns]')
xlabel('Pulse No.')
set(gcf,'Position',figSize);
if savePlots; savePlot(saveDir,'tPFFBridgeAcq'); end;

figure;
plot(diff(tPFFBridgeLV))
title('PFF Bridge Diff LabVIEW Time Stamp')
ylabel('Time Diff [ms]')
xlabel('Pulse No.')
set(gcf,'Position',figSize);
if savePlots; savePlot(saveDir,'tPFFBridgeLV'); end;

figure;
plot(diff(tScope))
title('Scope Diff Acq Time Stamp')
ylabel('Time Diff [ns]')
xlabel('Pulse No.')
set(gcf,'Position',figSize);
if savePlots; savePlot(saveDir,'tScope'); end;

figure;
plot(diff(tBPM))
title('BPM Diff Acq Time Stamp')
ylabel('Time Diff [ns]')
xlabel('Pulse No.')
set(gcf,'Position',figSize);
if savePlots; savePlot(saveDir,'tBPM'); end;

%%
figure;
plot(diff(tPFFBridgeAcq))
title('PFF Bridge Diff Acq Time Stamp')
ylabel('Time Diff [ns]')
xlabel('Pulse No.')
hold all;
plot([0 nPulses],[2.4e9 2.4e9],'k','LineWidth',2);
plot([0 nPulses],[0 0],'k','LineWidth',2);
plot([0 nPulses],[1.8e9 1.8e9],'k--','LineWidth',1.5);
plot([0 nPulses],[0.6e9 0.6e9],'k--','LineWidth',1.5);
set(gcf,'Position',figSize);
if savePlots; savePlot(saveDir,'tPFFBridgeAcq_withLims'); end;

figure;
ax1=subplot(3,1,1);
plot(diff(tPFFBridgeLV))
title('PFF Bridge Diff LabVIEW Time Stamp')
ylabel('Time Diff [ms]')
xlabel('Pulse No.')
hold all;
plot([0 nPulses],[2400 2400],'k','LineWidth',2);
plot([0 nPulses],[0 0],'k','LineWidth',2);
plot([0 nPulses],[1800 1800],'k--','LineWidth',1.5);
plot([0 nPulses],[600 600],'k--','LineWidth',1.5);
ax2=subplot(3,1,2);
% for i=duplicatePulses
%     plot([i i],[0 2400],'r','LineWidth',0.5)
% end
plot(dupADC,'r');
title('Duplicate Acquisition')
xlabel('Pulse No.')
ax3=subplot(3,1,3);
plot(errorADC,'g')
title('Bad/Empty Acquisition')
xlabel('Pulse No.')
% for i=errorPulses
%     plot([i i],[0 2400],'g','LineWidth',0.5)
%     hold all
% end
linkaxes([ax1,ax2,ax3],'x');
set(gcf,'Position',figSize);
if savePlots; savePlot(saveDir,'tPFFBridgeLV_withLims'); end;

figure;
subplot(3,1,1);
ax1=plot(diff(tScope));
title('Scope Diff Acq Time Stamp')
ylabel('Time Diff [ns]')
xlabel('Pulse No.')
hold all;
plot([0 nPulses],[2.4e9 2.4e9],'k','LineWidth',2);
plot([0 nPulses],[0 0],'k','LineWidth',2);
plot([0 nPulses],[1.8e9 1.8e9],'k--','LineWidth',1.5);
plot([0 nPulses],[0.6e9 0.6e9],'k--','LineWidth',1.5);
ax2=subplot(3,1,2);
plot(dupScope,'r');
title('Duplicate Acquisition')
xlabel('Pulse No.')
ax3=subplot(3,1,3);
plot(errorScope,'g')
title('Bad/Empty Acquisition')
xlabel('Pulse No.')
linkaxes([ax1,ax2,ax3],'x');
set(gcf,'Position',figSize);
if savePlots; savePlot(saveDir,'tScope_withLims'); end;

figure;
plot(diff(tBPM))
title('BPM Diff Acq Time Stamp')
ylabel('Time Diff [ns]')
xlabel('Pulse No.')
hold all;
plot([0 nPulses],[2.4e9 2.4e9],'k','LineWidth',2);
plot([0 nPulses],[0 0],'k','LineWidth',2);
plot([0 nPulses],[1.8e9 1.8e9],'k--','LineWidth',1.5);
plot([0 nPulses],[0.6e9 0.6e9],'k--','LineWidth',1.5);
set(gcf,'Position',figSize);
if savePlots; savePlot(saveDir,'tBPM_withLims'); end;

