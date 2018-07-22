close all; clearvars;
%%
datDir = '/home/jack/PhaseFeedforward/Analysis/201511';
saveDir = '/home/jack/PhaseFeedforward/Analysis/201511/20151123_otherJitterSources';

bend540Files = {...
'20151123_1804_CT540S_-2'
'20151123_1757_CT540S_-1'
'20151123_1754_CT540S_0' % commented
'20151123_1805_CT540S_0'
'20151123_1755_CT540S_+1'
'20151123_1802_CT540S_+2'
};
bends540Setting = 133;
bends540FilesToFit = 2:6;

bend670Files = {...
'20151123_1810_CT670S_-0.5'
'20151123_1812_CT670S_-0.25'
'20151123_1814_CT670S_0'
% '20151123_1806_CT670S_0' % commented
'20151123_1811_CT670S_+0.25'
'20151123_1808_CT670S_+0.5'
};
bends670Setting = 164;
bends670FilesToFit = 1:5;

septumFiles = {...
'20151123_1741_Septum_-2'
% '20151123_1749_Septum_-1.5' % commented
'20151123_1800_Septum_-1.5' % commented
'20151123_1736_Septum_-1'
'20151123_1751_Septum_0'
% '20151123_1721_Septum_Nominal' % commented
% '20151123_1725_Septum_Nominal' % commented
% '20151123_1731_Septum_Nominal' % commented
'20151123_1734_Septum_+1'
'20151123_1742_Septum_+1.5'
'20151123_1739_Septum_+2'
};
septumSetting = 890;
septumFilesToFit = 2:6;

bend205Files = {...
'20151123_1824_CR205S_-3'
'20151123_1821_CR205S_-2'
'20151123_1818_CR205S_-1'
'20151123_1816_CR205S_0'
'20151123_1825_CR205S_0' % commented
'20151123_1817_CR205S_+1'
'20151123_1820_CR205S_+2'
'20151123_1823_CR205S_+3'
};
bends205Setting = 155.6;
bends205FilesToFit = 3:6;

powerJit = 1e-4;

%%
n540Files = length(bend540Files);
n670Files = length(bend670Files);
nSeptumFiles = length(septumFiles);
n205Files = length(bend205Files);

bend540Data = cell(1,n540Files);
bend540Values = NaN(1,n540Files);
for i=1:n540Files
    bend540Data{i} = load([datDir '/' bend540Files{i} '.mat']);
    tmpVal = strsplit(bend540Files{i},'_');
    tmpVal = tmpVal{end};
    bend540Values(i) = str2double(tmpVal);
end

bend670Data = cell(1,n670Files);
bend670Values = NaN(1,n670Files);
for i=1:n670Files
    bend670Data{i} = load([datDir '/' bend670Files{i} '.mat']);
    tmpVal = strsplit(bend670Files{i},'_');
    tmpVal = tmpVal{end};
    bend670Values(i) = str2double(tmpVal);
end

septumData = cell(1,nSeptumFiles);
septumValues = NaN(1,nSeptumFiles);
for i=1:nSeptumFiles
    septumData{i} = load([datDir '/' septumFiles{i} '.mat']);
    tmpVal = strsplit(septumFiles{i},'_');
    tmpVal = tmpVal{end};
    septumValues(i) = str2double(tmpVal);
end
%septumValues(1) = 0;

bend205Data = cell(1,n205Files);
bend205Values = NaN(1,n205Files);
for i=1:n205Files
    bend205Data{i} = load([datDir '/' bend205Files{i} '.mat']);
    tmpVal = strsplit(bend205Files{i},'_');
    tmpVal = tmpVal{end};
    bend205Values(i) = str2double(tmpVal);
end

%%

meanPhase540 = NaN(1,n540Files);
meanPhase540_err = NaN(1,n540Files);
up540 = NaN(1,n540Files);
up540_err = NaN(1,n540Files);
for i=1:n540Files
    [meanPhase540(i),~,meanPhase540_err(i)] = nanMeanStdErr(bend540Data{i}.meanPulsePhase(3,:));
    [up540(i),~,up540_err(i)] = nanMeanStdErr(bend540Data{i}.meanPulsePhase(2,:));
end
diffPhase540 = meanPhase540-up540;
diffPhase540_err = sqrt(meanPhase540_err.^2 + up540_err.^2);

meanPhase670 = NaN(1,n670Files);
meanPhase670_err = NaN(1,n670Files);
up670 = NaN(1,n670Files);
up670_err = NaN(1,n670Files);
for i=1:n670Files
    [meanPhase670(i),~,meanPhase670_err(i)] = nanMeanStdErr(bend670Data{i}.meanPulsePhase(3,:));
    [up670(i),~,up670_err(i)] = nanMeanStdErr(bend670Data{i}.meanPulsePhase(2,:));
end
diffPhase670 = meanPhase670-up670;
diffPhase670_err = sqrt(meanPhase670_err.^2 + up670_err.^2);

meanPhaseSeptum = NaN(1,nSeptumFiles);
meanPhaseSeptum_err = NaN(1,nSeptumFiles);
upSeptum = NaN(1,nSeptumFiles);
upSeptum_err = NaN(1,nSeptumFiles);
for i=1:nSeptumFiles
    [meanPhaseSeptum(i),~,meanPhaseSeptum_err(i)] = nanMeanStdErr(septumData{i}.meanPulsePhase(3,:));
    [upSeptum(i),~,upSeptum_err(i)] = nanMeanStdErr(septumData{i}.meanPulsePhase(2,:));
end
diffPhaseSeptum = meanPhaseSeptum-upSeptum;
diffPhaseSeptum_err = sqrt(meanPhaseSeptum_err.^2 + upSeptum_err.^2);

meanPhase205 = NaN(1,n205Files);
meanPhase205_err = NaN(1,n205Files);
up205 = NaN(1,n205Files);
up205_err = NaN(1,n205Files);
for i=1:n205Files
    [meanPhase205(i),~,meanPhase205_err(i)] = nanMeanStdErr(bend205Data{i}.meanPulsePhase(3,:));
    [up205(i),~,up205_err(i)] = nanMeanStdErr(bend205Data{i}.meanPulsePhase(2,:));
end
diffPhase205 = meanPhase205-up205;
diffPhase205_err = sqrt(meanPhase205_err.^2 + up205_err.^2);

%%
[fit540,~,fit540_err] = nanpolyfit(bend540Values(bends540FilesToFit),diffPhase540(bends540FilesToFit),1,1./(diffPhase540_err(bends540FilesToFit).^2));
fit540_err = (fit540-fit540_err(:,1)')/2;
jit540 = abs(fit540(1).*bends540Setting.*powerJit);
jit540_err = fit540_err(1).*bends540Setting.*powerJit;
fprintf('540: setting=%.2f, gradient=%.2f, jit=%.3f%c%.3f\n',bends540Setting,fit540(1),jit540,char(177),jit540_err);

[fit670,~,fit670_err] = nanpolyfit(bend670Values(bends670FilesToFit),diffPhase670(bends670FilesToFit),1,1./(diffPhase670_err(bends670FilesToFit).^2));
fit670_err = (fit670-fit670_err(:,1)')/2;
jit670 = abs(fit670(1).*bends670Setting.*powerJit);
jit670_err = fit670_err(1).*bends670Setting.*powerJit;
fprintf('670: setting=%.2f, gradient=%.2f, jit=%.3f%c%.3f\n',bends670Setting,fit670(1),jit670,char(177),jit670_err);

[fitSeptum,~,fitSeptum_err] = nanpolyfit(septumValues(septumFilesToFit),diffPhaseSeptum(septumFilesToFit),1,1./(diffPhaseSeptum_err(septumFilesToFit).^2));
fitSeptum_err = (fitSeptum-fitSeptum_err(:,1)')/2;
jitSeptum = abs(fitSeptum(1).*septumSetting.*powerJit);
jitSeptum_err = fitSeptum_err(1).*septumSetting.*powerJit;
fprintf('Septum: setting=%.2f, gradient=%.2f, jit=%.3f%c%.3f\n',septumSetting,fitSeptum(1),jitSeptum,char(177),jitSeptum_err);

[fit205,~,fit205_err] = nanpolyfit(bend205Values(bends205FilesToFit),diffPhase205(bends205FilesToFit),1,1./(diffPhase205_err(bends205FilesToFit).^2));
fit205_err = (fit205-fit205_err(:,1)')/2;
jit205 = abs(fit205(1).*bends205Setting.*powerJit);
jit205_err = fit205_err(1).*bends205Setting.*powerJit;
fprintf('205: setting=%.2f, gradient=%.2f, jit=%.3f%c%.3f\n',bends205Setting,fit205(1),jit205,char(177),jit205_err);


%%
% figure;
% errorbar(bend540Values,meanPhase540,meanPhase540_err,'o','MarkerFaceColor','b');
% title('540')
% 
% figure;
% errorbar(bend670Values,meanPhase670,meanPhase670_err,'o','MarkerFaceColor','b');
% title('670')
% 
% figure;
% errorbar(septumValues,meanPhaseSeptum,meanPhaseSeptum_err,'o','MarkerFaceColor','b');
% title('Septum')
% 
% figure;
% errorbar(bend205Values,meanPhase205,meanPhase205_err,'o','MarkerFaceColor','b');
% title('205')
%%
% figure;
% plot(bend540Values,up540,'o','MarkerFaceColor','b');
% title('up 540')
% 
% figure;
% plot(bend670Values,up670,'o','MarkerFaceColor','b');
% title('up 670')
% 
% figure;
% plot(septumValues,upSeptum,'o','MarkerFaceColor','b');
% title('up Septum')
% 
% figure;
% plot(bend205Values,up205,'o','MarkerFaceColor','b');
% title('up 205')
%%
figure;
errorbar(bend540Values,diffPhase540,diffPhase540_err,'o','MarkerFaceColor','b');
hold all;
plot(bend540Values,polyval(fit540,bend540Values),'k');
title('Downstream Phase vs. TL1 540/710 Dipoles')
xlabel('Current Offset [A]');
ylabel('Downstream Phase [degrees]');
format_plots;
savePlot(saveDir,'tl1540');

figure;
errorbar(bend670Values,diffPhase670,diffPhase670_err,'o','MarkerFaceColor','b');
hold all;
plot(bend670Values,polyval(fit670,bend670Values),'k');
title('Downstream Phase vs. TL1 630/670 Dipoles')
xlabel('Current Offset [A]');
ylabel('Downstream Phase [degrees]');
format_plots;
savePlot(saveDir,'tl1670');

figure;
errorbar(septumValues,diffPhaseSeptum,diffPhaseSeptum_err,'o','MarkerFaceColor','b');
hold all;
plot(septumValues,polyval(fitSeptum,septumValues),'k');
title('Downstream Phase vs. CR Septum')
xlabel('Current Offset [A]');
ylabel('Downstream Phase [degrees]');
format_plots;
savePlot(saveDir,'crSeptum');

figure;
errorbar(bend205Values,diffPhase205,diffPhase205_err,'o','MarkerFaceColor','b');
hold all;
plot(bend205Values,polyval(fit205,bend205Values),'k');
title('Downstream Phase vs. CR Dipoles')
xlabel('Current Offset [A]');
ylabel('Downstream Phase [degrees]');
format_plots;
savePlot(saveDir,'crBends');
