%dataFile = '/home/jack/PhaseFeedforward/Analysis/CorrelationAnalysis/R56/20141215_1411_allForStability_R56_+0.3/correlationAnalysis.mat';%
%dataFile = '/home/jack/PhaseFeedforward/Analysis/CorrelationAnalysis/R56/20141205_1330_R56_+0.3/correlationAnalysis.mat';
%dataFile = '/home/jack/PhaseFeedforward/Analysis/CorrelationAnalysis/R56/20141205_1125_R56_-0.2_Gain26/correlationAnalysis.mat';
%dataFile = '/home/jack/PhaseFeedforward/Analysis/CorrelationAnalysis/R56/20141205_1125_R56_-0.2_Gain-26/correlationAnalysis.mat';
%dataFile = '/home/jack/PhaseFeedforward/Analysis/CorrelationAnalysis/R56/20141205_1125_R56_+0.4/correlationAnalysis.mat';
%dataFile = '/home/jack/PhaseFeedforward/Analysis/CorrelationAnalysis/R56/20141205_1125_R56_+0.5/correlationAnalysis.mat';
%dataFile = '/home/jack/PhaseFeedforward/Analysis/CorrelationAnalysis/Wiggle/20141216_1915_Gain0_Wiggle23Average/correlationAnalysis.mat';
dataFile = '/home/jack/PhaseFeedforward/Analysis/CorrelationAnalysis/Normal/20141215_1823_FFGain_0/correlationAnalysis.mat';

dataStruct = load(dataFile);

%residR56 = 0.11;
dispersion = 0.17;

saveDir = '/home/jack/PhaseFeedforward/Analysis/Simulation/OptimumPhaseEnergyMix/';

%%
close all

dataName = fileparts(dataFile);
dataName = strsplit(dataName,'/');
dataName = dataName{end};


%%
stdPhaseCT = nanstd(dataStruct.meanMon2Phase);%1.2624*(pi/180);
stdPhaseCB = nanstd(dataStruct.meanMon3Phase);%2.4962*(pi/180);
bpmCT285Index = ismember(dataStruct.bpmNames,'CT.SVBPM0285');
stdPosition = nanstd(dataStruct.meanBPMH(bpmCT285Index,:));%0.1263/1000;

stdPhaseCT = stdPhaseCT*(pi/180);
stdPhaseCB = stdPhaseCB*(pi/180);
stdPosition = stdPosition/1000;

residR56 = (dispersion^2)*((stdPhaseCB^2)-(stdPhaseCT^2));
residR56 = residR56/((stdPosition^2)*(((2*pi*12e9)/(3e8))^2));
residR56 = sqrt(residR56);

%%



% calculate what optimum gain should be
optiGain = (residR56*2*pi*12e9)/(dispersion*3e8);

phaseCT = dataStruct.meanMon2Phase-nanmean(dataStruct.meanMon2Phase);
phaseCB = dataStruct.meanMon3Phase-nanmean(dataStruct.meanMon3Phase);
bpmCT285Index = find(ismember(dataStruct.bpmNames,'CT.SVBPM0285'));
positionCT = dataStruct.meanBPMH(bpmCT285Index,:)-nanmean(dataStruct.meanBPMH(bpmCT285Index,:));

phaseCT = phaseCT';
phaseCB = phaseCB';
positionCT = positionCT./1000;

isGood = ~(isnan(phaseCT) | isnan(phaseCB) | isnan(positionCT));
phaseCT = phaseCT(isGood);
phaseCB = phaseCB(isGood);
positionCT = positionCT(isGood);

corrPhaseCTCB = corrcoef(phaseCT,phaseCB);
corrPhaseCTCB = corrPhaseCTCB(1,2);
signPhaseCorr = sign(corrPhaseCTCB);

corrPosCTPhaseCB = corrcoef(phaseCB,positionCT);
corrPosCTPhaseCB = corrPosCTPhaseCB(1,2);
signENCorr = sign(corrPosCTPhaseCB);

%theta = 0:0.01:(pi/2);
gainEnergy = 0:1:50000;%0:1:50000;%cos(theta);
gainPhase = ones(1,length(gainEnergy));%ones(1,length(gainEnergy));%sin(theta);ones(1,length(gainEnergy))

ffInput = NaN(length(gainPhase),length(phaseCT));
corrInputPhase = NaN(1,length(gainPhase));

for i=1:length(gainEnergy)
    ffInput(i,:) = signPhaseCorr.*gainPhase(i).*phaseCT + signENCorr.*gainEnergy(i).*positionCT;%(positionCT./dispCT).*((residR56*2*pi*12e9)./(3e8));
    %ffInput(i,:) = corrPhaseCTCB.*gainPhase(i).*phaseCT + corrPosCTPhaseCB.*optiGain.*gainEnergy(i).*positionCT;
    %isGood = ~(isnan(ffInput(i,:)) | isnan(phaseCB));
    %tmpCorr = corrcoef(ffInput(i,isGood),phaseCB(isGood));
    tmpCorr = corrcoef(ffInput(i,:),phaseCB);
    corrInputPhase(i) = tmpCorr(1,2);
end

[maxCorr, maxCorrInd] = max(corrInputPhase);

figure;
set(gcf,'PaperPositionMode','auto'); % forces figures to be saved with same aspect ratio
set(gcf,'units','normalized','position',[0.3181 0.2500 0.4506 0.6311])
plot(gainEnergy,corrInputPhase,'LineWidth',2.5);
hold all
%plot([min(fracEnergy) max(fracEnergy)],[dataFile.corrMon2Mon3 dataFile.corrMon2Mon3],'k');
%plot([min(fracEnergy) max(fracEnergy)],[abs(dataFile.corrMon3_BPMH(bpmCT285Index)) abs(dataFile.corrMon3_BPMH(bpmCT285Index))],'k');
plot([min(gainEnergy) max(gainEnergy)],[abs(corrPhaseCTCB) abs(corrPhaseCTCB)],'r','LineWidth',1.5);
plot([min(gainEnergy) max(gainEnergy)],[abs(corrPosCTPhaseCB) abs(corrPosCTPhaseCB)],'Color',[0 0.7 0],'LineWidth',1.5);
plot([min(gainEnergy) max(gainEnergy)],[maxCorr maxCorr],'k','LineWidth',1.5)
legend( 'FF Input - CB Phase',...
        sprintf('CT - CB Phase (%.2f)',corrPhaseCTCB),...
        sprintf('CT Position - CB Phase (%.2f)',corrPosCTPhaseCB),...
        sprintf('Max correlation (%.2f at gain %.2d)',maxCorr,gainEnergy(maxCorrInd)),...
        'Location','best');
grid  on;
xlabel('Energy Gain in FF Input')
ylabel('Correlation')
title({strrep(dataName,'_',' '),sprintf('(residual R56 = %.2f m)',residR56)})
yAxLim = get(gca,'YLim');
plot([gainEnergy(maxCorrInd) gainEnergy(maxCorrInd)],yAxLim,'k','LineWidth',1.5)
ylim(yAxLim);
set(gca,'FontSize',14)
set(findall(gcf,'type','text'),'FontSize',14)
saveStr = [saveDir dataName];
print([saveStr '.png'],'-dpng');
savefig([saveStr '.fig']);
