for i=1:length(data)
    
    up = data{i}.stdMeanPulsePhase(2,:);
    down = data{i}.stdMeanPulsePhase(3,:);
    up_err = data{i}.stdMeanPulsePhase_err(2,:);
    down_err = data{i}.stdMeanPulsePhase_err(3,:);
    
    ratio
    
end

%%

bpmCorr2 = NaN(1,length(data));
bpmCorr3 = NaN(1,length(data));
phasCorr = NaN(1,length(data));
bpmCorr2_err = NaN(1,length(data));
bpmCorr3_err = NaN(1,length(data));
phasCorr_err = NaN(1,length(data));

meanBPM = NaN(1,length(data));
stdBPM = NaN(1,length(data));
meanBPM_err = NaN(1,length(data));
stdBPM_err = NaN(1,length(data));

for i=1:length(data)
    
    bpmCorr2(i) = data{i}.corrMon2_BPMH(14);
    bpmCorr3(i) = data{i}.corrMon3_BPMH(14);
    phasCorr(i) = data{i}.corrMeanMix2Mix3;
    
    bpmCorr2_err(i) = data{i}.corrMon2_BPMH_err(14);
    bpmCorr3_err(i) = data{i}.corrMon3_BPMH_err(14);
    phasCorr_err(i) = data{i}.corrMeanMix2Mix3_err;
    
    tmpBPM = data{i}.meanBPMH{14};
    [meanBPM(i),~,meanBPM_err(i),~] = nanMeanStdErr(tmpBPM);
    stdBPM(i) = data{i}.stdMeanBPMH{14};
    stdBPM_err(i) = data{i}.stdMeanBPMH_err{14};

    
%     mean2 = data{i}.meanPulsePhase(2,:);
%     mean3 = data{i}.meanPulsePhase(3,:);
%     diffPhas = mean3-mean2;
%     
%     figure;
%     plot(tmpBPM,diffPhas,'o');
%     title(i)
end


figure;
plot(bpmCorr2,'bo')
hold all;
plot(bpmCorr3,'ro')
title('corr with bpm up and down')

figure;
plot(bpmCorr3-bpmCorr2);
title('diff bpm corr');

figure;
errorbar(bpmCorr3-bpmCorr2,phasCorr,phasCorr_err,'bo','MarkerFaceColor','b');
hold all;
herrorbar(bpmCorr3-bpmCorr2,phasCorr,sqrt(bpmCorr3_err.^2 + bpmCorr2_err.^2),'bo');
title('phas corr vs. diff bpm corr')

figure;
plot(meanBPM,phasCorr,'o');
title('mean bpm vs. phas corr')

figure;
plot(stdBPM,phasCorr,'o');

figure;
plot(phasCorr)

figure;
subplot(2,1,2)
errorbar((stdBPM/0.61)*1e-3,(stdBPM_err/0.6)*1e-3,'bo-','LineWidth',2,'MarkerFaceColor','b')
xlim([0.8 15.2])
set(gca,'XTick',2:2:14);
set(gca,'XTickLabel',dataSetLabels(2:2:14))
ylabel('Energy Jitter')
xlabel('Dataset')
subplot(2,1,1)
errorbar(((meanBPM-meanBPM(2))/-0.61)*1e-3,(sqrt(meanBPM_err.^2 + meanBPM_err(2).^2)./0.61)*1e-3,'bo-','LineWidth',2,'MarkerFaceColor','b')
xlim([0.8 15.2])
ylim([-2.2 2.2]*1e-3)
set(gca,'XTick',2:2:14);
set(gca,'XTickLabel',dataSetLabels(2:2:14))
xlabel('Dataset')
ylabel('Energy Offset')


