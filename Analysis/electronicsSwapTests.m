close all;

dataDir = '/home/jack/PhaseFeedforward/Analysis/201504/PhaseMonitors';
saveDir = '/home/jack/PhaseFeedforward/Analysis/201504/PhaseMonitors/Plots/Thesis';
saveData = 0;

dataSetNames = {...
'20150421_1152_Mix1Mon1_Mix2Mon2_10dbThenAmp'
'20150421_1401_Mix1Mon2_Mix2Mon1_NoAmp'
'20150421_1534_Mix1Mon1_Mix3Mon2_NoAmp'
'20150422_1031_Mix1Mon1_Mix2Mon2_NoAmp_PhShiftersSwapped'
'20150422_1331_Mix1Mon1PhShft3_Mix3Mon2PhShft1_NoAmp'
};

plotMix = ...
[1 1 0; 1 1 0; 1 0 1; 1 1 0; 1 0 1];

nDatSets = length(dataSetNames);

for i=1:nDatSets
    data = load([dataDir '/' dataSetNames{i} '.mat']);
    hands = [];
    figure;
    if (plotMix(i,1))
        tmp=shadedErrorBar(data.timeAxis,data.stdPhaseAlongPulse(1,:),data.stdPhaseAlongPulse_err(1,:),'b');
        hands = [hands tmp.mainLine];
        hold all;
    end
    if (plotMix(i,2))
        tmp=shadedErrorBar(data.timeAxis,data.stdPhaseAlongPulse(2,:),data.stdPhaseAlongPulse_err(2,:),'r');
        hands = [hands tmp.mainLine];
        hold all;
    end
    
    if (plotMix(i,3))
        tmp=shadedErrorBar(data.timeAxis,data.stdPhaseAlongPulse(3,:),data.stdPhaseAlongPulse_err(3,:),'g');
        hands = [hands tmp.mainLine];
        hold all;
    end
    xlim([30 224])
    ylim([0.5 2.5])
    xlabel('Time [ns]');
    ylabel('Jitter [degrees]')
    splitName = strsplit(dataSetNames{i},'_');
    if i==4
        legend(hands,'Mix1Mon1PhShft2','Mix2Mon2PhShft1')
    else
        legend(hands,splitName{3},splitName{4})
    end
    title([splitName{3} ' and ' splitName{4}])
    format_plots
    
    if (saveData)
        if i==4
            savePlot(saveDir,'Mix1Mon1PhShft2_Mix2Mon2PhShft1');
        else
            savePlot(saveDir,[splitName{3} '_' splitName{4}]);
        end
    end
    
    [~,statsStart] = min(abs(data.timeAxis-30));
    [~,statsEnd] = min(abs(data.timeAxis-224));
    statsRange = statsStart:statsEnd;
    [mean1,err1]=combineMeas(data.stdPhaseAlongPulse(1,statsRange),data.stdPhaseAlongPulse_err(1,statsRange));
    [mean2,err2]=combineMeas(data.stdPhaseAlongPulse(2,statsRange),data.stdPhaseAlongPulse_err(2,statsRange));
    [mean3,err3]=combineMeas(data.stdPhaseAlongPulse(3,statsRange),data.stdPhaseAlongPulse_err(3,statsRange));
    
    fprintf('%s: Mix1=%.2f%c%.2f, Mix2=%.2f%c%.2f,Mix3=%.2f%c%.2f\n',[splitName{3} '_' splitName{4}],mean1,char(177),err1,mean2,char(177),err2,mean3,char(177),err3);
    
    
end
