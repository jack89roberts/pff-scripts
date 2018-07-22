
dataDir = '/nfs/cs-ccr-nfs6/vol29/share/ctf/data/ctfmod/MatLab/PhaseFeedforward/data/20141119_1700_Kick1_GateScan';
bpmSectionPrefixes = {'CC','CM'}; % which section bpms to include in data

%% load data, extract settings from data set comment, extract bpms from data
CTFData = mergeMatMonData(dataDir);

[settingHeadings, pulseSettings] = extractSettingsFromDataSetComment(CTFData);
startGate = pulseSettings(1,:);
endGate = pulseSettings(2,:);
dacOut = pulseSettings(3,:);

[bpmNames, bpmH_MeansAtCursor, bpmH_Samples, bpmS_MeansAtCursor, bpmS_Samples] = extractSectionBPMsFromCTFData(bpmSectionPrefixes, CTFData);
[nBPMs,nPulses,nSamples] = size(bpmH_Samples);

saveDir = [dataDir '/plots'];
if (~isdir(saveDir))
    mkdir(saveDir);
end

%% group the data by start gate and DAC output

uniqueStartGates = unique(startGate);
uniqueDACOuts = unique(dacOut);

nUniqueStartGates = length(uniqueStartGates);
nUniqueDACOuts = length(uniqueDACOuts);

grouped_bpmH_MeansAtCursor = cell(nUniqueStartGates,nUniqueDACOuts);
grouped_bpmS_MeansAtCursor = cell(nUniqueStartGates,nUniqueDACOuts);
grouped_bpmH_Samples = cell(nUniqueStartGates,nUniqueDACOuts);
grouped_bpmS_Samples = cell(nUniqueStartGates,nUniqueDACOuts);

mean_bpmH_MeansAtCursor = NaN*ones(nUniqueStartGates,nUniqueDACOuts, nBPMs);
mean_bpmS_MeansAtCursor = NaN*ones(nUniqueStartGates,nUniqueDACOuts, nBPMs);
mean_bpmH_Samples = NaN*ones(nUniqueStartGates,nUniqueDACOuts, nBPMs, nSamples);
mean_bpmS_Samples = NaN*ones(nUniqueStartGates,nUniqueDACOuts, nBPMs, nSamples);

for g = 1:nUniqueStartGates
    gate = uniqueStartGates(g);
    for d = 1:nUniqueDACOuts
        dac = uniqueDACOuts(d);
        
        grouped_bpmH_MeansAtCursor{g,d} = bpmH_MeansAtCursor(:, startGate==gate & dacOut==dac);
        grouped_bpmS_MeansAtCursor{g,d} = bpmS_MeansAtCursor(:, startGate==gate & dacOut==dac);
        grouped_bpmH_Samples{g,d} = bpmH_Samples(:, startGate==gate & dacOut==dac, :);
        grouped_bpmS_Samples{g,d} = bpmS_Samples(:, startGate==gate & dacOut==dac, :);
        
        mean_bpmH_MeansAtCursor(g,d,:) = nanmean(grouped_bpmH_MeansAtCursor{g,d},2);
        mean_bpmS_MeansAtCursor(g,d,:) = nanmean(grouped_bpmS_MeansAtCursor{g,d},2);
        mean_bpmH_Samples(g,d,:,:) = squeeze(nanmean(grouped_bpmH_Samples{g,d},2));
        mean_bpmS_Samples(g,d,:,:) = squeeze(nanmean(grouped_bpmS_Samples{g,d},2));
        
    end
    
end

% just difference between 2nd and 1st DAC settings for now (but this was
% orbit difference to 0 kick in first data I took).
diff_bpmH_MeansAtCursor = squeeze(mean_bpmH_MeansAtCursor(:,2,:) - mean_bpmH_MeansAtCursor(:,1,:));
diff_bpmS_MeansAtCursor = squeeze(mean_bpmS_MeansAtCursor(:,2,:) - mean_bpmS_MeansAtCursor(:,1,:));
diff_bpmH_Samples = squeeze(mean_bpmH_Samples(:,2,:,:) - mean_bpmH_Samples(:,1,:,:));
diff_bpmS_Samples = squeeze(mean_bpmS_Samples(:,2,:,:) - mean_bpmS_Samples(:,1,:,:));

%%

bpmToPlot = 11;
for g=1:nUniqueStartGates
    subplot(3,1,1:2)
    toPlotNoKick = grouped_bpmH_Samples{g,1};
    [nBPMs,nPulses,nSamples] = size(toPlotNoKick);
    for p=1:nPulses
        plot(squeeze(toPlotNoKick(bpmToPlot,p,:)),'b');
        hold all;
    end
    toPlotKicked = grouped_bpmH_Samples{g,2};
    for p=1:nPulses
        plot(squeeze(toPlotKicked(bpmToPlot,p,:)),'r');
        hold all;
    end
    hold off;
    xlim([250 450]);
    ylim([-5 5]);
    title(sprintf('START GATE %d, BLUE=NoKick, RED=Kicked',uniqueStartGates(g)));
    xlabel('Sample No.');
    ylabel(sprintf('Position %s [mm]',bpmNames{bpmToPlot}));
    
    subplot(3,1,3);
    plot(squeeze(diff_bpmH_Samples(g,bpmToPlot,:)),'k');
    xlim([250 450]);
    ylim([-3 3]);
    grid on
    title('MEAN DIFFERENCE KICK VS. NO KICK');
    xlabel('Sample No.')
    ylabel(sprintf('Difference %s [mm]',bpmNames{bpmToPlot}));
    
    print([saveDir '/'  ], '-dpng');
    pause;
end