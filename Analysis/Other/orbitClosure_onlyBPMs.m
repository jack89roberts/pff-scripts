% orbitClosure
%
% Takes constant kick data and makes some orbit closure plots

close all; clearvars;
%% inputs

dataSetNames = {...
'20161208_1559_OrbClos_K1_4000_K2_-4000_Intlv',...
'20161208_1604_OrbClos_K1_-4000_K2_4000_Intlv',...
};
dataDir = '/home/jack/PhaseFeedforward/CTFData/201612';
saveDir = '/home/jack/PhaseFeedforward/Analysis/201612/Plots/20161208_OrbClos';
savePlots = 0;

bpmSampRange = 555:698;

ampInCounts = [3800 -3800];

pointsToFit = 1:2;

bpmNames = {...
    %'CC_SVBPM0275',...
    'CC_SVBPM0365',...
    'CC_SVBPM0435',...
    'CC_SVBPI0535',...
    'CC_SVBPI0645',...
    'CC_SVBPI0685',...
    'CC_SVBPI0735',...
    'CC_SVBPM0845',...
    'CC_SVBPM0930'...
    %'CB_SVBPM0150'...
    %'CB_SVBPS0210'...
    %'CB_SVBPS0250'...
    %'CB_SVBPS0310'...
};
bpmPos = [3.1831 6.2205 10.1413 13.4795 15.5532 17.0922 20.8066 25.3145];

%% load data
addpath('../../');

nDataSets = length(dataSetNames);
nBPMs = length(bpmNames);

ampInVolts = ampInCounts*(2/4096);

% figure;
for ds=1:nDataSets
    load([dataDir '/' dataSetNames{ds} '/Merged/' dataSetNames{ds} '.mat']);
    
    % Initialise arrays
    if (ds==1)
        nPulses = length(CTFData)/2; % divide by 2 - interleaved mode, splitting pulses in to different arrays later
                
        bpmNSamp = length(CTFData(1).CC_SVBPI0645V.Samples.samples.value);
        bpmSampInt = 5.2083;
        kickOnBPMH = NaN(nDataSets,nBPMs,nPulses,bpmNSamp);
        kickOffBPMH = NaN(nDataSets,nBPMs,nPulses,bpmNSamp);
        kickOnBPMV = NaN(nDataSets,nBPMs,nPulses,bpmNSamp);
        kickOffBPMV = NaN(nDataSets,nBPMs,nPulses,bpmNSamp);
        kickOnBPMS = NaN(nDataSets,nBPMs,nPulses,bpmNSamp);
        kickOffBPMS = NaN(nDataSets,nBPMs,nPulses,bpmNSamp);

    end
    
    
    % ===================== BPM ========================================    
    CTFData = removePulsesNoBeam(CTFData);
    
    bpmH = NaN(nBPMs,nPulses*2,bpmNSamp);
    bpmS = NaN(nBPMs,nPulses*2,bpmNSamp);
    bpmV = NaN(nBPMs,nPulses*2,bpmNSamp);
    for b=1:nBPMs
        [bpmH(b,:,:), bpmS(b,:,:), bpmV(b,:,:)] = extractBPMFromCTFData( bpmNames{b}, CTFData);
        [bpmS(b,:,:),others] = getAlignedXCorr(squeeze(bpmS(b,:,:)),'end',{squeeze(bpmH(b,:,:)) squeeze(bpmV(b,:,:))});
        bpmH(b,:,:) = others{1};
        bpmV(b,:,:) = others{2};
        
        [bpmS(b,1:2:end,:),others] = removeBadPulses(squeeze(bpmS(b,1:2:end,:)),bpmSampRange,{squeeze(bpmH(b,1:2:end,:)) squeeze(bpmV(b,1:2:end,:))});
        bpmH(b,1:2:end,:) = others{1};
        bpmV(b,1:2:end,:) = others{2};
        
        [bpmS(b,2:2:end,:),others] = removeBadPulses(squeeze(bpmS(b,2:2:end,:)),bpmSampRange,{squeeze(bpmH(b,2:2:end,:)) squeeze(bpmV(b,2:2:end,:))});
        bpmH(b,2:2:end,:) = others{1};
        bpmV(b,2:2:end,:) = others{2};

    end
    
    % =========== Split in odd even =======================
    kickOnBPMH(ds,:,:,:) = bpmH(:,1:2:end,:);
    kickOnBPMV(ds,:,:,:) = bpmV(:,1:2:end,:);
    kickOnBPMS(ds,:,:,:) = bpmS(:,1:2:end,:);

    kickOffBPMH(ds,:,:,:) = bpmH(:,2:2:end,:);
    kickOffBPMV(ds,:,:,:) = bpmV(:,2:2:end,:);
    kickOffBPMS(ds,:,:,:) = bpmS(:,2:2:end,:);

        
    clear CTFData ampLA ampLB ampRA ampRB bpmH bpmV bpmS;
end


%% BPMs: calculate means etc.

[bpmHOnSampMean,~,bpmHOnSampMean_err] = nanMeanStdErr(kickOnBPMH,3);
[bpmSOnSampMean,~,bpmSOnSampMean_err] = nanMeanStdErr(kickOnBPMS,3);
[bpmVOnSampMean,~,bpmVOnSampMean_err] = nanMeanStdErr(kickOnBPMV,3);
[bpmHOffSampMean,~,bpmHOffSampMean_err] = nanMeanStdErr(kickOffBPMH,3);
[bpmSOffSampMean,~,bpmSOffSampMean_err] = nanMeanStdErr(kickOffBPMS,3);
[bpmVOffSampMean,~,bpmVOffSampMean_err] = nanMeanStdErr(kickOffBPMV,3);

bpmHOnPulseMean = squeeze(nanmean(kickOnBPMH(:,:,:,bpmSampRange),4));
bpmVOnPulseMean = squeeze(nanmean(kickOnBPMV(:,:,:,bpmSampRange),4));
bpmSOnPulseMean = squeeze(nanmean(kickOnBPMS(:,:,:,bpmSampRange),4));
bpmHOffPulseMean = squeeze(nanmean(kickOffBPMH(:,:,:,bpmSampRange),4));
bpmVOffPulseMean = squeeze(nanmean(kickOffBPMV(:,:,:,bpmSampRange),4));
bpmSOffPulseMean = squeeze(nanmean(kickOffBPMS(:,:,:,bpmSampRange),4));

[bpmHOnMean,~,bpmHOnMean_err] = nanMeanStdErr(bpmHOnPulseMean,3);
[bpmVOnMean,~,bpmVOnMean_err] = nanMeanStdErr(bpmVOnPulseMean,3);
[bpmSOnMean,~,bpmSOnMean_err] = nanMeanStdErr(bpmSOnPulseMean,3);
[bpmHOffMean,~,bpmHOffMean_err] = nanMeanStdErr(bpmHOffPulseMean,3);
[bpmVOffMean,~,bpmVOffMean_err] = nanMeanStdErr(bpmVOffPulseMean,3);
[bpmSOffMean,~,bpmSOffMean_err] = nanMeanStdErr(bpmSOffPulseMean,3);

bpmHDiffSampMean = bpmHOnSampMean-bpmHOffSampMean;
bpmHDiffSampMean_err = sqrt(bpmHOnSampMean_err.^2 + bpmHOffSampMean_err.^2);
bpmHDiffPulseMean = bpmHOnPulseMean-bpmHOffPulseMean;
bpmHDiffMean = bpmHOnMean-bpmHOffMean;
bpmHDiffMean_err = sqrt(bpmHOnMean_err.^2 + bpmHOffMean_err.^2);

bpmVDiffSampMean = bpmVOnSampMean-bpmVOffSampMean;
bpmVDiffSampMean_err = sqrt(bpmVOnSampMean_err.^2 + bpmVOffSampMean_err.^2);
bpmVDiffPulseMean = bpmVOnPulseMean-bpmVOffPulseMean;
bpmVDiffMean = bpmVOnMean-bpmVOffMean;
bpmVDiffMean_err = sqrt(bpmVOnMean_err.^2 + bpmVOffMean_err.^2);

bpmSDiffSampMean = bpmSOnSampMean-bpmSOffSampMean;
bpmSDiffSampMean_err = sqrt(bpmSOnSampMean_err.^2 + bpmSOffSampMean_err.^2);
bpmSDiffPulseMean = bpmSOnPulseMean-bpmSOffPulseMean;
bpmSDiffMean = bpmSOnMean-bpmSOffMean;
bpmSDiffMean_err = sqrt(bpmSOnMean_err.^2 + bpmSOffMean_err.^2);


bpmHDiffFit = NaN(nBPMs,2);
bpmHDiffFit_rsq = NaN(1,nBPMs);
bpmHDiffFit_err = NaN(nBPMs,2);

bpmVDiffFit = NaN(nBPMs,2);
bpmVDiffFit_rsq = NaN(1,nBPMs);
bpmVDiffFit_err = NaN(nBPMs,2);

bpmSDiffFit = NaN(nBPMs,2);
bpmSDiffFit_rsq = NaN(1,nBPMs);
bpmSDiffFit_err = NaN(nBPMs,2);
for b=1:nBPMs
    [bpmHDiffFit(b,:),bpmHDiffFit_rsq(b),tmpConf] = nanpolyfit(ampInVolts(pointsToFit),bpmHDiffMean(pointsToFit,b),1,1./bpmHDiffMean_err(pointsToFit,b).^2);
    bpmHDiffFit_err(b,:) = (bpmHDiffFit(b,:)-tmpConf(1,:))/2;
    
    [bpmVDiffFit(b,:),bpmVDiffFit_rsq(b),tmpConf] = nanpolyfit(ampInVolts(pointsToFit),bpmVDiffMean(pointsToFit,b),1,1./bpmVDiffMean_err(pointsToFit,b).^2);
    bpmVDiffFit_err(b,:) = (bpmVDiffFit(b,:)-tmpConf(1,:))/2;
    
    [bpmSDiffFit(b,:),bpmSDiffFit_rsq(b),tmpConf] = nanpolyfit(ampInVolts(pointsToFit),bpmSDiffMean(pointsToFit,b),1,1./bpmSDiffMean_err(pointsToFit,b).^2);
    bpmSDiffFit_err(b,:) = (bpmSDiffFit(b,:)-tmpConf(1,:))/2;
end

%%
xTickLab = strrep(bpmNames,'_','.');
xTickLab = strrep(xTickLab,'SV','');
xTickLab = strrep(xTickLab,'BPM','');
xTickLab = strrep(xTickLab,'BPI','');
xTickLab = strrep(xTickLab,'.0','.');

figure
myCols = varycolor(nDataSets);
for i=1:nDataSets
    %shadedErrorBar(1:nBPMs,bpmHDiffMean(i,:),bpmHDiffMean_err(i,:),{'Color',myCols(i,:)},0.9);
    plot(1:nBPMs,((-1)^i)*bpmHDiffMean(i,:),'LineWidth',2,'Color',myCols(i,:),'LineWidth',2);
    hold all;
end
set(gca,'XTickLabel',xTickLab)
set(gca,'XTick',1:8)
xlim([1 8])
ylim([-1.5 1.5])
plot([2.5 2.5],[-1.5 1.5],'k');
plot([6.5 6.5],[-1.5 1.5],'k');
% rotateXLabels(gca,45)
xlabel('BPM Index')
ylabel('Position [mm]')
title('Horizontal Position')
set(gcf, 'Colormap', myCols);
figColBar = colorbar;
allTicks = linspace(0,1,nDataSets+1);
set(figColBar,'YTick',allTicks(1:2:nDataSets));
set(figColBar,'YTickLabel',[-2.0 -1.5 -1.0 -0.5 0.0 0.5 1.0 1.5 2.0]);
ylabel(figColBar,'Input [V]');
format_plots;
if savePlots
    savePlot(saveDir,'HOrbitVsInput');
end;


figure
myCols = varycolor(nDataSets);
for i=1:nDataSets
    plot(1:nBPMs,bpmVDiffMean(i,:),'Color',myCols(i,:),'LineWidth',2);
%     plot(1:nBPMs,bpmHDiffMean(i,:),'LineWidth',2);%,'Color',myCols(i,:),'LineWidth',2);
    hold all;
end
set(gca,'XTickLabel',xTickLab)
set(gca,'XTick',1:8)
xlim([1 8])
ylim([-1.5 1.5])
plot([2.5 2.5],[-1.5 1.5],'k');
plot([6.5 6.5],[-1.5 1.5],'k');
% rotateXLabels(gca,45)
xlabel('BPM Index')
ylabel('Position [mm]')
title('Vertical Position')
set(gcf, 'Colormap', myCols);
figColBar = colorbar;
allTicks = linspace(0,2,nDataSets+1);
set(figColBar,'YTick',allTicks(1:2:nDataSets));
set(figColBar,'YTickLabel',[-2.0 -1.5 -1.0 -0.5 0.0 0.5 1.0 1.5 2.0]);
ylabel(figColBar,'Input [V]');
format_plots;
if savePlots; savePlot(saveDir,'VOrbitVsInput'); end;

figure;
% i=1;errorbar(1:nBPMs,-bpmHDiffMean(i,:),-bpmHDiffMean_err(i,:),'bo-','MarkerFaceColor','b')
% hold all;
% i=2;errorbar(1:nBPMs,bpmHDiffMean(i,:),bpmHDiffMean_err(i,:),'ro-','MarkerFaceColor','r')
i=1;errorbar(bpmPos-bpmPos(1),-bpmHDiffMean(i,:),-bpmHDiffMean_err(i,:),'bo-','MarkerFaceColor','b','LineWidth',2)
hold all;
i=2;errorbar(bpmPos-bpmPos(1),bpmHDiffMean(i,:),bpmHDiffMean_err(i,:),'ro-','MarkerFaceColor','r','LineWidth',2)
xlim([0 22.2]);
ylim([-1.55 1.55]);
xlabel('s [m]')
ylabel('Horizontal Offset [mm]')
legend('+2 V',' -2 V')
format_plots;
grid off;
