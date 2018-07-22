dataSetName= '20150715_1414_FF_Gain63_Gate185_350_Int_Odd';
dataDir = '/home/jack/PhaseFeedforward/Analysis/201507';
saveDir = '/home/jack/PhaseFeedforward/Analysis/201507/MatchingSourceOfStaticCROffset';

% NB - this script only works for CR.155 and CR.195 as it is written
bpmIndex1 = 22; % CR.155
bpmIndex2 = 23; % CR.195

meanSubtractSamples = 550:680;%(pulseSampleRange{3}(1)+5):(pulseSampleRange{3}(end)-5);

startPeaks = [558 575 592 626 650 668];
endPeaks = [566 582 597 635 654 674];
 
figSize = [66 374 1529 445];
twissSTickMarks = 0:2:70;

%% load data
addpath('../');
addpath('../Analysis/');
addpath('/home/jack/Documents/MATLAB/madx2Matlab');

savePlotDir = sprintf('%s/%s',saveDir,dataSetName);

load([dataDir '/' dataSetName '.mat']);

%%
nPeaks = length(startPeaks);

bpmH1 = meanBPMHAlongPulse{bpmIndex1};
bpmH2 = meanBPMHAlongPulse{bpmIndex2};

bpmH1 = bpmH1-nanmean(bpmH1(meanSubtractSamples));
bpmH2 = bpmH2-nanmean(bpmH2(meanSubtractSamples));

figure;
set(gcf,'PaperPositionMode','auto');
plot(bpmH1,'b');
hold all;
plot(bpmH2,'r');
xlim([meanSubtractSamples(1) meanSubtractSamples(end)])
legend(bpmNames{bpmIndex1},bpmNames{bpmIndex2})
title([strrep(dataSetName,'_',' ') ' (amplitude slopes with thick lines used for matching)'])
xlabel('Sample No.');
ylabel('Horizontal Position [mm]');
grid
set(gcf,'Position',figSize);
yLim = get(gca,'YLim');
for i=1:nPeaks
    plot(startPeaks(i):endPeaks(i),bpmH1(startPeaks(i):endPeaks(i)),'b','LineWidth',4);
    plot(startPeaks(i):endPeaks(i),bpmH2(startPeaks(i):endPeaks(i)),'r','LineWidth',4);
%     plot(endPeaks(i),bpmH1(endPeaks(i)),'vb');
%     plot(startPeaks(i),bpmH2(startPeaks(i)),'^r');
%     plot(endPeaks(i),bpmH2(endPeaks(i)),'vr');
%     plot([startPeaks(i) startPeaks(i)],yLim,'k');
%     plot([endPeaks(i) endPeaks(i)],yLim,'k')
end
savePlot(savePlotDir,'slopesUsed');

peakAmplitudesBPM1 = NaN(1,nPeaks);
peakAmplitudesBPM2 = NaN(1,nPeaks);
for i=1:nPeaks    
    peakAmplitudesBPM1(i) = bpmH1(endPeaks(i))-bpmH1(startPeaks(i));
    peakAmplitudesBPM2(i) = bpmH2(endPeaks(i))-bpmH2(startPeaks(i));
    
    currentsFile = fopen('matchCROffset.currents','w');
    fprintf(currentsFile,'matchCR155 = %.8f;\n',peakAmplitudesBPM1(i)./1000);
    fprintf(currentsFile,'matchCR195 = %.8f;\n',peakAmplitudesBPM2(i)./1000);
    fclose(currentsFile);
    
    system('madx < matchCROffset.madx');
    twissDat = madx2matlab.parseTFSTableFromFile('matchedCROffset.twiss');
    
    twissName = [twissDat.DATA.NAME];
    twissS = [twissDat.DATA.S];
    twissX = [twissDat.DATA.X];
    twissPX = [twissDat.DATA.PX];
    
    if (i==1)
        devNames = twissName;
        devPositions = twissS;
        twissIndexBPM1 = find(strcmp('"CR.BPM0155"',devNames));
        twissIndexBPM2 = find(strcmp('"CR.BPM0195"',devNames));
        matchedPeakX = NaN(nPeaks,length(twissX));
        matchedPeakPX = NaN(nPeaks,length(twissPX));
    end
    
    matchedPeakX(i,:) = twissX;
    matchedPeakPX(i,:) = twissPX;
end

xFig = figure;
set(gcf,'PaperPositionMode','auto');
pxFig = figure;
set(gcf,'PaperPositionMode','auto');
legendNames = cell(1,nPeaks);
for i=1:nPeaks
    legendNames{i} = sprintf('Peak %d',startPeaks(i));

    figure(xFig);
    plot(twissS,matchedPeakX(i,:));
    hold all;
    
    figure(pxFig);
    plot(twissS,matchedPeakPX(i,:));
    hold all;   
end

figure(xFig);
set(gcf,'Position',figSize);
grid;
title('MATCHED HORIZONTAL ORBIT: X');
legend(legendNames);
set(gca,'XTick',devPositions(1:5:end));
set(gca,'XTickLabel',{devNames{1:5:end}});
rotateXLabels(gca,90);
ylabel('Horizontal Position [m]');
yLim = get(gca,'YLim');
plot([devPositions(twissIndexBPM1) devPositions(twissIndexBPM1)],yLim, 'k');
plot([devPositions(twissIndexBPM2) devPositions(twissIndexBPM2)],yLim, 'k');
savePlot(savePlotDir,'matchedX');

figure(pxFig);
set(gcf,'Position',figSize);
title('MATCHED HORIZONTAL ORBIT: PX');
grid;
legend(legendNames);
set(gca,'XTick',devPositions(1:5:end));
set(gca,'XTickLabel',{devNames{1:5:end}});
rotateXLabels(gca,90);
ylabel('PX [m]');
yLim = get(gca,'YLim');
plot([devPositions(twissIndexBPM1) devPositions(twissIndexBPM1)],yLim, 'k');
plot([devPositions(twissIndexBPM2) devPositions(twissIndexBPM2)],yLim, 'k');
savePlot(savePlotDir,'matchedPX');

figure;
set(gcf,'PaperPositionMode','auto');
plot(startPeaks,peakAmplitudesBPM1,'b');
hold all;
plot(startPeaks,peakAmplitudesBPM2,'r');
plot(startPeaks,matchedPeakX(:,twissIndexBPM1)*1000,'bx');
plot(startPeaks,matchedPeakX(:,twissIndexBPM2)*1000,'rx');
legend('Data Peaks BPM1','Data Peaks BPM2','Matched Peaks BPM1','Matched Peaks BPM2');
xlabel('Sample No. Peak Start');
ylabel('Peak amplitude [mm]');
grid;
savePlot(savePlotDir,'matchResults');