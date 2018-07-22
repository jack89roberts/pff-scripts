% load('/home/jack/PhaseFeedforward/Analysis/201603/20160316_1144_PosScanCT_H.mat');
% savePlotDir = '/home/jack/PhaseFeedforward/Analysis/201603/Plots/20160316_1144_PosScanCT_H';
% plane = 1;

load('/home/jack/PhaseFeedforward/Analysis/201603/20160316_1242_PosScanCT_V.mat');
savePlotDir = '/home/jack/PhaseFeedforward/Analysis/201603/Plots/20160316_1242_PosScanCT_V';
plane = 2;

mon1430Ratio = 0.29; %mm offset in mon1 due to kick at 360 per mm offset in CT.430
mon2430Ratio = 0.41; %mm offset in mon2 due to kick at 360 per mm offset in CT.430

nBins = 20;
%%
ct430H = meanBPMH{10};
ct430S = meanBPMS{10};
ct430V = meanBPMV{10};

ct430H = ct430H-nanmean(ct430H);
ct430V = ct430V-nanmean(ct430V);

mon1 = meanPulsePhase(1,:);
mon2 = meanPulsePhase(2,:);


binMon1 = NaN(1,nBins);
binMon2 = NaN(1,nBins);
binMon1_err = NaN(1,nBins);
binMon2_err = NaN(1,nBins);
binCT430H = NaN(1,nBins);
binCT430S = NaN(1,nBins);
binCT430V = NaN(1,nBins);
binCT430H_err = NaN(1,nBins);
binCT430S_err = NaN(1,nBins);
binCT430V_err = NaN(1,nBins);
binDiode1 = NaN(1,nBins);
binDiode2 = NaN(1,nBins);
binDiode1_err = NaN(1,nBins);
binDiode2_err = NaN(1,nBins);
binDiff12 = NaN(1,nBins);
binDiff12_err = NaN(1,nBins);
nPulsesBin = NaN(1,nBins);

if plane==1
    maxPos = max(ct430H);
    minPos = min(ct430H);
else
    maxPos = max(ct430V);
    minPos = min(ct430V);
end
binInt = (maxPos-minPos)/nBins;    
binStart = minPos+((0:(nBins-1))*binInt);
binEnd=minPos+((1:nBins)*binInt);

figure;
if plane==1
    plot(ct430H,mon2,'.');
else
    plot(ct430V,mon2,'.');
end
hold all;
yLim=get(gca,'YLim');
for i=1:nBins
    plot([binStart(i) binStart(i)],yLim,'k','LineWidth',2);
end
plot([binEnd(end) binEnd(end)],yLim,'k','LineWidth',2);
xlim([minPos maxPos])
title('DATA BINS')
if plane==1
    xlabel('CT430 H POSITION [mm]')
else
    xlabel('CT430 V POSITION [mm]')
end
ylabel('PHASE [degs]')

for i=1:nBins
    if (plane==1)
        binPulses = (ct430H>=binStart(i)) & (ct430H<binEnd(i));   
    else
        binPulses = (ct430V>=binStart(i)) & (ct430V<binEnd(i));           
    end
    nPulsesBin(i) = sum(binPulses);
    
    [binMon1(i),~,binMon1_err(i),~] = nanMeanStdErr(mon1(binPulses));
    [binMon2(i),~,binMon2_err(i),~] = nanMeanStdErr(mon2(binPulses));
    [binCT430H(i),~,binCT430H_err(i),~] = nanMeanStdErr(ct430H(binPulses));
    [binCT430V(i),~,binCT430V_err(i),~] = nanMeanStdErr(ct430V(binPulses));
    [binCT430S(i),~,binCT430S_err(i),~] = nanMeanStdErr(ct430S(binPulses));
    [binDiode1(i),~,binDiode1_err(i),~] = nanMeanStdErr(meanDiode(1,binPulses));
    [binDiode2(i),~,binDiode2_err(i),~] = nanMeanStdErr(meanDiode(2,binPulses));
    
    diff12 = mon2(binPulses)-mon1(binPulses);
    [binDiff12(i),~,binDiff12_err(i),~] = nanMeanStdErr(diff12);

end

binMon1PosH = binCT430H*mon1430Ratio;
binMon1PosH_err = binCT430H_err*mon1430Ratio;
binMon1PosV = binCT430V*mon1430Ratio;
binMon1PosV_err = binCT430V_err*mon1430Ratio;
binMon2PosH = binCT430H*mon2430Ratio;
binMon2PosH_err = binCT430H_err*mon2430Ratio;
binMon2PosV = binCT430V*mon2430Ratio;
binMon2PosV_err = binCT430V_err*mon2430Ratio;

if plane==1
    [fitDiff12,~,fitDiff12_err] = nanpolyfit(binCT430H,binDiff12,1,1./(binDiff12_err.^2));
else
    [fitDiff12,~,fitDiff12_err] = nanpolyfit(binCT430V,binDiff12,1,1./(binDiff12_err.^2));    
end
fitDiff12
fitDiff12_err = (fitDiff12-fitDiff12_err(1,:))./2

phasDep = fitDiff12./(mon2430Ratio-mon1430Ratio)
phaseDepErr = fitDiff12_err./(mon2430Ratio-mon1430Ratio)

if plane==1
    figure;
    errorbar(binCT430H,binDiff12,binDiff12_err,'bo','MarkerFaceColor','b');
    hold all;
    plot(binCT430H,polyval(fitDiff12,binCT430H),'k','LineWidth',2);
    legend('Data','Fit')
    herrorbar(binCT430H,binDiff12,binCT430H_err,'bo');
    xlabel('CT.BPM0430 H Position [mm]')
    ylabel('Phase Difference [degrees]')
    title('Mon2-Mon1 vs. Horizontal Position')
    format_plots
else
    figure;
    errorbar(binCT430V,binDiff12,binDiff12_err,'bo','MarkerFaceColor','b');
    hold all;
    plot(binCT430V,polyval(fitDiff12,binCT430V),'k','LineWidth',2);
    legend('Data','Fit')
    herrorbar(binCT430V,binDiff12,binCT430V_err,'bo');
    xlabel('CT.BPM0430 V Position [mm]')
    ylabel('Phase Difference [degrees]')
    title('Mon2-Mon1 vs. Vertical Position')
    format_plots    
end
%%

% figure;
% errorbar(binCT430S,binMon1,binMon1_err,'bo','MarkerFaceColor','b')
% hold all
% errorbar(binCT430S,binMon2,binMon2_err,'ro','MarkerFaceColor','r')
% legend('Mon1','Mon2')
% herrorbar(binCT430S,binMon1,binCT430S_err,'bo')
% herrorbar(binCT430S,binMon2,binCT430S_err,'ro')
% xlabel('CT 430 S POSITION [mm]')
% ylabel('Phase [degrees]')
% title('SUM')
% format_plots;
% 
if (plane==1)
    figure;
    errorbar(binCT430H,binMon1,binMon1_err,'bo','MarkerFaceColor','b')
    hold all
    errorbar(binCT430H,binMon2,binMon2_err,'ro','MarkerFaceColor','r')
    legend('Mon1','Mon2')
    herrorbar(binCT430H,binMon1,binCT430H_err,'bo')
    herrorbar(binCT430H,binMon2,binCT430H_err,'ro')
    xlabel('H POSITION [mm]')
    ylabel('Phase [degrees]')
    title('HORIZONTAL')
    format_plots;

    figure;
    errorbar(binCT430H,(binDiode1-nanmean(binDiode1))./nanmean(binDiode1),binDiode1_err./nanmean(binDiode1),'bo','MarkerFaceColor','b')
    hold all
    errorbar(binCT430H,(binDiode2-nanmean(binDiode2))./nanmean(binDiode2),binDiode2_err./nanmean(binDiode2),'ro','MarkerFaceColor','r')
    legend('Mon1','Mon2')
    herrorbar(binCT430H,(binDiode1-nanmean(binDiode1))./nanmean(binDiode1),binCT430H_err,'bo')
    herrorbar(binCT430H,(binDiode2-nanmean(binDiode2))./nanmean(binDiode2),binCT430H_err,'ro')
    xlabel('H POSITION [mm]')
    ylabel('Relative Offset')
    title('DIODE')
    format_plots;

    figure;
    errorbar(binCT430H,binCT430V,binCT430V_err,'bo','MarkerFaceColor','b')
    hold all
    herrorbar(binCT430H,binCT430V,binCT430H_err,'bo')
    xlabel('CT 430 H POSITION [mm]')
    ylabel('CT 430 V POSITION [mm]')
    title('VERTICAL VS. HORIZONTAL')
    format_plots;

    figure;
    errorbar(binCT430H,binCT430S,binCT430S_err,'bo','MarkerFaceColor','b')
    hold all
    herrorbar(binCT430H,binCT430S,binCT430H_err,'bo')
    xlabel('CT 430 H POSITION [mm]')
    ylabel('CT 430 S [A]')
    title('SUM VS. HORIZONTAL')
    format_plots;
else
    figure;
    errorbar(binCT430V,binMon1,binMon1_err,'bo','MarkerFaceColor','b')
    hold all
    errorbar(binCT430V,binMon2,binMon2_err,'ro','MarkerFaceColor','r')
    legend('Mon1','Mon2')
    herrorbar(binCT430V,binMon1,binCT430V_err,'bo')
    herrorbar(binCT430V,binMon2,binCT430V_err,'ro')
    xlabel('CT 430 V POSITION [mm]')
    ylabel('Phase [degrees]')
    title('VERTICAL')
    format_plots;
    
    figure;
    errorbar(binCT430V,(binDiode1-nanmean(binDiode1))./nanmean(binDiode1),binDiode1_err./nanmean(binDiode1),'bo','MarkerFaceColor','b')
    hold all
    errorbar(binCT430V,(binDiode2-nanmean(binDiode2))./nanmean(binDiode2),binDiode2_err./nanmean(binDiode2),'ro','MarkerFaceColor','r')
    legend('Mon1','Mon2')
    herrorbar(binCT430V,(binDiode1-nanmean(binDiode1))./nanmean(binDiode1),binCT430V_err,'bo')
    herrorbar(binCT430V,(binDiode2-nanmean(binDiode2))./nanmean(binDiode2),binCT430V_err,'ro')
    xlabel('V POSITION [mm]')
    ylabel('Relative Offset')
    title('DIODE')
    format_plots;
    
    figure;
    errorbar(binCT430V,binCT430H,binCT430H_err,'bo','MarkerFaceColor','b')
    hold all
    herrorbar(binCT430V,binCT430H,binCT430V_err,'bo')
    xlabel('CT 430 V POSITION [mm]')
    ylabel('CT 430 H POSITION [mm]')
    title('HORIZONTAL VS. VERTICAL')
    format_plots;

    figure;
    errorbar(binCT430V,binCT430S,binCT430S_err,'bo','MarkerFaceColor','b')
    hold all
    herrorbar(binCT430V,binCT430S,binCT430V_err,'bo')
    xlabel('CT 430 V POSITION [mm]')
    ylabel('CT 430 S [A]')
    title('SUM VS. VERTICAL')
    format_plots;
end