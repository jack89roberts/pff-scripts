allEven = data(1:17);
allOdd = data(18:34);
kick = -4000:500:4000;

meanEven = NaN(1,17);
meanOdd = NaN(1,17);
meanEven_err = NaN(1,17);
meanOdd_err = NaN(1,17);
for ds=1:17
    meanEven(ds) = nanmean(allEven{ds}.meanPulsePhase(3,:));
    meanOdd(ds) = nanmean(allOdd{ds}.meanPulsePhase(3,:));
    meanEven_err(ds) = nanmean(allEven{ds}.meanPulsePhase_err(3,:));
    meanOdd_err(ds) = nanmean(allOdd{ds}.meanPulsePhase_err(3,:));

end
diffKick = abs(meanOdd-meanEven);
diffKick(1:8) = -diffKick(1:8);
diffKick_err = sqrt(meanEven_err.^2 + meanOdd_err.^2);
fitDiffKick = nanpolyfit(kick,diffKick,1);

figure;
errorbar(kick,meanEven,meanEven_err,'bo');
hold all;
errorbar(kick,meanOdd,meanOdd_err,'ro');

figure;
errorbar(kick,diffKick,diffKick_err,'bo','MarkerSize',markerSizeBig,'MarkerFaceColor','b');
xlim([-4500 4500])
xlabel('DAC Output [counts]')
ylabel('Phase Offset [degrees]')
hold all;
xVals = -4500:500:4500;
plot(xVals,polyval(fitDiffKick,xVals),'b')
title({dataDescr 'Phase Offset vs. DAC Output'})
legend('Data',sprintf('Fit (%.1f^o per 1000 cnts',fitDiffKick(1)*1000),'Location','NorthWest')
format_plots;
hold off;
savePlot(saveDir,'KickVsDAC')

figure;
diffToFit = diffKick-polyval(fitDiffKick,kick);
plot(kick,diffToFit);