% load all interleaved data first, with all even datasets first half of
% dataSetNames, all odd second half of dataSetNames (using
% plotData_compareDatasets)

%% to copy in to plotData_compareDatasets
% dataSetNames = {
%%% ON
% '20161103_1300_TimChk_K1d_0_K2d_1_gain1000_Odd'
% '20161103_1313_TimChk_K1d_4_K2d_5_gain1000_Odd'
% '20161103_1255_TimChk_K1d_7_K2d_8_Even'
% '20161103_1316_TimChk_K1d_10_K2d_11_gain1000_Odd'
% '20161103_1306_TimChk_K1d_14_K2d_15_gain1000_Even'
% '20161103_1309_TimChk_K1d_21_K2d_22_gain1000_Odd'
% '20161103_1311_TimChk_K1d_28_K2d_29_gain1000_Even'
%%% OFF
% '20161103_1300_TimChk_K1d_0_K2d_1_gain1000_Even'
% '20161103_1313_TimChk_K1d_4_K2d_5_gain1000_Even'
% '20161103_1255_TimChk_K1d_7_K2d_8_Odd'
% '20161103_1316_TimChk_K1d_10_K2d_11_gain1000_Even'
% '20161103_1306_TimChk_K1d_14_K2d_15_gain1000_Odd'
% '20161103_1309_TimChk_K1d_21_K2d_22_gain1000_Even'
% '20161103_1311_TimChk_K1d_28_K2d_29_gain1000_Odd'
% };
%%

k1D = [0 4 7 10 14 21 28];

upPeakLoc = 379;

downPeakBeam = [506 505 506 504 506 506 506];
downPeakKick = [500 502 504 503 508 511 515];

saveDir = '/home/jack/PhaseFeedforward/Analysis/201611/Plots/20161103_TimChk';
saveFiles=0;

close all;
%%

myCols = varycolor(nDataSets/2);

figure(1)
figure(2)
figure(3)
for i=1:(nDataSets/2)
    diffPhaseAlong(i,:) = data{i}.meanPhaseAlongPulse(3,:)-data{i+(nDataSets/2)}.meanPhaseAlongPulse(3,:);
    if (diffPhaseAlong(i,430)<0)
        diffPhaseAlong(i,:) = -diffPhaseAlong(i,:);
    
        figure;
        plot(data{i+(nDataSets/2)}.meanPhaseAlongPulse(3,:)-data{i+(nDataSets/2)}.meanPhaseAlongPulse(3,505),'r','LineWidth',2);
        hold all
        plot(diffPhaseAlong(i,:)-diffPhaseAlong(i,505),'k','LineWidth',2);
        xlim([480 540])
        yLim = get(gca,'YLim');
        plot([downPeakBeam(i) downPeakBeam(i)],yLim,'r--','LineWidth',2)
        plot([downPeakKick(i) downPeakKick(i)],yLim,'k--','LineWidth',2)
        title(strrep(dataSetNames{i},'_','.'))
        ylim(yLim)
        xlabel('Sample No.')
        ylabel('Phase [degrees]')
        legend('PFF Off','PFF On - PFF Off')
        if saveFiles; savePlot(saveDir,sprintf('K1D_%d',k1D(i))); end;
        
     else        
        figure;
        plot(data{i}.meanPhaseAlongPulse(3,:)-data{i}.meanPhaseAlongPulse(3,505),'r','LineWidth',2);
        hold all
        plot(diffPhaseAlong(i,:)-diffPhaseAlong(i,505),'k','LineWidth',2);
        title(strrep(dataSetNames{i},'_','.'))
        xlim([480 540])
        yLim = get(gca,'YLim');
        plot([downPeakBeam(i) downPeakBeam(i)],yLim,'r--','LineWidth',2)
        plot([downPeakKick(i) downPeakKick(i)],yLim,'k--','LineWidth',2)
        xlabel('Sample No.')
        ylabel('Phase [degrees]')
        legend('PFF Off','PFF On - PFF Off')
        if saveFiles; savePlot(saveDir,sprintf('K1D_%d',k1D(i))); end;
                


    end


end

for i=1:(nDataSets/2)
     if (diffPhaseAlong(i,430)<0)
       figure(2)
        plot(data{i+(nDataSets/2)}.meanPhaseAlongPulse(2,:),'LineWidth',2,'Color',myCols(i,:))
        hold all;   
        figure(3)
        plot(data{i+(nDataSets/2)}.meanPhaseAlongPulse(3,:),'LineWidth',2,'Color',myCols(i,:))
        hold all;  

     else
         figure(2)
        plot(data{i}.meanPhaseAlongPulse(2,:),'LineWidth',2,'Color',myCols(i,:))
        hold all;   
        figure(3)
        plot(data{i}.meanPhaseAlongPulse(3,:),'LineWidth',2,'Color',myCols(i,:))
        hold all;
     end
     
             
    figure(1)
    plot(diffPhaseAlong(i,:),'LineWidth',2,'Color',myCols(i,:))
    hold all;
end

figure(1)
xlim([380 620])
ylim([-4 4])
xlabel('Sample No.')
ylabel('Phase [degrees]')
title('Downstream Phase: PFF On - PFF Off')
set(gcf, 'Colormap', myCols);
figColBar = colorbar;
allTicks = linspace(0,1,(nDataSets/2)+1);
set(figColBar,'YTick',allTicks(1:end-1));
set(figColBar,'YTickLabel',k1D);   
ylabel(figColBar,'K1 Delay');
if saveFiles; savePlot(saveDir,'Mon3Diff'); end;

figure(2)
xlim([201 539])
xlabel('Sample No.')
ylabel('Phase [degrees]')
title('Upstream Phase')
set(gcf, 'Colormap', myCols);
figColBar = colorbar;
allTicks = linspace(0,1,(nDataSets/2)+1);
set(figColBar,'YTick',allTicks(1:end-1));
set(figColBar,'YTickLabel',k1D);   
ylabel(figColBar,'K1 Delay');  
if saveFiles; savePlot(saveDir,'Mon2'); end;

figure(3)
xlim([349 655])
% ylim([-4 4])
xlabel('Sample No.')
ylabel('Phase [degrees]')
title('Downstream Phase')
set(gcf, 'Colormap', myCols);
figColBar = colorbar;
allTicks = linspace(0,1,(nDataSets/2)+1);
set(figColBar,'YTick',allTicks(1:end-1));
set(figColBar,'YTickLabel',k1D);   
ylabel(figColBar,'K1 Delay');
if saveFiles; savePlot(saveDir,'Mon3'); end;

figure(1)
xlim([480 540])
if saveFiles; savePlot(saveDir,'Mon3Diff_zoom'); end;

figure(2)
xlim([350 410])
if saveFiles; savePlot(saveDir,'Mon2_zoom'); end;

figure(3)
xlim([480 540])
if saveFiles; savePlot(saveDir,'Mon3_zoom'); end;


[coeffs,rsquare,confints]=nanpolyfit(k1D,downPeakKick-downPeakBeam,1);
optK1 = -coeffs(2)/coeffs(1)

figure;
plot(k1D,downPeakKick-downPeakBeam,'bo','MarkerFaceColor','b')
hold all;
plot(k1D,polyval(coeffs,k1D),'k','LineWidth',2)
yLim = get(gca,'YLim');
legend('Data','Fit','Location','NorthWest')
plot([optK1 optK1],[yLim(1) 0],'k--')
plot([0 optK1],[0 0],'k--')
title({'Kick Time Offset vs. K1 Delay' sprintf('Opt K1 Delay = %.1f clk cyc (%.1f ns)',optK1,optK1*2.8)})
xlabel('K1 Delay [clock cycles]')
ylabel('Kick-Beam Sample Offset')
if saveFiles; savePlot(saveDir,'OptDelayFit'); end;

