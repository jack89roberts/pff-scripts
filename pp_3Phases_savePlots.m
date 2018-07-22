plotsDesc = 'R56Scan';
baseDir = '/home/jack/PhaseFeedforward/Analysis/';

normSamp = 50;%98; %90
plotStartSamp = 40;
plotEndSamp = plotStartSamp + 85;
%%
close all;

baseDate = regexprep(DataSetBaseTitle,'[^\w'']','_');
saveDir = [baseDir baseDate '_' plotsDesc '/'];

if (~exist(saveDir,'dir'))
    mkdir(saveDir);
end

dataSetsFile = fopen([saveDir 'dataSetList.txt'],'w');
for i=1:nDataSets
    fprintf(dataSetsFile,'%s, %s\n', dataSetNames{i}, dataSets{i});
end
fclose(dataSetsFile);
%% Mean pulse phase
figure;
set(gcf,'PaperPositionMode','auto'); % forces figures to be saved with same aspect ratio

for i=1:nDataSets
  tmp3 = nanmean(MIX3a{i});
  tmp3 = tmp3-tmp3(normSamp);
  timeAxis = 4*(1:length(tmp3))-4;
  plot(timeAxis,tmp3,'-','Color',kolors{i},'LineWidth',2);
  hold all;
end
legend(rawNames,'Location','SouthEast');
for i=1:nDataSets
    tmp1 = nanmean(MIX1a{i});
    tmp1 = tmp1-tmp1(normSamp);
    timeAxis = 4*(1:length(tmp1))-4;
    plot(timeAxis,tmp1,'--','Color',kolors{i},'LineWidth',2);
end
xlim([timeAxis(plotStartSamp) timeAxis(plotEndSamp)])
set(gcf,'units','normalized','position',[0 0 0.55 0.5])
xlabel('Time [ns]');
ylabel('Mean Phase [degrees]');
title({dataDescription, DataSetBaseTitle});
grid on;
annotation(gcf,'textbox',...
    [0.143863636363634 0.14070139168596 0.219682167087077 0.116076386091817],...
    'String',{'Solid Lines: CLEX','Dashed Lines: CT'},...
    'FontSize',14,...
    'FontName','Helvetica',...
    'FitBoxToText','off',...
    'BackgroundColor',[1 1 1]);
set(gca,'FontSize',14)
set(findall(gcf,'type','text'),'FontSize',14)
yAxLim = get(gca,'YLim');
plot([timeAxis(min(sampleRange)) timeAxis(min(sampleRange))],yAxLim,'k','LineWidth',2);
plot([timeAxis(max(sampleRange)) timeAxis(max(sampleRange))],yAxLim,'k','LineWidth',2);
saveStr = sprintf('%sMeanPulsePhase',saveDir);
print([saveStr '.png'],'-dpng');
savefig([saveStr '.fig']);
clf;
%% Phase jitter
stdwinsize = 0; %1
for i=1:nDataSets
   winstdMIX3a{i} = zeros(201,1);
   winstdMIX1a{i} = zeros(201,1);
   for j=40:250-stdwinsize-1
     winstdMIX3a{i}(j,1) = std2(MIX3a{i}(:,j:j+stdwinsize));
     winstdMIX1a{i}(j,1) = std2(MIX1a{i}(:,j:j+stdwinsize));
   end
  timeAxis = 4*(1:length(winstdMIX3a{i}))-4;
  plot(timeAxis,winstdMIX3a{i},'-x','Color', kolors{i},'LineWidth',2);
  hold all;
end
annotation(gcf,'textbox',...
    [0.145 0.747368058352627 0.219682167087078 0.116076386091817],...
    'String',{'Solid Lines: CLEX','Dashed Lines: CT'},...
    'FontSize',14,...
    'FontName','Helvetica',...
    'BackgroundColor',[1 1 1]);
legend(rawNames);
for i=1:nDataSets
    plot(timeAxis,winstdMIX1a{i},'--','Color', kolors{i},'LineWidth',2);
end
xlim([timeAxis(plotStartSamp) timeAxis(plotEndSamp)])
title({dataDescription, DataSetBaseTitle});
xlabel('Time [ns]')
ylabel('Phase Jitter [degrees]');
grid on;
set(gcf,'units','normalized','position',[0 0 0.55 0.5])
set(gca,'FontSize',16)
set(findall(gcf,'type','text'),'FontSize',16)
yAxLim = get(gca,'YLim');
plot([timeAxis(min(sampleRange)) timeAxis(min(sampleRange))],yAxLim,'k','LineWidth',2);
plot([timeAxis(max(sampleRange)) timeAxis(max(sampleRange))],yAxLim,'k','LineWidth',2);

saveStr = sprintf('%sPhaseJitter',saveDir);
print([saveStr '.png'],'-dpng');
savefig([saveStr '.fig']);
clf;

%% Mean phase vs time
set(gcf,'units','normalized','position',[0 0 0.5 0.5])

for i=1:nDataSets
    tmp3 = meanMIX3a{i} - nanmean(meanMIX3a{i});
    std3 = nanstd(tmp3);
    plot(tmp3,'-r','LineWidth',2);
    hold all;
    tmp1 = meanMIX1a{i} - nanmean(meanMIX1a{i});
    std1 = nanstd(tmp1);
    plot(tmp1,'-b','LineWidth',2);

    legend(sprintf('CLEX (std=%.2f^o)',std3),sprintf('CT (std=%.2f^o)',std1));
    title({sprintf('%s (%s)',dataDescription, DataSetBaseTitle),rawNames{i}});
    grid on;
    %ylim([-10 10]);
    %xlim([0 250]);

    xlabel('Time [Pulse No.]');
    ylabel('Mean Phase [degrees]');

    set(gca,'FontSize',16)
    set(findall(gcf,'type','text'),'FontSize',16)

    saveStr = sprintf('%sMeanPhaseVsTime_%s',saveDir,strrep(strrep(rawNames{i},',','_'),' ',''));
    print([saveStr '.png'],'-dpng');
    savefig([saveStr '.fig']);
    hold off;
end

clf;

%%
min3 = 1000;
min1 = 1000;
max3 = -1000;
max1 = -1000;
for i=1:nDataSets
    tmpMin3 = min(meanMIX3a{i});
    tmpMax3 = max(meanMIX3a{i});
    tmpMin1 = min(meanMIX1a{i});
    tmpMax1 = max(meanMIX1a{i}); 
    if(tmpMin3 < min3)
        min3 = tmpMin3;
    end
    if(tmpMax3 > max3)
        max3 = tmpMax3;
    end
    if(tmpMin1 < min1)
        min1 = tmpMin1;
    end
    if(tmpMax1 > max1)
        max1 = tmpMax1;
    end
end
min3 = floor(min3);
min1 = floor(min1);
max3 = ceil(max3);
max1 = ceil(max1);

for i=1:nDataSets
    scatter(meanMIX1a{i}, meanMIX3a{i},'o');
  
    [cf,pval] = corr(meanMIX1a{i}, meanMIX3a{i});
  
    title({sprintf('%s (%s)',dataDescription, DataSetBaseTitle),sprintf('%s (corr = %.2f)',rawNames{i},cf)});
    xlabel('Phase CT [12 GHz deg.]')
    ylabel('Phase CLEX [12 GHz deg.]')
    xlim([min1 max1]);
    ylim([min3 max3]);

    grid on;
    set(gca,'FontSize',14)
    set(findall(gcf,'type','text'),'FontSize',14)
    
    saveStr = sprintf('%sCorrelation_%s',saveDir,strrep(strrep(rawNames{i},',','_'),' ',''));
    print([saveStr '.png'],'-dpng');
    savefig([saveStr '.fig']);
    hold off;

end
