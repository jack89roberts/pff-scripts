sampleRange = 45:125;%111:130;%40:130;%50:129;%111:140; % 140:165 % JACK CHANGED RANGE

% for R56 scan 20141205
if (exist('backupMergedData','var')) 
    mergedData = backupMergedData;
end
backupMergedData = mergedData;
mergedData = {mergedData{1},mergedData{2},mergedData{3},mergedData{4},[mergedData{7} mergedData{8}]};
dataSetNames = {
    'R56 = +0.5',...
    'R56 = +0.4',...
    'R56 = +0.3',...
    'R56 = +0.2',...
    'R56 = -0.2',...
};
rawNames = dataSetNames;
nDataSets = length(mergedData);


bprCCCal = 1;
close all;
   kolors={'k','b','r',[0,0.7,0],[1 0.8 0],[0.8 0.2 0.8],[0.3 0.9 0.9],'y'};
   calMIX1 = 1;
   calMIX2 = 5;
   calMIX3 = 1;
  
   MIX1 = cell(1,nDataSets);
   mix2 = cell(1,nDataSets);
   MIX3 = cell(1,nDataSets);
   
   MIX1a = cell(1,nDataSets);
   MIX3a = cell(1,nDataSets);
   meanMIX3a = cell(1,nDataSets);
   winmeanMIX3a = cell(1,nDataSets);
   winstdMIX3a = cell(1,nDataSets);
   
   meanbpmTL1_HDriftSubMIX3a= cell(1,nDataSets);
   bpmTL1_H = cell(1,nDataSets);
   bpmCT_H = cell(1,nDataSets);
   bpmTL1_Ha = cell(1,nDataSets);
   bpmCT_Ha = cell(1,nDataSets);
   
   bprCC_raw = cell(1,nDataSets);
   bprCR_raw = cell(1,nDataSets);
   bprCC_phase = cell(1,nDataSets);
   bprCR_phase = cell(1,nDataSets);
   bprCCa = cell(1,nDataSets);
   
for i=1:nDataSets
    fprintf(1,'DataSet %d\n',i);
% %     allProcessedPhaseData{i} = extractProcessedPhaseData(mergedData{i},calibrationConstants);
% %     
% %     allPhases{i} = allProcessedPhaseData{i}.phases;
% %     allDiodes{i} = allProcessedPhaseData{i}.diodes;
% %     20141216_1845_Gate325to450_Gain-30_WiggleAmplitudeAlong3
% %     
% %     allMeanPhases{i} = allProcessedPhaseData{i}.meanPhases;
% %     allMeanPulsePhases{i} = allProcessedPhaseData{i}.meanPulsePhases;
% %     allStdPhases{i} = allProcessedPhaseData{i}.stdPhases;
% %     
% %     allMeanMixers{i} = allProcessedPhaseData{i}.meanMixers;
% %     allMeanPulseMixers{i} = allProcessedPhaseData{i}.meanPulseMixers;
    
    %pp_CleanAndSynchr(mergedData);
    bprCC_raw{i} =  extractCTFSignalFromMergedData('CC.STBPR0915S.Samples.samples.value',mergedData{i});
 %   bprCR_raw{i} =  extractCTFSignalFromMergedData('CR.STBPR0505S.Samples.samples.value',mergedData{i});
   
    bprCC_phase{i} = asind(bprCC_raw{i}./bprCCCal);
 %   bprCR_phase{i} = asind(bprCR_raw{i}./bprCRCal);
    
%    bpmCT_H{i} =  extractCTFSignalFromMergedData('CT.SVBPM0285H.Samples.samples.value',mergedData{i});
%    bpmTL1_H{i} = extractCTFSignalFromMergedData('CT.SVBPI0608H.Samples.samples.value',mergedData{i});
    
    %MIX1{i} = squeeze(mixers(1,:,:));
    
    
    if (useDiode)
      [mixers,diodes] = extractMixerDiode(mergedData{i},0);
      MIX1{i} = getPhaseMixerDiode(squeeze(mixers(2,:,:)),squeeze(diodes(2,:,:)),calibrationConstants(2,1),calibrationConstants(2,4));
      MIX3{i} = getPhaseMixerDiode(squeeze(mixers(3,:,:)),squeeze(diodes(3,:,:)),calibrationConstants(3,1),calibrationConstants(3,4));
    elseif(~useDiode)
      [mixers,diodes] = extractMixerDiode(mergedData{i},0);
      MIX1{i} = getPhaseMixerDiode(squeeze(mixers(2,:,:)),[],calibrationConstants(2,1),calibrationConstants(2,4));
      MIX3{i} = getPhaseMixerDiode(squeeze(mixers(3,:,:)),[],calibrationConstants(3,1),calibrationConstants(3,4));
    else
       [mixers,diodes] = extractMixerDiode(mergedData{i},1);
       MIX1{i} = squeeze(mixers(2,:,:));
       MIX3{i} = squeeze(mixers(3,:,:)); 
    end
     
%     MIX1a{i} = getAligned2(MIX1{i},10,300);
%     MIX3a{i} = getAligned2(MIX3{i},25,300); 
    
    bprCCa{i}  = getAligned2(bprCC_raw{i},40,300);  
    %bpmTL1_Ha{i} = getAligned2(bpmTL1_H{i},20,200);  
    bpmTL1_Ha{i} = bpmTL1_H{i};
    %bpmCT_Ha{i} = getAlignedBPM(bpmCT_S{i},bpmCT_S{i}, bpmCT_S{i}, 20,200);   
    
    % Jack added bad pulse removal based on diode signal level
    removeThreshold = 0.02;
    
    % first remove pulses which never go above a few millivolts, which can
    % cause alignment to fail.
    isBadDIO1 = ~sum(abs(squeeze(diodes(2,:,:)))>removeThreshold,2);
    isBadDIO3 = ~sum(abs(squeeze(diodes(3,:,:)))>removeThreshold,2);
    isBadPulse = isBadDIO1 | isBadDIO3;
    tmpDIO1 = squeeze(diodes(2,~isBadPulse,:));
    tmpMIX1 = MIX1{i}(~isBadPulse,:);
    tmpDIO3 = squeeze(diodes(3,~isBadPulse,:));
    tmpMIX3 = MIX3{i}(~isBadPulse,:);
    
    [DIO1a,b] = getAligned2(tmpDIO1,20,300,0,{tmpMIX1}); % Jack - align signals using diode rather than phase directly, which can have large noise spikes in baseline
    MIX1a{i} = b{1};
    meanDIO1a = nanmean(DIO1a(:,sampleRange),2);
    isGoodDIO1a = meanDIO1a-min(meanDIO1a) < removeThreshold;
    
    [DIO3a,b] = getAligned2(tmpDIO3,20,300,0,{tmpMIX3});
    if (i==3)
    
    end
    MIX3a{i} = b{1};
    meanDIO3a = nanmean(DIO3a(:,sampleRange),2);
    isGoodDIO3a = meanDIO3a-min(meanDIO3a) < removeThreshold;
    
    isGoodPulse = isGoodDIO1a & isGoodDIO3a;
    MIX1a{i} = MIX1a{i}(isGoodPulse,:);
    MIX3a{i} = MIX3a{i}(isGoodPulse,:);
    %bprCCa{i} = bprCCa{i}(isGoodPulse,:);
    %bpmTL1_Ha{i} = bpmTL1_H{i}(isGoodPulse,:);

end

if(0)
% figure(94);
% pp_DataSets(bprCT_raw);
% title([DataSetBaseTitle 'CT-BPR RAW SIGNAL']);
% 
% figure(95);
% pp_DataSets(bprCR_raw);
% title([DataSetBaseTitle 'CR-BPR RAW SIGNAL']);

figure(96);
pp_DataSets(bprCC_raw);
title([DataSetBaseTitle 'CC-BPR RAW PHASE']);


figure(97);
pp_DataSets(bprCCa);
title([DataSetBaseTitle 'CC-BPR aligned PHASE']);

figure(98);
pp_DataSets(bpmCT_H);
title([DataSetBaseTitle ' RAW DATA' ]);

figure(99);
pp_DataSets(MIX1);
title([DataSetBaseTitle ' RAW DATA' ]);

figure(990);
pp_DataSets(MIX1a);
title([DataSetBaseTitle ' RAW DATA' ]);

figure(100);
pp_DataSets(MIX3);
title([DataSetBaseTitle ' RAW DATA' ]);

figure(101);
pp_DataSets(MIX3a);
title([DataSetBaseTitle ' MIX CLEX All traces cleaned and aligned' ]);

figure(102);
pp_DataSets(MIX1a);
title([DataSetBaseTitle ' MIX CT All traces cleaned and aligned' ]);

return;
end

%%

figure(104)

for i=1:nDataSets
   meantrace =  nanmean(MIX3a{i});
  plot(meantrace - meantrace(1),'LineWidth',2)
  hold all;
end

title([DataSetBaseTitle ' Mean traces (over all pulses) normalized to 1st point' ]);

legend(dataSetNames);

%%
figure(105)

for i=1:nDataSets
  plot(nanmean(MIX3a{i}),'LineWidth',2)
  %shadedErrorBar(1:length(MIX3a{i}),nanmean(MIX3a{i}),nanstd(MIX3a{i}),{'Color',[rand() rand() rand()]},1)
  hold all;
end

for i=1:nDataSets
  plot(nanmean(MIX1a{i}),'-k','LineWidth',1)
  hold all;
end
title([DataSetBaseTitle ' Mean traces (over all pulses)' ]);

legend(dataSetNames);

%%
% figure(200)
% for i=1:nDataSets
%   meanbprCT{i} = nanmean(bprCT_phase{i}(:,215:240),2);
%   meanbprCR{i} = nanmean(bprCR_phase{i}(:,260:300),2);
%   
%   
%   
%   plot(meanbprCT{i},'LineWidth',2);
%   hold all;
%   plot(meanMIX1a{i},'LineWidth',2);
%   plot(meanMIX3a{i},'LineWidth',2);
% end
% title([DataSetBaseTitle ' Evolution of mean CT.BPR' ]);
% legend('CT BPR', 'CT Frascati', 'CLEX Frascati');

%%



figure(106)

legenda = dataSetNames;
for i=1:nDataSets
  meanMIX3a{i} = nanmean(MIX3a{i}(:,sampleRange),2); 
  meanMIX1a{i} = nanmean(MIX1a{i}(:,sampleRange),2);

  meanBprCC{i} = nanmean(bprCCa{i}(:,90:100),2);

  
  %meanMIX3a{i} = mean(MIX3a{i}(:,150:180),2);
  %meanMIX1a{i} = mean(MIX1a{i}(:,150:180),2);

  plot(meanMIX3a{i},'LineWidth',2);
  legenda{i} = [dataSetNames{i} '  stddev = ' num2str(std(meanMIX3a{i}))];
  hold all;
end
legend(legenda);
title([DataSetBaseTitle ' Evolution of mean pulse (over samples)' ]);

avgwinsize = 5;
for i=1:nDataSets
  
  winmeanMIX3a{i} = meanMIX3a{i};
  
  winmeanMIX3a{i}(1:avgwinsize) = mean(meanMIX3a{i}(1:avgwinsize));
  
  [n,m] = size(winmeanMIX3a{i});
  
  for j=avgwinsize+1:n
      winmeanMIX3a{i}(j) = mean(meanMIX3a{i}(j-avgwinsize:j));
  end
  
  plot(winmeanMIX3a{i},'--','LineWidth',1);
  
end

%%
figure(107)
legenda = dataSetNames;
stdwinsize = 0; %1
for i=1:nDataSets
   winstdMIX3a{i} = zeros(201,1);
   winstdMIX1a{i} = zeros(201,1);
   for j=40:250-stdwinsize-1
     winstdMIX3a{i}(j,1) = std2(MIX3a{i}(:,j:j+stdwinsize));
     winstdMIX1a{i}(j,1) = std2(MIX1a{i}(:,j:j+stdwinsize));
   end
  
  plot(winstdMIX3a{i},'-x','Color', kolors{i},'LineWidth',2);
  legenda{2*i-1} = [rawNames{i} ' CLEX Phase'];%[dataSetNames{i} ' CLEX Phase'];
  hold all;
  plot(winstdMIX1a{i},'--','Color', kolors{i},'LineWidth',2);
  legenda{2*i} = [rawNames{i} ' CT Phase'];%[dataSetNames{i} ' CT Phase']; 
end
legend(legenda);
title([DataSetBaseTitle ' Sliding std dev' ]);
%%
figure(108)

for i=1:nDataSets
  meanDriftSubMIX3a{i} = meanMIX3a{i} - winmeanMIX3a{i};
  plot(meanDriftSubMIX3a{i},'LineWidth',2);
  legenda{i} = [dataSetNames{i} '  stddev = ' num2str(std(meanDriftSubMIX3a{i}))];
  hold all;
end
title('Means over pulses drift subtracted');
legend(legenda);

%% Rejects outlayers for correlation computation 

figure(109);
subplot(1,nDataSets,1);
for i=1:nDataSets
  subplot(1,nDataSets,i);
  
  mi1 = nanmean(meanMIX1a{i});
  st1 = nanstd(meanMIX1a{i});
  d_m1 = mi1 - st1*3;
  d_p1 = mi1 + st1*3;
  mi3 = nanmean(meanMIX3a{i});
  st3 = nanstd(meanMIX3a{i});
  d_m3 = mi3 - st3*3;
  d_p3 = mi3 + st3*3;
  
  arr1 = zeros(1,1);
  arr3 = zeros(1,1);
  arrCTH = zeros(1,1);
  arrCCR = zeros(1,1);
  
  narr = 0;
  for j=1:min(length(meanMIX3a{i}), length(meanMIX1a{i}) )
      if (meanMIX3a{i}(j,1) < d_m3)
         continue; 
      end
      if (meanMIX3a{i}(j,1) > d_p3) 
         continue;
      end
      if (meanMIX1a{i}(j,1) < d_m1)
         %display('skipped  ');
         continue; 
      end
      if (meanMIX1a{i}(j,1) > d_p1) 
         continue; 
      end
      narr = narr + 1;
      arr3(narr,1) = meanMIX3a{i}(j,1);
      arr1(narr,1) = meanMIX1a{i}(j,1);
%      arrCTH(narr,1) = meanbpmCT_H{i}(j,1);
% can be shorter array
%      arrCCR(narr,1) = meanBprCC{i}(j,1);

     % if (i == 2)
     %     fprintf(1,'%d %f %f %f \n',j,meanMIX1a{i}(j,1), d_m1, d_p1);
     % end
  end
  
  meanMIX1a{i} = arr1;
  meanMIX3a{i} = arr3;
  meanbpmCT_H{i} = arrCTH;
%  meanBprCC{i} = arrCCR;
  
  %scatterColourProgression(meanMIX1a{i}, meanMIX3a{i},'.');
  scatter(meanMIX1a{i}, meanMIX3a{i},'o');
  
  [cf,pval] = corr(meanMIX3a{i}, meanMIX1a{i});
    
  legenda{i} = [dataSetNames{i} '  cf = ' num2str(cf)];
  
  title({ dataSetNames{i}, ['  cf = ' num2str(cf)] });
  xlabel('Phase CT [12 GHz deg.]')
  ylabel('Phase CLEX [12 GHz deg.]')
  %xlim([-10 10]);
  %ylim([2 10]);
end
%%


%%
figure(110);
subplot(1,nDataSets,1);
for i=1:nDataSets
  subplot(1,nDataSets,i);
  plot(meanMIX3a{i} - nanmean(meanMIX3a{i}),'-r');
  hold all;

%  plot(meanMIX3a{i}(1:2:end) - nanmean(meanMIX3a{i}));
%  hold on;20141216_1845_Gate325to450_Gain-30_WiggleAmplitudeAlong3
%  plot(meanMIX3a{i}(2:2:end) - nanmean(meanMIX3a{i}));
  plot(meanMIX1a{i} - nanmean(meanMIX1a{i}),'-b');
  plot(meanbpmCT_H{i} - nanmean(meanbpmCT_H{i}),'-k');
  legend('CLEX','CT','CT H pos');
  title(dataSetNames{i});
end


return; 
%% Continue if you have CT.BPM Hor




figure(111);
subplot(1,nDataSets,1);
for i=1:nDataSets
  subplot(1,nDataSets,i);
  scatterColourProgression(meanMIX3a{i}, meanbpmCT_H{i},'.');
  [cf,pval] = corr(meanMIX3a{i}, meanbpmCT_H{i});
  title({'phase in CLEX vs CT.BPM285H' ,dataSetNames{i}, ['  cf = ' num2str(cf)] });
end;


figure(112);
subplot(1,nDataSets,1);
for i=1:nDataSets
  subplot(1,nDataSets,i);
  scatterColourProgression(meanMIX1a{i}, meanbpmCT_H{i},'.');
  [cf,pval] = corr(meanMIX1a{i}, meanbpmCT_H{i});
  title({'phase in CT vs CT.BPM285H ' ,dataSetNames{i}, ['  cf = ' num2str(cf)] });
  
end;

%%

figure(201);
subplot(1,nDataSets,1);
for i=1:nDataSets
  subplot(1,nDataSets,i);
  scatterColourProgression(meanMIX1a{i}, meanBprCC{i},'.');
  [cf,pval] = corr(meanMIX1a{i}, meanbprCC{i});
  title({'phase in CT: frascati vs bpr' ,dataSetNames{i}, ['  cf = ' num2str(cf)] });
  
end;

%%
figure(1060)
for i=1:nDataSets
  meanbpmCT_H{i} = mean(bpmCT_H{i}(:,150:180),2);
  plot(meanbpmCT_H{i},'LineWidth',2);
  hold all;
end
title([DataSetBaseTitle ' Evolution of mean CT.BPM285H' ]);

return;














figure(106)

plot(nanmean(MIX3a{2}) - nanmean(MIX3a{1}),'LineWidth',2)
hold all;
plot(nanmean(MIX3a{3}) - nanmean(MIX3a{1}),'LineWidth',2)

plot(nanmean(MIX3a{3}) - nanmean(MIX3a{2}),'LineWidth',2)


figure(107)

plot(nanmean(MIX3a{5}) - nanmean(MIX3a{4}),'LineWidth',2)
hold all;
plot(nanmean(MIX3a{6}) - nanmean(MIX3a{4}),'LineWidth',2)

plot(nanmean(MIX3a{6}) - nanmean(MIX3a{5}),'LineWidth',2)


%%
%gain = [0,20,40,63,-20,-40,-63];
gain = [0.5,0.4,0.3,0.2,-0.2,-0.2];
subtractRange = 60:110;

figure;
subplot(1,2,1)
for i=1:nDataSets
    plot(gain(i),nanmean(winstdMIX3a{i}(sampleRange)),'bo'); %-
    hold all;
    %plot(gain(i),nanmean(winstdMIX3a{i}(subtractRange)),'ro');
end
title('MEAN of STDs');
%ylim([0.5 1.3])
%xlim([-70 70])
hold off;
subplot(1,2,2)
for i=1:nDataSets
    [cf,pval] = corr(meanMIX3a{i}, meanMIX1a{i});
    plot(gain(i),cf,'bo');
    hold all;
    %[cf,pval] = corr(nanmean(MIX3a{i}(:,subtractRange),2), nanmean(MIX1a{i}(:,subtractRange),2));
    %plot(gain(i),cf,'ro');
end
title('MEAN CORR');
%ylim([-0.7 0.7])
%xlim([-70 70])
hold off;

figure;
for i=1:nDataSets
    %subplot(1,2,1)
    plot(gain(i),std(meanMIX3a{i}),'bo')
    hold all;
    %plot(gain(i),std(nanmean(MIX3a{i}(:,subtractRange),2)),'ro');
    %subplot(1,2,2)
    %plot(gain(i),std(meanMIX3a{i})-std(nanmean(MIX3a{i}(:,subtractRange),2)),'bo');
    %hold all;
end
title('STD of MEANS')

figure;
for i=1:nDataSets
    dat = MIX3a{i}(:,sampleRange);
    dat = dat(:);
    plot(gain(i),nanstd(dat),'bo'); %-
    hold all;
    dat = MIX3a{i}(:,subtractRange);
    dat = dat(:);
    %plot(gain(i),nanstd(dat),'ro'); %-
end
title('STD of SAMPLES');

figure;
sampleRange = 100:150; % 140:165 % JACK CHANGED RANGE
%subtractRange = 70:110;
for s=sampleRange
    subplot(1,2,1)
    for i=1:nDataSets
        plot(gain(i),winstdMIX3a{i}(s),'bo'); %-nanmean(winstdMIX3a{i}(subtractRange))
        hold all;
    end
    title(sprintf('STD %d',s));
    %ylim([0.5 1.3])
    xlim([-70 70])
    hold off;
    
    subplot(1,2,2)
    for i=1:nDataSets
        [cf,pval] = corr(MIX3a{i}(:,s), MIX1a{i}(:,s));
        %cf = cf - corr(nanmean(MIX3a{i}(:,subtractRange),2), nanmean(MIX1a{i}(:,subtractRange),2));
        plot(gain(i),cf,'bo');
        hold all;
    end
    title(sprintf('CORR %d',s));
    %ylim([-0.7 0.7])
    xlim([-70 70])
    hold off;
    
    input('continue...');
end

%%
gain = [0.5,0.4,0.3,0.2,-0.2];
saveNames = {'StaticWiggle_NoFF','TimeWiggle_NoFF', 'TimeWiggle_FFGain+30','TimeWiggle_FFGain-30'};
order = 1:nDataSets;
iter = 1;
figure;
for i=order
  %subplot(1,nDataSets,iter);
  iter = iter+1;
  scatter(meanMIX1a{i}, meanMIX3a{i},'o');
  hold all;
  
%   p = polyfit(meanMIX1a{i},meanMIX3a{i},1);
%   x = 2:10;
%   y = polyval(p,x);
%   plot(x,y,'k','LineWidth',1.5)
%   
  [cf,pval] = corr(meanMIX3a{i}, meanMIX1a{i});
  
  title(sprintf('%s: correl = %.2f',rawNames{i},cf));
  xlabel('Phase CT [degrees]')
  ylabel('Phase CLEX [degrees]')
  xlim([-20 0]);
  ylim([-10 25]);
  grid on;
  hold off;
  set(gca,'FontSize',18)
  set(findall(gcf,'type','text'),'FontSize',18)
  
  saveStr = sprintf('/home/jack/PhaseFeedforward/Analysis/20141217_Wiggle3_%s',saveNames{i});
  print([saveStr '.png'],'-dpng');
  savefig([saveStr '.fig']);
end

% figure;
% for i=1:nDataSets
%     p = polyfit(meanMIX1a{i},meanMIX3a{i},1);
%     plot(gain(i),p(2),'bo');
%     hold all;
% end

%%
figure
normSamp = 50;%98; %90
for i=1:nDataSets
  tmp3 = nanmean(MIX3a{i});
  tmp3 = tmp3-tmp3(normSamp);
  timeAxis = 4*(1:length(tmp3))-4
  plot(timeAxis,tmp3,'-','Color',kolors{i},'LineWidth',2);
  hold all;
  %plot(timeAxis,tmp3-tmp1,'Color',kolors{i},'LineWidth',2);
  %shadedErrorBar(1:length(MIX3a{i}),nanmean(MIX3a{i}),nanstd(MIX3a{i}),{'Color',[rand() rand() rand()]},1)
  hold all;
end
legend(rawNames,'Location','SouthEast');
for i=1:nDataSets
    tmp1 = nanmean(MIX1a{i});
    tmp1 = tmp1-tmp1(normSamp);
    timeAxis = 4*(1:length(tmp1))-4
    plot(timeAxis,tmp1,'--','Color',kolors{i},'LineWidth',2);
end
xlim([timeAxis(40) timeAxis(40+160)])
set(gcf,'units','normalized','position',[0 0 0.55 0.5])
yAxLim = get(gca,'YLim');
plot([timeAxis(min(sampleRange)) timeAxis(min(sampleRange))],yAxLim,'k','LineWidth',2);
plot([timeAxis(max(sampleRange)) timeAxis(max(sampleRange))],yAxLim,'k','LineWidth',2);
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


%%
figure
legenda = dataSetNames;
stdwinsize = 0; %1
for i=1:nDataSets
   winstdMIX3a{i} = zeros(201,1);
   winstdMIX1a{i} = zeros(201,1);
   for j=40:250-stdwinsize-1
     winstdMIX3a{i}(j,1) = std2(MIX3a{i}(:,j:j+stdwinsize));
     winstdMIX1a{i}(j,1) = std2(MIX1a{i}(:,j:j+stdwinsize));
   end
  timeAxis = 4*(1:length(winstdMIX3a{i}))-4
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
yAxLim = get(gca,'YLim');
plot([timeAxis(min(sampleRange)) timeAxis(min(sampleRange))],yAxLim,'k','LineWidth',2);
plot([timeAxis(max(sampleRange)) timeAxis(max(sampleRange))],yAxLim,'k','LineWidth',2);
xlim([timeAxis(40) timeAxis(40+160)])
title({dataDescription, DataSetBaseTitle});
xlabel('Time [ns]')
ylabel('Phase Jitter [degrees]');
grid on;
set(gcf,'units','normalized','position',[0 0 0.55 0.5])
set(gca,'FontSize',16)
set(findall(gcf,'type','text'),'FontSize',16)

%%
saveNames = {'StaticWiggle_NoFF','TimeWiggle_NoFF', 'TimeWiggle_FFGain+30','TimeWiggle_FFGain-30'};

figure; 
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
  ylim([-10 10]);
  xlim([0 250]);
  
  xlabel('Time [Pulse No.]');
  ylabel('Mean Phase [degrees]');
  
  set(gca,'FontSize',16)
  set(findall(gcf,'type','text'),'FontSize',16)

%   saveStr = sprintf('/home/jack/PhaseFeedforward/Analysis/20141217_Wiggle3_MeanPhaseVsTime_%s',saveNames{i});
%   print([saveStr '.png'],'-dpng');
%   savefig([saveStr '.fig']);
%   hold off;
end

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

figure;
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

end
