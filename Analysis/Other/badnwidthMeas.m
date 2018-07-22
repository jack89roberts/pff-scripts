dataSetNames = {...
    '20150731_1444_BandwidthMeas',...
    '20150731_1504_BandwidthMeas',...
    '20150731_1520_BandwidthMeas',...
    '20150731_1526_BandwidthMeas',...
    '20150731_1534_BandwidthMeas',...
    '20150731_1541_BandwidthMeas',...
    '20150731_1554_BandwidthMeas',...
    '20150731_1602_BandwidthMeas',...
    '20150731_1607_BandwidthMeas',...
    '20150731_1620_BandwidthMeas',...
    '20150731_1624_BandwidthMeas',...
    '20150731_1628_BandwidthMeas',...
    '20150731_1630_BandwidthMeas',...
    '20150731_1633_BandwidthMeas',...
    '20150731_1637_BandwidthMeas',...
    '20150731_1640_BandwidthMeas',...
    '20150731_1644_BandwidthMeas'...
};

%%
nDataSets = length(dataSetNames);
downstreamPulses = NaN(nDataSets,3,384);
for i=1:nDataSets
    [CTFData,~,~] = loadMergedData( dataSetNames{i} );
    downstreamPulses(i,1,:) = CTFData(1).CT_SCOPE01_CH02.Acquisition.value.value;
    downstreamPulses(i,2,:) = CTFData(1).CT_SCOPE01_CH04.Acquisition.value.value;
    downstreamPulses(i,3,:) = CTFData(1).CT_SCOPE01_CH06.Acquisition.value.value;

end

figure;
for i=1:nDataSets
    plot(squeeze(downstreamPulses(i,2,:)));
    hold all;
end

%%
datSets={
    '/home/jack/PhaseFeedforward/Analysis/201507/Bandwidth/20150731_1444_BandwidthMeas.mat'
    '/home/jack/PhaseFeedforward/Analysis/201507/Bandwidth/20150731_1504_BandwidthMeas.mat'
    '/home/jack/PhaseFeedforward/Analysis/201507/Bandwidth/20150731_1520_BandwidthMeas.mat'
    '/home/jack/PhaseFeedforward/Analysis/201507/Bandwidth/20150731_1526_BandwidthMeas.mat'
    '/home/jack/PhaseFeedforward/Analysis/201507/Bandwidth/20150731_1534_BandwidthMeas.mat'
    '/home/jack/PhaseFeedforward/Analysis/201507/Bandwidth/20150731_1541_BandwidthMeas.mat'
    '/home/jack/PhaseFeedforward/Analysis/201507/Bandwidth/20150731_1554_BandwidthMeas.mat'
    '/home/jack/PhaseFeedforward/Analysis/201507/Bandwidth/20150731_1602_BandwidthMeas.mat'
    '/home/jack/PhaseFeedforward/Analysis/201507/Bandwidth/20150731_1607_BandwidthMeas.mat'
    '/home/jack/PhaseFeedforward/Analysis/201507/Bandwidth/20150731_1620_BandwidthMeas.mat'
    '/home/jack/PhaseFeedforward/Analysis/201507/Bandwidth/20150731_1624_BandwidthMeas.mat'
    '/home/jack/PhaseFeedforward/Analysis/201507/Bandwidth/20150731_1628_BandwidthMeas.mat'
    '/home/jack/PhaseFeedforward/Analysis/201507/Bandwidth/20150731_1630_BandwidthMeas.mat'
    '/home/jack/PhaseFeedforward/Analysis/201507/Bandwidth/20150731_1633_BandwidthMeas.mat'
    '/home/jack/PhaseFeedforward/Analysis/201507/Bandwidth/20150731_1637_BandwidthMeas.mat'
    '/home/jack/PhaseFeedforward/Analysis/201507/Bandwidth/20150731_1640_BandwidthMeas.mat'
    '/home/jack/PhaseFeedforward/Analysis/201507/Bandwidth/20150731_1644_BandwidthMeas.mat'
};

% datSets={
%     '20150731_1444_BandwidthMeas'
%     '20150731_1504_BandwidthMeas'
%     '20150731_1520_BandwidthMeas'
%     '20150731_1526_BandwidthMeas'
%     '20150731_1534_BandwidthMeas'
%     '20150731_1541_BandwidthMeas'
%     '20150731_1554_BandwidthMeas'
%     '20150731_1602_BandwidthMeas'
%     '20150731_1607_BandwidthMeas'
%     '20150731_1620_BandwidthMeas'
%     '20150731_1624_BandwidthMeas'
%     '20150731_1628_BandwidthMeas'
%     '20150731_1630_BandwidthMeas'
%     '20150731_1633_BandwidthMeas'
%     '20150731_1637_BandwidthMeas'
%     '20150731_1640_BandwidthMeas'
%     '20150731_1644_BandwidthMeas'
% };

figure;
for datSet=1:length(datSets)
    load(datSets{datSet});
%     CTFData=loadMergedData(datSets{i});
    
%     plot(CTFData(10).CT_SVBPI0608S.Samples.samples.value)
%     hold all
    subplot(1,2,1)
    plot(squeeze(nanmean(mixers(3,:,:),2)));
    subplot(1,2,2)
    plot(meanBPMSAlongPulse{14})
    input(sprintf('%d...',datSet))
end

% figure;
% for b=1:nBPMs
%     try
%         a=figure;
%         plot(((1:bpmNSamples(b))-bpmPulseRange{b}(1)).*bpmTimePerSample(b),meanBPMSAlongPulse{b})
%         title(bpmNames{b})
%         hold all
%     catch
%         close(a)
%     end
% end

figure;
for datSet=1:length(datSets)
    load(datSets{datSet});
    plot(meanBPMSAlongPulse{14})
    hold all
end