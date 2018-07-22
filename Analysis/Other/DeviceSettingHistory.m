addpath('../../')

[dev,val,time]=loadRefSettings('/home/jack/PhaseFeedforward/CTFData/201511/CTF3Monitor/20151120');

%%
meanVal = nanmean(val);
normVal = val;
for d=1:length(dev)
    normVal(:,d) = normVal(:,d)-meanVal(d);
    normVal(:,d) = normVal(:,d)./max(abs(normVal(:,d)));
end

stdDev = nanstd(val);
changedDevs = dev(stdDev>0);

plot(normVal(:,stdDev>0))

startTime = datenum(2015,11,20,15,00,00);
endTime = datenum(2015,11,20,18,30,00);

filtInds = (startTime<time)&(time<endTime);
filtTime = time(filtInds);
filtVal = val(filtInds,:);
filtNormVal = normVal(filtInds,:);

filtStdDev = nanstd(filtVal);
filtChangedInds = find(filtStdDev>0);
filtChangedDevs = dev(filtChangedInds)

for i=1:length(filtChangedDevs)
    plot(filtTime,filtVal(:,filtChangedInds(i)));
    title(filtChangedDevs{i})
    datetick('x');
    input('Next plot...')
end