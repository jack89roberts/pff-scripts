
%% with localtimestamp
auxSize = length(myDat);
auxTiming = zeros(7,auxSize);
for i=2:auxSize
    auxTiming(:,i) = myDat(i).localTimeStamp_ms - myDat(i-1).localTimeStamp_ms + (1:7)*1000;
end

plot(auxTiming')

%% 
auxSize = length(myDat);

meanPulse = NaN(1,auxSize);
timeStamp = NaN(1,auxSize);
for i=2:auxSize
    meanPulse(i) = nanmean(myDat(i).CT.SVBPI0608S.Samples.samples.value(220:290))-...
        nanmean(myDat(i-1).CT.SVBPI0608S.Samples.samples.value(220:290));
    timeStamp(i) = (myDat(i).CT.SVBPI0608S.Samples.samples.timeStamp-myDat(i-1).CT.SVBPI0608S.Samples.samples.timeStamp)/1000000;
end
%%
plot(meanPulse==0)
hold on
%plot(timeStamp,'r-')
hold off
grid
%%
plot(myDat(60).CT.SVBPI0608S.Samples.samples.value,'x-')
hold on
plot(myDat(61).CT.SVBPI0608S.Samples.samples.value,'r.-')
hold off