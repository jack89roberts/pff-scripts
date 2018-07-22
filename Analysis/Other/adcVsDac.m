clearvars;close all;
%% look at adc vs dac

datFile='/home/jack/PhaseFeedforward/FONTData/201612/Extracted/201612_1800_straightIntlv_gain_-1600_091216.mat';%'/home/jack/PhaseFeedforward/FONTData/201612/Extracted/201612_1800_straightIntlv_gain_-1600_091216.mat';%'/home/jack/PhaseFeedforward/FONTData/201612/Extracted/201612_1800_straightIntlv_gain_-400_091216.mat';
gain=1600;%400;

%params
downOffset=-235; % downstream samps wrt upstream samps (before avg)
dacOffset=-1; % dac samps wrt upstream samps
pulseRange=83:277; % pulse range on dac (after dacOffset)
dacSign=-1;


%%
%load
load(datFile);

%extract
[ FONTDataFFOn, FONTDataFFOff ] = splitInterleavedFONTData( FONTData );
onUp=squeeze(FONTDataFFOn.ADCs(2,:,:));
onDown=squeeze(FONTDataFFOn.ADCs(6,:,:));
offDown=squeeze(FONTDataFFOff.ADCs(6,:,:));
onDac=dacSign*squeeze(FONTDataFFOn.DACs(1,:,:));
[~,nPulses,~]=size(FONTDataFFOn.ADCs);

% align all sigs with upstream
onDown=delaySignal(onDown,downOffset);
offDown=delaySignal(offDown,downOffset);
onDac=delaySignal(onDac,dacOffset);

% average adc samples to match dac
avgUp=NaN(nPulses,450);
avgDown=NaN(nPulses,450);
avgDOff=NaN(nPulses,450);

sSamp = 1:2:899;
for i=1:450
    avgUp(:,i)=(onUp(:,sSamp(i))+onUp(:,sSamp(i)+1))/2;
    avgDown(:,i)=(onDown(:,sSamp(i))+onDown(:,sSamp(i)+1))/2;
    avgDOff(:,i)=(offDown(:,sSamp(i))+offDown(:,sSamp(i)+1))/2;
end
% avgUp=ceil(avgUp);
% avgDown=ceil(avgDown);
% avgDOff=ceil(avgDOff);

avgUp=removeBadPulses(avgUp,pulseRange);
avgDown=removeBadPulses(avgDown,pulseRange);
avgDOff=removeBadPulses(avgDOff,pulseRange);
onDac=removeBadPulses(onDac,pulseRange);

dacConv=gain/64;

figure;
plot(avgUp','b')
hold all
plot(onDac'/dacConv,'r');
plot(avgDOff','g');
plot(avgDown','m');

figure;
plot(avgUp(:,pulseRange),onDac(:,pulseRange),'o')

% some arrays without saturated dac points
notSatDac=onDac(:,pulseRange);
notSatUp=avgUp(:,pulseRange);
notSatDown=avgDown(:,pulseRange);
notSatDOff=avgDOff(:,pulseRange);
notSatDac=notSatDac(:);
notSatUp=notSatUp(:);
notSatDown=notSatDown(:);
notSatDOff=notSatDOff(:);

satDac=notSatDac>4095 | notSatDac<-4094;
notSatDac=notSatDac(~satDac);
notSatUp=notSatUp(~satDac);
notSatDown=notSatDown(~satDac);
notSatDOff=notSatDOff(~satDac);

figure
plot(notSatUp,notSatDac,'o')
estGain=nanpolyfit(notSatDac,notSatUp,1);
estGain=64/estGain(1)

figure;
plot(onDac(:,pulseRange)',avgDown(:,pulseRange)'-avgDOff(:,pulseRange)','o');

figure;
plot(avgDown'-avgDOff');

% figure;
% plot((-avgDown(2,:)+avgDOff(2,:))*dacConv);
% hold all
% plot(onDac(2,:))
% plot(avgUp(2,:)*dacConv)

figure;
plot(notSatDac,notSatDown-notSatDOff,'o');

% sort, bin plot adc vs dac
[sortDac,order]=sort(notSatDac);
sortDDiff=notSatDown-notSatDOff;
sortDDiff=sortDDiff(order);

nBins=50;
binCount=round(length(sortDDiff)/nBins);
binDac=NaN(1,nBins);
binDDiff=NaN(1,nBins);
binDac_err=NaN(1,nBins);
binDDiff_err=NaN(1,nBins);
for i=1:nBins
    startRange = ((i-1)*binCount)+1;
    if i==nBins
        myRange=startRange:length(sortDDiff);
    else
        myRange=startRange:(startRange+binCount-1);
    end
    [binDac(i),~,binDac_err(i)]=nanMeanStdErr(sortDac(myRange));
    [binDDiff(i),~,binDDiff_err(i)]=nanMeanStdErr(sortDDiff(myRange));
end
figure;
errorbar(binDac,-binDDiff,binDDiff_err,'bo');
hold all;
herrorbar(binDac,-binDDiff,binDDiff_err,'bo');

% pulse mean adc vs. dac
x=onDac(:,180);%nanmean(onDac(:,180),2);
y=onUp(:,180);%nanmean(avgDown(:,180)-avgDOff(:,180),2);
figure;
plot(x,y,'o')