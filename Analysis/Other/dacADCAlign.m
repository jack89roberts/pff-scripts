%% look at adc vs dac

datFile=
gain=400;

%params
downOffset=-235; % downstream samps wrt upstream samps (before avg)
dacOffset=-1; % dac samps wrt upstream samps
pulseRange=80:280; % pulse range on dac (after dacOffset)
dacConv=gain/65;
dacSign=-1;


%%
[ FONTDataFFOn, FONTDataFFOff ] = splitInterleavedFONTData( FONTData );


%extract
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
avgUp=ceil(avgUp);
avgDown=ceil(avgDown);
avgDOff=ceil(avgDOff);


figure;
plot(onDac'/dacConv,'r');
hold all;
plot(avgUp','b');

figure;
plot((onDac./avgUp)')

figure;
plot(avgUp','b')
hold all
plot(onDac'/dacConv,'r');
plot(avgDOff','g');
plot(avgDown','m');

figure;
plot(avgUp(:,pulseRange),onDac(:,pulseRange),'o')

notSatDac=onDac(:,pulseRange);
notSatUp=avgUp(:,pulseRange);
notSatDac=notSatDac(:);
notSatUp=notSatUp(:);
satDac=notSatDac>4095;
notSatDac=notSatDac(~satDac);
notSatUp=notSatUp(~satDac);
figure
plot(notSatUp,notSatDac,'o')
estGain=nanpolyfit(notSatDac,notSatUp,1);
estGain=64/estGain(1)

figure;
plot(onDac(:,pulseRange)',avgDown(:,pulseRange)'-avgDOff(:,pulseRange)','o');
%%
