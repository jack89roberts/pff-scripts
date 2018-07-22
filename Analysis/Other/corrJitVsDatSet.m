clearvars; close all
%% 201612 scan

datFFOff = {
'/home/jack/PhaseFeedforward/Analysis/201612/201612_1800_straightIntlv_gain_-400_091216_FFOff.mat'

'/home/jack/PhaseFeedforward/Analysis/201612/201612_1800_straightIntlv_gain_-400_091216_FFOff.mat'
'/home/jack/PhaseFeedforward/Analysis/201612/201612_1745_straightIntlv_gain_-600_091216_FFOff.mat'
'/home/jack/PhaseFeedforward/Analysis/201612/201612_1745_straightIntlv_gain_-800_091216_FFOff.mat'
'/home/jack/PhaseFeedforward/Analysis/201612/201612_1745_straightIntlv_gain_-1000_091216_FFOff.mat'
'/home/jack/PhaseFeedforward/Analysis/201612/201612_1745_straightIntlv_gain_-1200_091216_FFOff.mat'
'/home/jack/PhaseFeedforward/Analysis/201612/201612_1745_straightIntlv_gain_-1400_091216_FFOff.mat'
'/home/jack/PhaseFeedforward/Analysis/201612/201612_1800_straightIntlv_gain_-1600_091216_FFOff.mat'
};

datFFOn = {
    '/home/jack/PhaseFeedforward/Analysis/201612/201612_1745_straightIntlv_gain_-600_091216_FFOff.mat'


'/home/jack/PhaseFeedforward/Analysis/201612/201612_1800_straightIntlv_gain_-400_091216_FFOn.mat'
'/home/jack/PhaseFeedforward/Analysis/201612/201612_1745_straightIntlv_gain_-600_091216_FFOn.mat'
'/home/jack/PhaseFeedforward/Analysis/201612/201612_1745_straightIntlv_gain_-800_091216_FFOn.mat'
'/home/jack/PhaseFeedforward/Analysis/201612/201612_1745_straightIntlv_gain_-1000_091216_FFOn.mat'
'/home/jack/PhaseFeedforward/Analysis/201612/201612_1745_straightIntlv_gain_-1200_091216_FFOn.mat'
'/home/jack/PhaseFeedforward/Analysis/201612/201612_1745_straightIntlv_gain_-1400_091216_FFOn.mat'
'/home/jack/PhaseFeedforward/Analysis/201612/201612_1800_straightIntlv_gain_-1600_091216_FFOn.mat'
};

upMon = 1;

gain = [0 400 600 800 1000 1200 1400 1600];
gainScale = (1/624)*0.7;

%% 201608 scan

% datFFOff = {...
% % off
% '/home/jack/PhaseFeedforward/Analysis/201608/20160811_1658_GScan_-2500_Intlv_Odd.mat'
% '/home/jack/PhaseFeedforward/Analysis/201608/20160811_1704_GScan_-2000_Intlv_Odd.mat'
% '/home/jack/PhaseFeedforward/Analysis/201608/20160811_1709_GScan_-1500_Intlv_Even.mat'
% '/home/jack/PhaseFeedforward/Analysis/201608/20160811_1714_GScan_-1000_Intlv_Odd.mat'
% '/home/jack/PhaseFeedforward/Analysis/201608/20160811_1720_GScan_-500_Intlv_Even.mat'
% '/home/jack/PhaseFeedforward/Analysis/201608/20160811_1803_GScan_1250_Intlv_Odd.mat'
% '/home/jack/PhaseFeedforward/Analysis/201608/20160811_1725_GScan_500_Intlv_Odd.mat'
% '/home/jack/PhaseFeedforward/Analysis/201608/20160811_1759_GScan_750_Intlv_Even.mat'
% '/home/jack/PhaseFeedforward/Analysis/201608/20160811_1733_GScan_1000_Intlv_Odd.mat'
% '/home/jack/PhaseFeedforward/Analysis/201608/20160811_1803_GScan_1250_Intlv_Odd.mat'
% '/home/jack/PhaseFeedforward/Analysis/201608/20160811_1738_GScan_1500_Intlv_Even.mat'
% '/home/jack/PhaseFeedforward/Analysis/201608/20160811_1743_GScan_2000_Intlv_Even.mat'
% '/home/jack/PhaseFeedforward/Analysis/201608/20160811_1751_GScan_2500_Intlv_Odd.mat'
% };
% 
% datFFOn = {...
% '/home/jack/PhaseFeedforward/Analysis/201608/20160811_1658_GScan_-2500_Intlv_Even.mat'
% '/home/jack/PhaseFeedforward/Analysis/201608/20160811_1704_GScan_-2000_Intlv_Even.mat'
% '/home/jack/PhaseFeedforward/Analysis/201608/20160811_1709_GScan_-1500_Intlv_Odd.mat'
% '/home/jack/PhaseFeedforward/Analysis/201608/20160811_1714_GScan_-1000_Intlv_Even.mat'
% '/home/jack/PhaseFeedforward/Analysis/201608/20160811_1720_GScan_-500_Intlv_Odd.mat'
% '/home/jack/PhaseFeedforward/Analysis/201608/20160811_1759_GScan_750_Intlv_Even.mat'
% '/home/jack/PhaseFeedforward/Analysis/201608/20160811_1725_GScan_500_Intlv_Even.mat'
% '/home/jack/PhaseFeedforward/Analysis/201608/20160811_1759_GScan_750_Intlv_Odd.mat'
% '/home/jack/PhaseFeedforward/Analysis/201608/20160811_1733_GScan_1000_Intlv_Even.mat'
% '/home/jack/PhaseFeedforward/Analysis/201608/20160811_1803_GScan_1250_Intlv_Even.mat'
% '/home/jack/PhaseFeedforward/Analysis/201608/20160811_1738_GScan_1500_Intlv_Odd.mat'
% '/home/jack/PhaseFeedforward/Analysis/201608/20160811_1743_GScan_2000_Intlv_Odd.mat'
% '/home/jack/PhaseFeedforward/Analysis/201608/20160811_1751_GScan_2500_Intlv_Even.mat'
% };
% 
% upMon = 2;
% 
% gain = [-2500 -2000 -1500 -1000 -500 0 500 750 1000 1250 1500 2000 2500];
% gainScale = (1/624)*0.7;

%%
gain = gain*gainScale;
nDatSets = length(datFFOff);

stdFFOff = NaN(3,nDatSets);
stdFFOn = NaN(3,nDatSets);
stdFFOff_err = NaN(3,nDatSets);
stdFFOn_err = NaN(3,nDatSets);
corrFFOff = NaN(1,nDatSets);
corrFFOn = NaN(1,nDatSets);
corrFFOff_err = NaN(1,nDatSets);
corrFFOn_err = NaN(1,nDatSets);


for i=1:nDatSets
    myData=load(datFFOff{i});
    stdFFOff(:,i) = myData.stdMeanPulsePhase;
    stdFFOff_err(:,i) = myData.stdMeanPulsePhase_err;
    if upMon==1
        corrFFOff(i) = myData.corrMeanMix1Mix3;
        corrFFOff_err(i) = myData.corrMeanMix1Mix3_err;
    else
        corrFFOff(i) = myData.corrMeanMix2Mix3;
        corrFFOff_err(i) = myData.corrMeanMix2Mix3_err;
    end
    
    myData=load(datFFOn{i});
    stdFFOn(:,i) = myData.stdMeanPulsePhase;
    stdFFOn_err(:,i) = myData.stdMeanPulsePhase_err;
    if upMon==1
        corrFFOn(i) = myData.corrMeanMix1Mix3;
        corrFFOn_err(i) = myData.corrMeanMix1Mix3_err;
    else
        corrFFOn(i) = myData.corrMeanMix2Mix3;
        corrFFOn_err(i) = myData.corrMeanMix2Mix3_err;
    end
end
clear myData

%% theoretical

if upMon==1
    stdUpOff = stdFFOff(1,:);
    stdUpOn = stdFFOn(1,:);
    stdUpOff_err = stdFFOff_err(1,:);
    stdUpOn_err = stdFFOn_err(1,:);
else
    stdUpOff = stdFFOff(2,:);
    stdUpOn = stdFFOn(2,:);
    stdUpOff_err = stdFFOff_err(2,:);
    stdUpOn_err = stdFFOn_err(2,:);
end
stdDownOff = stdFFOff(3,:);
stdDownOn = stdFFOn(3,:);
stdDownOff_err = stdFFOff_err(3,:);
stdDownOn_err = stdFFOn_err(3,:);

stdOpt = stdUpOff.*sqrt(1-corrFFOff.^2);
stdTheoret = sqrt(...
                stdDownOff.^2 ...
                +(gain.^2).*(stdUpOff.^2) ...
                -2.*gain.*corrFFOff.*stdUpOff.*stdDownOff...
                );

stdTheoret_err = sqrt(...
                    ((2*stdDownOff - 2*gain.*corrFFOff.*stdUpOff).*stdDownOff_err).^2 ...
                   +((2*(gain.^2).*stdUpOff - 2*gain.*corrFFOff.*stdDownOff).*stdUpOff_err).^2 ...
                   +((-2*gain.*stdUpOff.*stdDownOff).*corrFFOff_err).^2 ...
                 ); % error on std
 stdTheoret_err = stdTheoret_err./(2*stdTheoret); % error std^2 to error on std
                   
%%
figure;
errorbar(gain,stdFFOff(3,:),stdFFOff_err(3,:),'bo','MarkerFaceColor','b');
hold all
errorbar(gain,stdFFOn(3,:),stdFFOff_err(3,:),'ro','MarkerFaceColor','r');
title('std')

figure;
errorbar(gain,corrFFOff,corrFFOff_err,'bo','MarkerFaceColor','b');
hold all
errorbar(gain,corrFFOn,corrFFOn_err,'ro','MarkerFaceColor','r');
title('corr')

figure;
theLine=shadedErrorBar(gain,stdTheoret,stdTheoret_err,'r');
hold all
datLine=errorbar(gain,stdFFOn(3,:),stdFFOff_err(3,:),'bo','MarkerFaceColor','b');
legend([datLine theLine.patch],{'dat','the'})
