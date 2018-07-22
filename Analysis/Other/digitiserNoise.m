sisDat=loadMergedData('20151120_1538_Gain-800_R56_0.1_Interleaaved');
fontDat=saveExtractedFONTData('1540_PFFOptics_150414','/home/jack/PhaseFeedforward/FONTData/20140415');
%%
sisADC = double(sisDat(10).CT_SCOPE01_CH02.Acquisition.value.value);
sisSens = sisDat(10).CT_SCOPE01_CH04.Acquisition.sensitivity.value;
sisOff = sisDat(10).CT_SCOPE01_CH04.Acquisition.offset.value;
sisADC = sisADC.*sisSens - sisOff;
sisADC = sisADC-mean(sisADC);
sisCntPerVolt = 2048/2.5;
sisADC = sisADC*sisCntPerVolt;

fontADC = squeeze(fontDat.ADCs(1,178,2:end));
fontADC = fontADC-mean(fontADC);
% fontVoltPerCnt = 0.5/4096;
% fontADC = fontADC.*fontVoltPerCnt;

sisTime = ((1:length(sisADC))-1).*sisDat(10).CT_SCOPE01_CH02.Acquisition.sampleInterval.value;
fontTime = ((1:length(fontADC))-1).*(1000/357);
%%
figure;
plot(sisTime,sisADC,'b','LineWidth',1.5);
[~,stdSiS,~,stdSiS_err]=nanMeanStdErr(sisADC*1000)

hold all;
plot(fontTime,fontADC,'r','LineWidth',1.5);
[~,stdFONT,~,stdFONT_err]=nanMeanStdErr(fontADC*1000)

xlabel('Time [ns]')
ylabel('Output [V]')
title('Digitiser Noise')
xlim([147.3 1971.3]);
legend('SiS Digitisers','FONT5a Board')
format_plots;


%%

adcJit=1.47;
adcErr=0.04;
calA=2870;
calErr=43;

clc
phasJit = asind(adcJit/calA)
phasErr = (180/pi)*sqrt( ( (adcJit*calErr)^2 + (calA*adcErr)^2 ) / ((calA^4))*(1- (adcJit/calA)^2) )
