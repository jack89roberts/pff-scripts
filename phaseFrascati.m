% frascatiMon1_Diode = lastDataStruct.CT.SCOPE01.CH01.Acquisition.value.value;
% frascatiMon1_Mixer = lastDataStruct.CT.SCOPE01.CH02.Acquisition.value.value;
% 
% frascatiMon1_DiodeSensitivity = lastDataStruct.CT.SCOPE01.CH01.Acquisition.sensitivity.value;
% frascatiMon1_MixerSensitivity = lastDataStruct.CT.SCOPE01.CH02.Acquisition.sensitivity.value;
% 
% frascatiMon1_CalibrationConstant = 1.4299;
% 
% frascatiMon1_Phase = getPhaseMixerDiode(frascatiMon1_Mixer, frascatiMon1_Diode, frascatiMon1_MixerSensitivity, frascatiMon1_DiodeSensitivity, frascatiMon1_CalibrationConstant);
% 
% %subplot(4,2,3);
% subplot(1,3,1);
% plot(frascatiMon1_Phase);
% title('Frascati Monitor (CT)');
% ylim([-40 10]);
% %hold all;

mon1Mix = double(lastDataStruct.CT_SCOPE01_CH02.Acquisition.value.value).*lastDataStruct.CT_SCOPE01_CH02.Acquisition.sensitivity.value;
mon1Dio = double(lastDataStruct.CT_SCOPE01_CH01.Acquisition.value.value).*lastDataStruct.CT_SCOPE01_CH01.Acquisition.sensitivity.value;
mon2Mix = double(lastDataStruct.CT_SCOPE01_CH04.Acquisition.value.value).*lastDataStruct.CT_SCOPE01_CH04.Acquisition.sensitivity.value;
mon2Dio = double(lastDataStruct.CT_SCOPE01_CH03.Acquisition.value.value).*lastDataStruct.CT_SCOPE01_CH03.Acquisition.sensitivity.value;
mon3Mix = double(lastDataStruct.CT_SCOPE01_CH06.Acquisition.value.value).*lastDataStruct.CT_SCOPE01_CH06.Acquisition.sensitivity.value;
mon3Dio = double(lastDataStruct.CT_SCOPE01_CH05.Acquisition.value.value).*lastDataStruct.CT_SCOPE01_CH05.Acquisition.sensitivity.value;


% subplot(1,2,1)
% plot(mon1Mix)
% hold all;
% plot(mon2Mix)
% plot(mon3Mix)
% legend('Mon1Mix (CH2)','Mon2Mix (CH4)','Mon3Mix (CH6')
% hold off;
% subplot(1,2,2)
% plot(mon1Dio)
% hold all;
% plot(mon2Dio)
% plot(mon3Dio);
% legend('Mon1Dio (CH1)','Mon2Dio (CH3)','Mon3Dio (CH5')
% hold off;

mon1Phase = getPhaseMixerDiode(mon1Mix,[],0.436051,0.046733);
mon2Phase = getPhaseMixerDiode(mon2Mix,[],0.513644,0.007459);

plot(mon1Phase);
hold all;
plot(mon2Phase);
hold off;

% mon1Mix = mon1Mix - mean(mon1Mix(1:50));
% mon1Dio = mon1Dio - mean(mon1Dio(1:50));
% mon2Mix = mon2Mix - mean(mon2Mix(1:50));
% mon2Dio = mon2Dio - mean(mon2Dio(1:50));
% mon3Mix = mon3Mix - mean(mon3Mix(1:50));
% mon3Dio = mon3Dio - mean(mon3Dio(1:50));
% 
% 
% phaseMon1 = getPhaseMixerDiode(mon1Mix,mon1Dio,1.211,0.06);
% phaseMon2 = getPhaseMixerDiode(mon2Mix,mon2Dio,0.200,0.06);
% phaseMon3 = getPhaseMixerDiode(mon3Mix,mon3Dio,0.354,0.06);
% 
% 
% 
% meanPhaseMon1 = mean(phaseMon1(140:240));
% meanPhaseMon2 = mean(phaseMon2(140:240));
% meanPhaseMon3 = mean(phaseMon3(300:400));
% 
% phaseMon1 = phaseMon1 - mean(phaseMon1(140:240));
% phaseMon2 = phaseMon2 - mean(phaseMon2(140:240));
% phaseMon3 = phaseMon3 - mean(phaseMon3(300:400));
% 
% 
% subplot(1,2,1)
% plot(phaseMon1);
% hold all;
% plot(phaseMon2);
% plot(phaseMon3);
% legend('Mon1','Mon2');%,'Mon3');
% hold off;
% subplot(1,2,2)
% plot(meanPhaseMon1,meanPhaseMon2,'bo');
% hold all;
% plot(meanPhaseMon1,meanPhaseMon3,'ro');
% legend('Mon1-Mon2','Mon3-Mon2');
% 

