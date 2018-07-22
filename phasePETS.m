tblScope3_I = lastDataStruct.CE.SCOPE03.CH01.Acquisition.value.value;
tblScope3_Q = lastDataStruct.CE.SCOPE03.CH02.Acquisition.value.value;

[tblScope3_Phase,tblScope3_Power] = getPhaseIQ(tblScope3_I, tblScope3_Q);


%subplot(4,2,6);
%subplot(1,2,1);
plot(tblScope3_Phase);
title('PETS Monitor TBL - Phase');
%ylim([-90 90])

% subplot(1,2,2);
% plot(tblScope3_Power);
% title('PETS Monitor TBL - Power');
% ylim([0 16000])
% %hold all;