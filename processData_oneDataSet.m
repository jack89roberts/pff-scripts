dataDir = '/home/jack/PhaseFeedforward/CTFData/201412/20141205_1125_R56_+0.4_TimingScans_FF_K1Gain_0_K2Gain_0_noGate';
frascatiCal = '/home/jack/PhaseFeedforward/CTFData/201412/frascatiCalibrations/frascatiCalibrationConstants_20141205_0944';
saveIfBPMDataBad = 0;
saveIfPhaseDataBad = 0;
%%
dataName = regexp(dataDir,'/','split');
dataName = dataName{length(dataName)};

fprintf(1,'-----------------------------------\n');
fprintf('DATA SET: %s\n',dataName);
fprintf(1,'-----------------------------------\n');

processedData = processData(dataDir,frascatiCal,saveIfBPMDataBad,saveIfPhaseDataBad);

if (~isempty(processedData))
    fprintf(1,'Saving data from %s\n', dataDir);
    saveDir = [dataDir '/processed/'];
    if (~exist(saveDir,'dir'))
        mkdir(saveDir);
    end

    save([saveDir dataName '.mat'],'processedData');
else
    fprintf(1,'All files in %s were bad, nothing saved.\n', dataDir);
end
