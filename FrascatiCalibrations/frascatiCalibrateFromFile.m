close all; clearvars;
%% input file names, type
calTimeStamp = '20161215_1945';%'20151116_1651';
useAlignment = 1; % if true uses alignment function (all monitors can use same sample range). If false no alignment and must define sample ranges for each monitor.
monIndices = []; % if empty uses all monitors. Otherwise uses monitors defined here.
saveData = 1; % do you want to save the results or not (will overwrite any pre-existing files)

useOffset = 2; % set zero if continuous signal rather than pulse

isFONTCal = 1; % 0 if calibration using SiS data, 1 if calibration using FONT data
mixerADCs = []; % if empty assumes default ADCs for mixers
diodeIsPresent = 1; % 0 if diodes not present in data and should not be used for removal of bad pulses etc.
diodeADCs = []; % if empty assumes default ADCs for diodes (or does not use diodes if ~diodeIsPresent)
%% load calibration data
if (isFONTCal) % load FONT cal data
    addpath('../FONT/');
    fontBaseDir = ['/home/jack/PhaseFeedforward/FONTData/' calTimeStamp(1:6) '/Calibration/' calTimeStamp];
    fontCalDir = [fontBaseDir '/Extracted'];
    
    if (~exist(fontCalDir,'dir')) % extract the calibration files if they haven't been already
        filesToExtract = dir([fontBaseDir '/*.dat']);
        nFiles = length(filesToExtract);
        for f=1:nFiles
            saveExtractedFONTData( strrep(filesToExtract(f).name,'.dat',''),fontBaseDir );
        end
    end
    
    [ rawMixers, rawDiodes, scanPhShiftValues ] = packageFONTCalData(fontCalDir, mixerADCs, diodeIsPresent, diodeADCs);
    rawMixerSensitivities = '';
    rawDiodeSensitivities = '';
    [nMons,~,~,~] = size(rawMixers);
    savePath = [fontBaseDir '/Processed/'];
    if (~exist(savePath,'dir'))
        mkdir(savePath);
    end    
    phShiftNames = {'Mon1','Mon2','Mon3'};
    monMixerNames = {'Mon1','Mon2','Mon3'};

else % load SiS data
    if (strcmp(getenv('USER'),'jack') == 1) % If username is 'jack' assume environment is Jack's laptop
        dataDir = ['/home/jack/PhaseFeedforward/CTFData/' calTimeStamp(1:6) '/FrascatiCalibrations/'];
        dataName = [dataDir 'frascatiCalibration_' calTimeStamp '.mat']; 
        load(dataName);
        savePath = dataDir;
    else % otherwise assume control room environment
        %load(sprintf('data/frascatiCalibration_%s.mat',calTimeStamp));
        load(sprintf('/user/ctf3op/PhaseFeedforward/FrascatiCalibrations/frascatiCalibration_%s.mat',calTimeStamp));
        savePath = '/user/ctf3op/PhaseFeedforward/FrascatiCalibrations/';
    end

end

%% calibration parameters: comment anything you don't want to overwrite from loaded data

calSampleRange = [];%735:737;%[]% leave empty to prompt for sample range to use
%calSampleRange = {1:2,3:4,5:6};
useDiode = 0; % 0=mixerOnly, 1=mixer/sqrtDiode, 2=mixer/Diode
phaseShiftFreq = 4;

%scanPhShiftValues = 0:2:16;
%scanPhShiftValues = [scanPhShiftValues 16.75];
  
%%
savePlotName = sprintf('%s%s',savePath,calTimeStamp);

% process signals
fprintf('Processing signals...\n');
[alignedMixers, alignedDiodes, calSampleRange] = ...
    frascatiCalProcessSignals(...
        rawMixers,...
        rawDiodes,...
        rawMixerSensitivities,...
        rawDiodeSensitivities,...
        calSampleRange,...
        monIndices,...
        useAlignment,...
        useOffset);
fprintf('Finished processing signals!\n');   

% fit data
fprintf('Fitting data...\n');
[calibrationFactors, zeroCrossings, calRSq, calConstErr] = ...
    frascatiCalCalculateCalibrations(...
        scanPhShiftValues,...
        alignedMixers,...
        alignedDiodes,...
        calSampleRange,...
        useDiode,...
        phaseShiftFreq,...
        savePlotName);
 calConstErr = calConstErr/2; % conf interval to approx std error   

fprintf('Finished fitting data!\n');

% calculate optimal sensitivities - to be implemented

% save data
if (saveData)
    fprintf('Saving data... ')

    zeroCrossName = sprintf('%s%s%s',savePath,'frascatiZeroCrossings_',calTimeStamp);
    zeroCrossFileBackup = fopen(zeroCrossName,'w');
    calibConstName = sprintf('%s%s%s',savePath,'frascatiCalibrationConstants_',calTimeStamp);
    calibConstFileBackup = fopen(calibConstName,'w');
    for mon=1:nMons 
        fprintf(zeroCrossFileBackup,'%s, %.2f\n',...
            phShiftNames{mon},...
            zeroCrossings(mon)...
        );
        fprintf(calibConstFileBackup,'%s, %f, %f, %f, %f\n',...
            monMixerNames{mon},...
            calibrationFactors(mon,1),...
            calibrationFactors(mon,2),...
            calibrationFactors(mon,3),...
            calibrationFactors(mon,4)...
        );
    end
    fclose(zeroCrossFileBackup);

    fprintf(calibConstFileBackup, 'useDiode, %d', useDiode);
    fclose(calibConstFileBackup);
end

for mon=1:nMons
    fprintf('---------------------------------\n')
    fprintf('MON %d\n',mon)
    fprintf('---------------------------------\n')
    
    fprintf('A = %.3f %c %.3f\n',calibrationFactors(mon,1),char(177),calConstErr(mon,1));
    fprintf('b = %.4f %c %.4f\n',calibrationFactors(mon,2),char(177),calConstErr(mon,2));
    fprintf('c = %.3f %c %.3f\n',calibrationFactors(mon,3),char(177),calConstErr(mon,3));
    fprintf('d = %.3f %c %.3f\n',calibrationFactors(mon,4),char(177),calConstErr(mon,4));

end

fprintf('---------------------------------\n')
fprintf('done!\n');

% Set values in machine if requested
if (~strcmp(getenv('USER'),'jack')) % assume in the control room if not Jack's laptop
    setCalInMachine = input('Set these values as reference (0 or 1)?');
    if (setCalInMachine)
        fprintf('Setting new shifter values and saving reference files... ')

        % files for saving reference calibrations
        zeroCrossFile = fopen('frascatiZeroCrossings','w');
        calibConstFile = fopen('frascatiCalibrationConstants','w');

        refNameFile = fopen('frascatiRefCalName','w');
        fprintf(refNameFile,'%s',calTimeStamp);
        fclose(refNameFile);

        for mon=1:nMons
            % always want to set the zero crossing to the falling slope so check
            % this here
            dydxAtZero =...
                calibrationFactors(mon,1)*calibrationFactors(mon,2)*...
                cos(calibrationFactors(mon,2)*zeroCrossings(mon) + calibrationFactors(mon,3)); % Ab*cos(bx+c)
            if (dydxAtZero > 0)
                zeroCrossingToSet = zeroCrossings(mon) + (180/(12/phaseShiftFreq)); 
            else
                zeroCrossingToSet = zeroCrossings(mon);
            end
            
            if (strcmp(scanType,'auto'))
                setPhaseShifter(mon,zeroCrossingToSet);
            end
            
            fprintf(zeroCrossFile,'%s, %.2f\n',...
                phShiftNames{mon},...
                zeroCrossingToSet...
            );
            fprintf(calibConstFile,'%s, %f, %f, %f, %f\n',...
                monMixerNames{mon}, calibrationFactors(mon,1),...
                calibrationFactors(mon,2),...
                calibrationFactors(mon,3),...
                calibrationFactors(mon,4)...
            );

        end

        fclose(zeroCrossFile);

        fprintf(calibConstFile, 'useDiode, %d', useDiode);
        fclose(calibConstFile);

    end
    fprintf('done!\n')
end



