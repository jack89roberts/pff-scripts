%% processData.m - May 2015
% General data processing/analysis script for creating stats and
% correlations from the Fracati monitors and BPMs.

close all; clearvars;

%% quick hack to loop through datasets
% In current implementation, only one sample range etc. defined that
% applies to all data sets.

dataSetNames = {...
%'201612_2048_straightIntlv_wiggleStep_gain_-1000_091216',...
%'201612_2030_straightIntlv_wiggle_gain_-1000_091216',...
'201612_2005_straightIntlv_wiggle_gain_-1000_091216',...
'201612_1920_straightIntlv_wiggle_gain_-1000_091216'...
};

for datSetInd=1:length(dataSetNames)
%% INPUT PARAMETERS

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
loadParametersFromOldFile = false; % If true and processed data already exists in saveDir for this data set, loads the parameters in this section from the old file.
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
parametersFile = '/home/jack/PhaseFeedforward/Analysis/201511/20151120_1538_Gain-800_R56_0.1_Interleaaved_Even.mat'; % file to load parameters from. If empty searches for file from previous processing of current dataset.

% added possibility to process data if diode channels were not
% digitised
diodeIsPresent = true;

% FONT data: If true only do phase analysis, plus different functions for
% loading data etc.
isFONTData = true;
mixerADCs = []; % if empty assumes default [2 4 6]
diodeADCs = []; % If empty assumes default [1 3 5]
sampInterval = []; % sampling rate in ns. If empty assumes 357MHz - entered as 1000/357.
stripChanOffset = false; % if true, adjust PFF on downstream phase to take in to drifts caused by changes made to channel offset

% Interleaved data: if true saves processes odd and even pulses separately.
isInterleaved = true;

% Options for processing bad pulses and drifts
stripNoBeam = true;
stripDuplicates = true;
stripOutliers = true;
stripDrift = false;
driftNAvg = 25; % use average of this many pulses in drift removal

% Phase subtraction: 
% 'none': no phase subtraction
% 'mean': remove mean phase of each monitor
% 'file': use the subtractPhase in the file subtractFile.
% 'value': specify phases to subtract. Below subtractPhase = [mon1,mon2,mon3];
subtractType = 'mean'; 
subtractPhase = [0,0,0]; % only used for value subtraction, otherwise overwritten
subtractFile = '/home/jack/PhaseFeedforward/Analysis/201511/20151123_1921_R5ScanGunWiggle_0.1.mat';

frascatiCalTimeStamp = '20161209_2105';%'20161209_2105';%[]% if this is empty, will take it from data set info.

alignTo = 'end'; % 'start' or 'end'. Whether to align devices to the start or the end of the pulse.
useStaticFrascatiDelay = false; % if true, avoid alignment issues for Frascati monitors by using sample delays below
staticFrascatiDelayType = 'value'; % 'value': use values in staticFrascatiDelay below, 'plot': plot signals and pick range manually, 'calculate': use values returned from calculateDiodeDelay for this data set. 
staticFrascatiDelay = [0 0 0];%[0 0 -163];

sampleRange = 638:738;%638:738;%784:826;%638:738;%565:680[]; % For Fracati monitors. Used to calculate means etc. If this is left empty, will prompt to ask for sample range
pulsesToProcess = [];%[]; % If empty processes all pulses, or define a range of pulses to use.

saveData = true;
saveDir = '/home/jack/PhaseFeedforward/Analysis/201612';

runBPM = false;
refBPM1 = 'CL_SVBPM0502'; % additional correlations etc. calculated for these BPMs (e.g. correlations between BPMs)
refBPM2 = 'CT_SVBPI0608';
transmissionBPM = 'CL_SVBPM0402S';
%bpmSampleRange = 30:40; %50:64;%[];%30:40; % For BPM. Used to calculate means etc. If this is left empty, will prompt to ask for sample range

runBPR = false;
runPETS = false;

dataSetName = dataSetNames{datSetInd};%'20150422_1849_Mix1Mon1PhShft1_Mix2Mon2PhShft2_6dBRemovedPhShft12';

%% load data
addpath('../');
if (isFONTData)
    addpath('../FONT/');
    runBPM = false;
    runBPR = false;
    runPETS = false;
end

if (loadParametersFromOldFile)
    if (isempty(parametersFile))
        saveName = [saveDir '/' dataSetName '.mat'];
    else
        saveName = parametersFile;
    end
    if (exist(saveName,'file'))
        oldData = load(saveName);
        try
            refBPM1 = oldData.refBPM1;
        catch
            refBPM1 = oldData.bpmName;
        end
        try
            frascatiCalTimeStamp = oldData.frascatiCalTimeStamp;
            sampleRange = oldData.sampleRange;
            baselineStartSamples = oldData.baselineStartSamples;
            baselineEndSamples = oldData.baselineEndSamples;
            refBPM2 = oldData.refBPM2;
            isInterleaved = oldData.isInterleaved;
            stripDuplicates = oldData.stripDuplicates;
            stripOutliers = oldData.stripOutliers;
            stripDrift = oldData.stripDrift;
            stripNoBeam = oldData.stripNoBeam;
            driftNAvg = oldData.driftNAvg; % use average of this many pulses in drift removal
            alignTo = oldData.alignTo;
            subtractType = oldData.subtractType; 
            subtractPhase = oldData.subtractPhase; % only used for value subtraction, otherwise overwritten
            runBPM = oldData.runBPM;
            runBPR = oldData.runBPR;
            runPETS = oldData.runPETS;
            subtractFile = oldData.subtractFile;
        catch
        end
        %bpmSampleRange = oldData.bpmSampleRange;
        clear oldData;
    end
end

if (isFONTData)
    [ FONTData, calStruct, dataDir ] = loadFONTData( dataSetName, frascatiCalTimeStamp);
else
    [CTFData, calStruct, dataDir] = loadMergedData(dataSetName,frascatiCalTimeStamp);
end
frascatiCalTimeStamp = calStruct.frascatiCalTimeStamp;
calibrationConstants =  calStruct.calibrationConstants;
useMixerOverSqrtDiode = calStruct.useMixerOverSqrtDiode;

if (~isFONTData)
    dataComment = CTFData(1).comment;
end

if (stripNoBeam)
    if (isFONTData)
        % JACK - to implement. Removal of no beam pulses in FONT data based
        % on no output above noise on diode  (ADC1).
    else
        CTFData = removePulsesNoBeam(CTFData);
    end
end

if (isInterleaved)
    if (isFONTData)
        [FONTDataFFOn,FONTDataFFOff] = splitInterleavedFONTData(FONTData);
        dataToProcess = {FONTDataFFOn FONTDataFFOff};
       
        saveNames = {...
            sprintf('%s_FFOn',dataSetNames{datSetInd}),...
            sprintf('%s_FFOff',dataSetNames{datSetInd})
        };

    else
        CTFDataOdd = CTFData(1:2:end);
        CTFDataEven = CTFData(2:2:end);
        dataToProcess = {CTFDataOdd CTFDataEven};

        saveNames = {...
            sprintf('%s_Odd',dataSetNames{datSetInd}),...
            sprintf('%s_Even',dataSetNames{datSetInd})
        };

    end

else
    if (isFONTData)
        dataToProcess = {FONTData};
    else
        dataToProcess = {CTFData};
    end
    saveNames = {dataSetNames{datSetInd}};
end

for datPartInd = 1:length(dataToProcess) % if interleaved, runs separately for odd and even pulses
    if (isFONTData)
        FONTData = dataToProcess{datPartInd};
        if (~isempty(pulsesToProcess))
            FONTData = extractSubsetFONTData(FONTData,pulsesToProcess);
        end
    else
        CTFData = dataToProcess{datPartInd};
        if (~isempty(pulsesToProcess))
            CTFData = CTFData(pulsesToProcess);
        end
    end
    
    %% 1) Frascati monitors processing/analysis (should always be first - required for other sections to work!)
    fprintf('Processing Frascati phase data...\n');
    processData_Frascati;
    fprintf('Finished processing Frascati phase data!\n');

    %% 2) BPM processing/analysis (has some dependencies on processData_Frascati.m)
    if (runBPM)
        fprintf('Processing BPM data...\n');
        processData_BPM;
        fprintf('Finished processing BPM data!\n');
    end

    %% 3) BPR processing/analysis (to implement)
    if (runBPR)
        fprintf('Processing BPR data...');
        processData_BPR;
        fprintf('done!\n');
    end

    %% 4) PETS processing/analysis (to implement)
    if (runPETS)
        fprintf('Processing PETS data...\n');
        processData_PETS;
        fprintf('done!\n');
    end

    %% Save data

    if (saveData)
        clear CTFData CTFDataOdd CTFDataEven FONTData FONTDataFFOn FONTDataFFOff
        if (~exist(saveDir,'dir'))
            mkdir(saveDir);
        end

        saveName = [saveDir '/' saveNames{datPartInd} '.mat'];
        save(saveName,'-regexp','^(?!(dataToProcess)$).'); % should save all variables except dataToProcess

    end
end

end