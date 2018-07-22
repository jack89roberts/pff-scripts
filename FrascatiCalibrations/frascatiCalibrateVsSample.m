close all; clearvars;
%%

% calTimeStamp = '20150707_1036';
% calTimeStamp = '20150707_1018';

% calTimeStamp = '20150714_1650';
% calTimeStamp = '20150714_1701';
% calTimeStamp = '20150714_1717';

% monCalTimeStamps = {...
%     '20151001_1504',...
%     '20151001_1517',...
%     '20151001_1504'...
% };
% newCalSampleRange = 530:705;%550:714;%260:330; %[]% Samples to calculate calibration for. Leave empty to prompt for sample range to use0

isFONTCal = 0;

% monCalTimeStamps = {...
%     '20151202_1255',...
%     '20151202_1255',...
%     '20151202_1255'...
% };

monCalTimeStamps = {...
    '20151116_1651',...
    '20151116_1651',...
    '20151116_1651'...
};

newCalSampleRange = 500:710;%545:705;%550:714;%260:330; %[]% Samples to calculate calibration for. Leave empty to prompt for sample range to use0
newUseMixerOverSqrtDiode = 0;
newPhaseShiftFreqs = [4 4 4];%[0.67 1 0.74]; %0.866

statsRange = 600:690; % sample range over which to calculate stats on variation in parameters etc.

saveBaseDir = '/home/jack/PhaseFeedforward/Analysis/Calibrations';
savePlots = 0;

monColours = {'b','r','g'};

%% plot raw data (can be used to check files correspond to correct monitors)

if (false)
    fprintf('Plotting raw data...\n');
    for f=1:length(monCalTimeStamps)
        calTimeStamp = monCalTimeStamps{f};
        datTitle = strrep(monCalTimeStamps{f},'_',' ');

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
        clear calibrationFactors zeroCrossings calibrationRSquare;

        if (~isempty(newCalSampleRange))
            calSampleRange = newCalSampleRange;
        end

        for mon=1:nMons
            [alignedMixers, alignedDiodes, calSampleRange] = ...
                frascatiCalProcessSignals(...
                    rawMixers,...
                    rawDiodes,...
                    rawMixerSensitivities,...
                    rawDiodeSensitivities,...
                    calSampleRange,...
                    mon);
            set(gcf,'Position',[225 372 1196 422]);
            set(gcf,'PaperPositionMode','auto');
            subplot(1,2,1)
            title({sprintf('Mon %d (%s)',mon,datTitle),'Diode'});
            subplot(1,2,2)
            hold all;
            title({sprintf('Mon %d (%s)',mon,datTitle),'Mixer'});
        end
    end
    return;
end

%% calculate fits
addpath('../');
saveDir = [saveBaseDir '/' monCalTimeStamps{1}];

nMons = length(monCalTimeStamps);
nSamplesToUse = length(newCalSampleRange);
calFactors = NaN(nMons,nSamplesToUse,4);
zeroCross = NaN(nMons,nSamplesToUse);
calRSquare = NaN(nMons,nSamplesToUse);
fitErr = NaN(nMons,nSamplesToUse,4); % fit parameter confidence intervals

statsRange = statsRange-newCalSampleRange(1);

for mon=1:nMons
    fprintf('-----------MON %d-----------------\n',mon);

    calTimeStamp = monCalTimeStamps{mon};
    
    fprintf('Loading data...\n');
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

        [ rawMixers, rawDiodes, scanPhShiftValues ] = packageFONTCalData(fontCalDir);
        rawMixerSensitivities = '';
        rawDiodeSensitivities = '';
        [nMons,~,~,~] = size(rawMixers);
        savePath = [fontBaseDir '/Processed/'];
        if (~exist(savePath,'dir'))
            mkdir(savePath);
        end    
        phShiftNames = {'Mon1','Mon2','Mon3'};
        monMixerNames = {'Mon1','Mon2','Mon3'};
        
    else
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
    clear calibrationFactors zeroCrossings calibrationRSquare;
    fprintf('Finished loading data.\n');
    
    if (~isempty(newCalSampleRange))
        calSampleRange = newCalSampleRange;
    end
    if (~isempty(newUseMixerOverSqrtDiode))
        useMixerOverSqrtDiode = newUseMixerOverSqrtDiode;
    end
    if (~isempty(newPhaseShiftFreqs))
        phaseShiftFreq = newPhaseShiftFreqs(mon);
    end
    
    % process signals
    fprintf('Processing signals...\n');
    [alignedMixers, alignedDiodes, calSampleRange] = ...
        frascatiCalProcessSignals(...
            rawMixers,...
            rawDiodes,...
            rawMixerSensitivities,...
            rawDiodeSensitivities,...
            calSampleRange,...
            mon);
    fprintf('Finished processing signals!\n');   
    set(gcf,'Position',[225 372 1196 422]);
    set(gcf,'PaperPositionMode','auto');
    datTitle = strrep(monCalTimeStamps{mon},'_',' ');
    subplot(1,2,1)
    title({sprintf('Mon %d (%s)',mon,datTitle),'Diode'});
    subplot(1,2,2)
    hold all;
    title({sprintf('Mon %d (%s)',mon,datTitle),'Mixer'});
    if savePlots; savePlot(saveDir,sprintf('rawMon%d',mon)); end;
        
    % fit data
    fprintf('Fitting data...\n');

    for s=1:nSamplesToUse
        calSamp = calSampleRange(s);
        fprintf('Sample %d...\n',calSamp);

        [tmpCalFactors, tmpZeroCross, tmpRSquare, tmpConf] = ...
            frascatiCalCalculateCalibrations(...
                scanPhShiftValues,...
                alignedMixers,...
                alignedDiodes,...
                calSamp,...
                useMixerOverSqrtDiode,...
                phaseShiftFreq,...
                '',...
                mon);
        
        calFactors(mon,s,:) = tmpCalFactors(mon,:);
        zeroCross(mon,s) = tmpZeroCross(mon);
        calRSquare(mon,s) = tmpRSquare(mon);
        fitErr(mon,s,:) = tmpConf(mon,:)/2;
        
    end
    fprintf('Finished fitting data!\n');

  
end;

%% calculate some stats


%% make plots

fprintf('Plotting fit results...\n');
figAmplitude = figure;
handlesAmplitude = cell(1,nMons);
figOffset = figure;
handlesOffset = cell(1,nMons);
figFreq = figure;
handlesFreq = cell(1,nMons);
figShift = figure;
handlesShift = cell(1,nMons);

for mon=1:nMons
    figure(figAmplitude);
    handlesAmplitude{mon} = shadedErrorBar(calSampleRange,squeeze(calFactors(mon,:,1)),squeeze(fitErr(mon,:,1)),monColours{mon},0.25);
    hold all;
    
    figure(figOffset);
    handlesOffset{mon} = shadedErrorBar(calSampleRange,squeeze(calFactors(mon,:,4)),squeeze(fitErr(mon,:,4)),monColours{mon},0.25);
    hold all;
    
    figure(figFreq);
    handlesFreq{mon} = shadedErrorBar(calSampleRange,squeeze(calFactors(mon,:,2)),squeeze(fitErr(mon,:,2)),monColours{mon},0.25);
    hold all;
    
    figure(figShift);
    handlesShift{mon} = shadedErrorBar(calSampleRange,squeeze(calFactors(mon,:,3)),squeeze(fitErr(mon,:,3)),monColours{mon},0.25);
    hold all;


   
end

datTitle = strrep(monCalTimeStamps{1},'_',' ');

figure(figAmplitude);
xlabel('Sample No.');
ylabel('Calibration Amplitude');
xlim([newCalSampleRange(1) newCalSampleRange(end)]);
legend([handlesAmplitude{1}.mainLine handlesAmplitude{2}.mainLine handlesAmplitude{3}.mainLine],'Mon1','Mon2','Mon3','Location','best');
title(sprintf('Calibration Amplitude (%s)',datTitle));
if savePlots; savePlot(saveDir,'calAmpVsSample'); end;

figure(figOffset);
xlabel('Sample No.');
ylabel('Calibration Offset');
xlim([newCalSampleRange(1) newCalSampleRange(end)]);
title(sprintf('Calibration Offset (%s)',datTitle));
legend([handlesOffset{1}.mainLine handlesOffset{2}.mainLine handlesOffset{3}.mainLine],'Mon1','Mon2','Mon3','Location','best');
if savePlots; savePlot(saveDir,'calOffsetVsSample'); end;

figure(figFreq);
xlabel('Sample No.');
ylabel('Calibration Frequency');
xlim([newCalSampleRange(1) newCalSampleRange(end)]);
title(sprintf('Calibration Frequency (%s)',datTitle));
legend([handlesFreq{1}.mainLine handlesFreq{2}.mainLine handlesFreq{3}.mainLine],'Mon1','Mon2','Mon3','Location','best');
if savePlots; savePlot(saveDir,'calFreqVsSample'); end;

figure(figShift);
xlabel('Sample No.');
ylabel('Calibration Phase Shift');
xlim([newCalSampleRange(1) newCalSampleRange(end)]);
title(sprintf('Calibration Phase Shift (%s)',datTitle));
legend([handlesShift{1}.mainLine handlesShift{2}.mainLine handlesShift{3}.mainLine],'Mon1','Mon2','Mon3','Location','best');
if savePlots; savePlot(saveDir,'calShiftVsSample'); end;


% pick some interesting samples to plot full calibration results:
% samples where fit is the best for each parameter
% samples where fit is the worst for each parameter
% 0.0 0.25, 0.5 and 0.75 1.0 along the sample range
% max fit amplitude
% min fit amplitude
fprintf('Plotting interesting samples...\n');

[~, worstFitsInd] = max(fitErr,[],2);
[~, bestFitsInd] = min(fitErr,[],2);
worstFitsInd = squeeze(worstFitsInd);
bestFitsInd = squeeze(bestFitsInd);
[~, maxFitsInd] = max(calFactors,[],2);
[~, minFitsInd] = min(calFactors,[],2);
maxFitsInd = squeeze(maxFitsInd);
minFitsInd = squeeze(minFitsInd);
sampleQuarters = 1 + ([0 0.25.*nSamplesToUse 0.5.*nSamplesToUse 0.75.*nSamplesToUse (nSamplesToUse-1)]);
sampleQuarters = round(sampleQuarters);

for mon=1:nMons
    fprintf('Mon %d...\n',mon);
    
    samplesToPlot = unique([sampleQuarters bestFitsInd(mon,:) worstFitsInd(mon,:) maxFitsInd(mon,:) minFitsInd(mon,:)]);
    
    calTimeStamp = monCalTimeStamps{mon};
    if (isFONTCal) % load FONT cal data
        addpath('../FONT/');
        fontBaseDir = ['/home/jack/PhaseFeedforward/FONTData/' calTimeStamp(1:6) '/Calibration/' calTimeStamp];
        fontCalDir = [fontBaseDir '/Extracted'];
        [ rawMixers, rawDiodes, scanPhShiftValues ] = packageFONTCalData(fontCalDir);
        rawMixerSensitivities = '';
        rawDiodeSensitivities = '';
        [nMons,~,~,~] = size(rawMixers);
        savePath = [fontBaseDir '/Processed/'];    
    else
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
    clear calibrationFactors zeroCrossings calibrationRSquare;
    if (~isempty(newCalSampleRange))
        calSampleRange = newCalSampleRange;
    end
    if (~isempty(newUseMixerOverSqrtDiode))
        useMixerOverSqrtDiode = newUseMixerOverSqrtDiode;
    end
    if (~isempty(newPhaseShiftFreqs))
        phaseShiftFreq = newPhaseShiftFreqs(mon);
    end
    
    [alignedMixers, alignedDiodes, calSampleRange] = ...
    frascatiCalProcessSignals(...
        rawMixers,...
        rawDiodes,...
        rawMixerSensitivities,...
        rawDiodeSensitivities,...
        calSampleRange,...
        mon);
    close(gcf)
    
    for s=samplesToPlot
        calsToPlot = squeeze(calFactors(:,s,:));
        errsToPlot = squeeze(fitErr(:,s,:));
        sampleToPlot = calSampleRange(s);
        plotFitResult(mon,scanPhShiftValues,alignedMixers,alignedDiodes,useMixerOverSqrtDiode,sampleToPlot,calsToPlot,phaseShiftFreq);
        titleLine1 = sprintf('Mon%d, Sample %d (%s)',mon,sampleToPlot,strrep(calTimeStamp,'_',' '));
        titleLine2 = sprintf('Fit params: a=%.2f, b=%.2f, c=%.3f, d=%.3f',calsToPlot(mon,1),calsToPlot(mon,2),calsToPlot(mon,3),calsToPlot(mon,4));
        titleLine3 = sprintf('Fit conf int: a:%.2f, b:%.2f, c:%.3f, d:%.3f',errsToPlot(mon,1),errsToPlot(mon,2),errsToPlot(mon,3),errsToPlot(mon,4));
        title({titleLine1, titleLine2, titleLine3});
        if savePlots; savePlot(saveDir,sprintf('Mon%d_Samp%d',mon,sampleToPlot)); end;
        close(gcf);
    end
    
end

%%
return;
sampleToPlot = 634;
sampIndex = find(calSampleRange==sampleToPlot);
calsToPlot = squeeze(calFactors(:,sampIndex,:));
plotFitResult(scanPhShiftValues,alignedMixers,alignedDiodes,useMixerOverSqrtDiode,sampleToPlot,calsToPlot);