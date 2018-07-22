close all;
%%
dataSetNames = {...
    '20141210_2035_R56_+0.3_Gate_375_435_GainK1_0_GainK2_0',...
    '20141210_2035_R56_+0.3_Gate_375_435_GainK1_20_GainK2_-20',...
    '20141210_2035_R56_+0.3_Gate_375_435_GainK1_40_GainK2_-40',...
    '20141210_2035_R56_+0.3_Gate_375_435_GainK1_63_GainK2_-63_END',...
    '20141210_2035_R56_+0.3_Gate_375_435_GainK1_-20_GainK2_20',...
    '20141210_2035_R56_+0.3_Gate_375_435_GainK1_-40_GainK2_40',...
    '20141210_2035_R56_+0.3_Gate_375_435_GainK1_-63_GainK2_63'...
    };

% dataSetNames = {...
%     '20141210_1915_R56_+0.3_Gate_350_410_ConstK1_0_ConstK2_0',...
%     '20141210_1915_R56_+0.3_Gate_375_435_ConstK1_800_ConstK2_-800',...
%     '20141210_1915_R56_+0.3_Gate_375_435_ConstK1_-800_ConstK2_800'...
% };


dataSetLabels = {...
    'Gain = 0',...
    'Gain = +20',...
    'Gain = +40',...
    'Gain = +63',...
    'Gain = -20',...
    'Gain = -40',...
    'Gain = -63',...
    };

dataSetColours = {'k','b','r','g','b','r','g'};
dataSetStyles = {'-','-','-','-','--','--','--'};

removeThreshold = 0.01; % remove pulses that are removeThreshold % away from max on diode
skipPercent = 20;
legendString = {'Mon1','Mon2','Mon3','PETS'};
monColours = {'b','r','g','m'};
%mon12StartSample = 150;
%mon12EndSample = 240;
mon12StartSample = 126;
mon12EndSample = 163;
mon3Delay = 648;
petsDelay = -742; % delay between PETS and Mon3 (will change depending on oasis setup) in ns.

plotData = false;

%%
fprintf(1,'Loading processed data files...\n');
processedData = loadProcessedData(dataSetNames);
nDataSets = length(dataSetNames);

% create arrays to contain everything
diodes = cell(1,nDataSets);
mixers = cell(1,nDataSets);
phases =  cell(1,nDataSets);
pets_Phase = cell(1,nDataSets);
pets_I = cell(1,nDataSets);
pets_Q = cell(1,nDataSets);
nMonitors  = cell(1,nDataSets);
nPulses = cell(1,nDataSets);
nSamples = cell(1,nDataSets);
petsNSamples = cell(1,nDataSets);
frascatiTimePerSample = cell(1,nDataSets);
petsTimePerSample = cell(1,nDataSets);
startSamples = cell(1,nDataSets);
endSamples = cell(1,nDataSets);
frascatiTimeAxes = cell(1,nDataSets);
petsStartSample = cell(1,nDataSets);
petsEndSample = cell(1,nDataSets);
petsTimeAxis = cell(1,nDataSets);
meanDiodes = cell(1,nDataSets);
meanMixers = cell(1,nDataSets);
meanPhases = cell(1,nDataSets);% mean per sample
stdPhases = cell(1,nDataSets);% std per sample
meanStdPhase = cell(1,nDataSets);
diffPhases12 = cell(1,nDataSets);
diffPhases13 = cell(1,nDataSets);
diffPhases23 = cell(1,nDataSets);
resolution12 = cell(1,nDataSets);
resolution13 = cell(1,nDataSets);
resolution23 = cell(1,nDataSets);
meanResolution12 = cell(1,nDataSets);
meanResolution13 = cell(1,nDataSets);
meanResolution23 = cell(1,nDataSets);
meanPulseDiodes = cell(1,nDataSets);
meanPulseMixers = cell(1,nDataSets);
meanPulsePhases = cell(1,nDataSets); % mean per pulse
stdPulseDiodes = cell(1,nDataSets);
stdPulseMixers = cell(1,nDataSets);
stdPulsePhases = cell(1,nDataSets); % error bars on mean pulse phase
corr12_meanPulsePhase = cell(1,nDataSets);
corr13_meanPulsePhase = cell(1,nDataSets);
corr23_meanPulsePhase = cell(1,nDataSets);
fit12_meanPulsePhase = cell(1,nDataSets);
fit13_meanPulsePhase = cell(1,nDataSets);
fit23_meanPulsePhase = cell(1,nDataSets);
meanDiffPhase12 = cell(1,nDataSets);
meanDiffPhase13 = cell(1,nDataSets);
meanDiffPhase23 = cell(1,nDataSets);
phases_diffToMean = cell(1,nDataSets);
corr12_samples = cell(1,nDataSets);
corr13_samples = cell(1,nDataSets);
corr23_samples = cell(1,nDataSets);
corr12_diffToMean = cell(1,nDataSets);
corr13_diffToMean = cell(1,nDataSets);
corr23_diffToMean = cell(1,nDataSets);
fit12_samples = cell(1,nDataSets);
fit13_samples = cell(1,nDataSets);
fit23_samples = cell(1,nDataSets);
fit12_diffToMean = cell(1,nDataSets);
fit13_diffToMean = cell(1,nDataSets);
fit23_diffToMean = cell(1,nDataSets);
pets_meanPhases = cell(1,nDataSets); % mean per sample
pets_meanPulsePhases = cell(1,nDataSets); % mean per pulse
pets_stdPhases = cell(1,nDataSets); % std per sample
pets_meanStdPhase = cell(1,nDataSets); % mean sample std
pets_stdPulsePhase = cell(1,nDataSets);
corr_PetsFrascati3_meanPulsePhase = cell(1,nDataSets);
fit_PetsFrascati3_meanPulsePhase = cell(1,nDataSets);
corr_PetsFrascati1_meanPulsePhase = cell(1,nDataSets);
fit_PetsFrascati1_meanPulsePhase = cell(1,nDataSets);
corr_PetsFrascati2_meanPulsePhase = cell(1,nDataSets);
fit_PetsFrascati2_meanPulsePhase = cell(1,nDataSets);
skipPhase = cell(1,nDataSets);
skipBPM = cell(1,nDataSets);
bpmCC930_H = cell(1,nDataSets);
meanBPMCC930_H = cell(1,nDataSets);

for dataSetIndex = 1:nDataSets
    fprintf('Analysing %s...\n',dataSetNames{dataSetIndex});
    tmpProcessedData = processedData{dataSetIndex};
    dataPath = tmpProcessedData.dataDir;
    
    saveDir = [dataPath '/plots/'];
    if (~exist(saveDir,'dir'))
        mkdir(saveDir);
    end
    
    tmpSkipPhase = 0;
    tmpSkipBPM = 0;

    try
        tmpDiodes = tmpProcessedData.phaseData.frascatiDiodes;
        tmpMixers = tmpProcessedData.phaseData.frascatiMixers;
        tmpPhases =  tmpProcessedData.phaseData.frascatiPhases;
        tmpPets_Phase = tmpProcessedData.phaseData.petsPhase;
        tmpPets_I = tmpProcessedData.phaseData.petsI;
        tmpPets_Q = tmpProcessedData.phaseData.petsQ;
        tmpFrascatiTimePerSample = tmpProcessedData.phaseData.frascatiTimePerSample;
        tmpPetsTimePerSample = tmpProcessedData.phaseData.petsTimePerSample;
    catch
        tmpSkipPhase = 1;
        warning('Error accessing phase data, skipping it');
    end
    
    [tmpNMonitors,tmpNPulses,tmpNSamples] = size(tmpMixers);
    [~,tmpPetsNSamples] = size(tmpPets_I);
    
    try
        tmpBPMCC930_H = tmpProcessedData.bpmData.CC.SVBPM0930H.Samples.samples;
    catch
        tmpSkipBPM = 1;
        warning('Error accessing BPM data, skipping it')
    end
    
    % % to help find mon3 delay
    % crosscorr(squeeze(meanPhases(1,:)),squeeze(meanPhases(3,:)),400)
    % 
    % % to help find pets Delay
    % tmpFrascati = meanPhases(3,:);%-mean(meanPhases(3,startSamples(3):endSamples(3)));
    % tmpPETS = pets_meanPhases; %- mean(pets_meanPhases(petsStartSample:petsEndSample));
    modSamples = linspace(1,tmpNSamples,tmpNSamples.*(tmpFrascatiTimePerSample./tmpPetsTimePerSample));
    % yy=spline(1:nSamples,tmpFrascati,linspace(1,nSamples,nSamples.*(frascatiTimePerSample./petsTimePerSample)));
    % crosscorr(yy,tmpPETS,900)
    % [~,maxInd] = max(crosscorr(yy,tmpPETS,900));
    % delaySamples = -900:900;
    % delaySamples = delaySamples(maxInd);
    % 
    % figure;
    % plot(4*(1:nSamples),meanPhases(3,:));%-mean(meanPhases(3,startSamples(3):endSamples(3))))
    % hold all;
    % plot((1:petsNSamples)-delaySamples,pets_meanPhases);%-mean(pets_meanPhases(petsStartSample:petsEndSample)))
    % plot((1:petsNSamples),pets_meanPhases);%-mean(pets_meanPhases(petsStartSample:petsEndSample)))

    if (~tmpSkipPhase)
    
        transmission = NaN*ones(tmpNMonitors,tmpNPulses);
        isGoodPulseMon = false(tmpNMonitors,tmpNPulses);

        for mon = 1:tmpNMonitors
            [startSamp,endSamp] = getSignalRange(squeeze(nanmean(tmpDiodes(mon,:,:))),20);
            transmission(mon,:) = squeeze(nanmean(tmpDiodes(mon,:,startSamp:endSamp),3));
            isGoodPulseMon(mon,:) = (transmission(mon,:)-min(transmission(mon,:))) < removeThreshold;
        end

        isGoodPulseAll = ones(1,tmpNPulses);
        for mon = 1:tmpNMonitors
            isGoodPulseAll = isGoodPulseAll & isGoodPulseMon(mon,:);
        end

        tmpMixers = tmpMixers(:,isGoodPulseAll,:);
        tmpDiodes = tmpDiodes(:,isGoodPulseAll,:);
        tmpPhases = tmpPhases(:,isGoodPulseAll,:);
        tmpPets_I = tmpPets_I(isGoodPulseAll,:);
        tmpPets_Q = tmpPets_Q(isGoodPulseAll,:);
        tmpPets_Phase = tmpPets_Phase(isGoodPulseAll,:);

        [~,tmpNPulses,~] = size(tmpMixers);
        
        if (~tmpSkipBPM)
            tmpBPMCC930_H = tmpBPMCC930_H(isGoodPulseAll,:);
        end
        
        pulseWidth = mon12EndSample-mon12StartSample;
        mon3StartSample = mon12StartSample + round(mon3Delay./tmpFrascatiTimePerSample);
        mon3EndSample = mon3StartSample + pulseWidth;
        if (mon3EndSample > tmpNSamples)
            mon3EndSample = tmpNSamples;
            warning('Mon3 end sample set to last sample as out of bounds index calculated');
        end
        tmpStartSamples = [mon12StartSample mon12StartSample mon3StartSample];
        tmpEndSamples = [mon12EndSample mon12EndSample mon3EndSample];
        tmpFrascatiTimeAxes = NaN*ones(tmpNMonitors,tmpNSamples);
        for mon=1:tmpNMonitors
            tmpFrascatiTimeAxes(mon,:) = ((1:tmpNSamples)-tmpStartSamples(mon)).*tmpFrascatiTimePerSample;
        end

        [~,mod3StartInd] = min(abs(modSamples-tmpStartSamples(3)));
        tmpPetsStartSample = mod3StartInd+petsDelay;%mon3StartSample - round(petsDelay./petsTimePerSample);
        if (tmpPetsStartSample<1)
            tmpPetsStartSample = 1;
            warning('PETS start sample set to 1 as negative value calculated');
        end
        tmpPetsEndSample = tmpPetsStartSample + round(pulseWidth.*(tmpFrascatiTimePerSample./tmpPetsTimePerSample));
        if (tmpPetsEndSample > tmpPetsNSamples)
            tmpPetsEndSample = tmpPetsNSamples;
            warning('PETS start sample set to 1 as negative value calculated');
        end
        tmpPetsTimeAxis = ((1:tmpPetsNSamples)-tmpPetsStartSample).*tmpPetsTimePerSample;
    
    end


    % pulseStartSamples = NaN*ones(nMonitors,nPulses);
    % pulseEndSamples = NaN*ones(nMonitors,nPulses);
    % for mon=1:nMonitors
    %     for p=1:nPulses
    %         [pulseStartSamples(mon,p),pulseEndSamples(mon,p)] = getSignalRange(squeeze(diodes(mon,p,:)),skipPercent);
    %     end
    % end
    % startSamples = round(nanmean(pulseStartSamples,2));
    % endSamples = round(nanmean(pulseEndSamples,2));
    % frascatiTimeAxes = NaN*ones(nMonitors,nSamples);
    % for mon=1:nMonitors
    %     frascatiTimeAxes(mon,:) = ((1:nSamples)-startSamples(mon)).*frascatiTimePerSample;
    % end

    % 
    % petsPulseStartSamples = NaN*ones(nPulses,1);
    % petsPulseEndSamples = NaN*ones(nPulses,1);
    % for p=1:nPulses
    %     [iStart,iEnd] = getSignalRange(squeeze(pets_I(p,:)),skipPercent);
    %     %[qStart,qEnd] = getSignalRange(squeeze(pets_Q(p,:)),skipPercent);
    %     %petsPulseStartSamples(p) = nanmean([iStart qStart]);
    %     %petsPulseEndSamples(p) = nanmean([iEnd qEnd]);
    %     petsPulseStartSamples(p) = iStart;
    %     petsPulseEndSamples(p) = iEnd;
    % end
    % petsStartSample = round(nanmean(petsPulseStartSamples));
    % petsEndSample = round(nanmean(petsPulseEndSamples));
    % 
    % petsTimeAxis = ((1:petsNSamples)-petsStartSample).*petsTimePerSample;

    % 
    % figure;
    % plotHandles = NaN*ones(1,4);
    % for mon=1:nMonitors
    %     plotHandles(mon) = plot(pulseStartSamples(mon,:),plotColours{mon});
    %     hold all;
    %     plot(pulseEndSamples(mon,:),plotColours{mon});
    % end
    % plotHandles(4) = plot(petsPulseStartSamples,plotColours{4});
    % plot(petsPulseEndSamples,plotColours{4});
    % title('Calculated Pulse Start and End Points');
    % xlabel('Pulse No.');
    % ylabel('Start/End Sample');
    % legend(plotHandles,legendString);

    %% Frascati
    if (~tmpSkipPhase)
    
        tmpMeanDiodes = squeeze(nanmean(tmpDiodes,2));
        tmpMeanMixers = squeeze(nanmean(tmpMixers,2));
        tmpMeanPhases = squeeze(nanmean(tmpPhases,2)); % mean per sample
        tmpStdPhases = squeeze(nanstd(tmpPhases,0,2)); % std per sample

        
        % outlier removal
        %for 
        
        
        for mon=1:tmpNMonitors
            tmpSubtractMean = nanmean(tmpMeanPhases(mon,tmpStartSamples(mon):tmpEndSamples(mon)));
            tmpMeanPhases(mon,:) = tmpMeanPhases(mon,:); %- tmpSubtractMean;

            for p=1:tmpNPulses
                tmpPhases(mon,p,:) = tmpPhases(mon,p,:);%-tmpSubtractMean;
            end
        end

        tmpMeanStdPhase = NaN*ones(1,tmpNMonitors);
        for mon=1:tmpNMonitors
            tmpMeanStdPhase(mon) = squeeze(nanmean(tmpStdPhases(mon,tmpStartSamples(mon):tmpEndSamples(mon)),2)); % mean sample std
        end
        % 
        tmpDiffPhases12 = tmpPhases(1,:,tmpStartSamples(1):tmpEndSamples(1)) - tmpPhases(2,:,tmpStartSamples(2):tmpEndSamples(2));
        tmpDiffPhases13 = tmpPhases(1,:,tmpStartSamples(1):tmpEndSamples(1)) - tmpPhases(3,:,tmpStartSamples(3):tmpEndSamples(3));
        tmpDiffPhases23 = tmpPhases(2,:,tmpStartSamples(2):tmpEndSamples(2)) - tmpPhases(3,:,tmpStartSamples(3):tmpEndSamples(3));
        % 
        tmpResolution12 = squeeze(nanstd(tmpDiffPhases12,0,3))./sqrt(2);
        tmpResolution13 = squeeze(nanstd(tmpDiffPhases13,0,3))./sqrt(2);
        tmpResolution23 = squeeze(nanstd(tmpDiffPhases23,0,3))./sqrt(2);
        tmpMeanResolution12 = mean(tmpResolution12);
        tmpMeanResolution13 = mean(tmpResolution13);
        tmpMeanResolution23 = mean(tmpResolution23);
        % 
        tmpMeanPulseDiodes = NaN*ones(tmpNMonitors,tmpNPulses);
        tmpMeanPulseMixers = NaN*ones(tmpNMonitors,tmpNPulses);
        tmpMeanPulsePhases = NaN*ones(tmpNMonitors,tmpNPulses); % mean per pulse
        for mon = 1:tmpNMonitors
            tmpMeanPulseDiodes(mon,:) = squeeze(nanmean(tmpDiodes(mon,:,tmpStartSamples(mon):tmpEndSamples(mon)),3));
            tmpMeanPulseMixers(mon,:) = squeeze(nanmean(tmpMixers(mon,:,tmpStartSamples(mon):tmpEndSamples(mon)),3));
            tmpMeanPulsePhases(mon,:) = squeeze(nanmean(tmpPhases(mon,:,tmpStartSamples(mon):tmpEndSamples(mon)),3));
        end
        tmpStdPulseDiodes =std(tmpMeanPulseDiodes,0,2); 
        tmpStdPulseMixers =std(tmpMeanPulseMixers,0,2);
        tmpStdPulsePhases =std(tmpMeanPulsePhases,0,2); % error bars on mean pulse phase
        % 
        tmpCorr12_meanPulsePhase = corrcoef(tmpMeanPulsePhases(1,:),tmpMeanPulsePhases(2,:));
        tmpCorr13_meanPulsePhase = corrcoef(tmpMeanPulsePhases(1,:),tmpMeanPulsePhases(3,:));
        tmpCorr23_meanPulsePhase = corrcoef(tmpMeanPulsePhases(2,:),tmpMeanPulsePhases(3,:));
        tmpCorr12_meanPulsePhase = tmpCorr12_meanPulsePhase(1,2);
        tmpCorr13_meanPulsePhase = tmpCorr13_meanPulsePhase(1,2);
        tmpCorr23_meanPulsePhase = tmpCorr23_meanPulsePhase(1,2);

        tmpFit12_meanPulsePhase = polyfit(tmpMeanPulsePhases(1,:),tmpMeanPulsePhases(2,:),1);
        tmpFit13_meanPulsePhase = polyfit(tmpMeanPulsePhases(1,:),tmpMeanPulsePhases(3,:),1);
        tmpFit23_meanPulsePhase = polyfit(tmpMeanPulsePhases(2,:),tmpMeanPulsePhases(3,:),1);
        tmpFit12_meanPulsePhase = tmpFit12_meanPulsePhase(1);
        tmpFit13_meanPulsePhase = tmpFit13_meanPulsePhase(1);
        tmpFit23_meanPulsePhase = tmpFit23_meanPulsePhase(1);
        % 
        tmpMeanDiffPhase12 = tmpMeanPhases(1,tmpStartSamples(1):tmpEndSamples(1)) - tmpMeanPhases(2,tmpStartSamples(2):tmpEndSamples(2));
        tmpMeanDiffPhase13 = tmpMeanPhases(1,tmpStartSamples(1):tmpEndSamples(1)) - tmpMeanPhases(3,tmpStartSamples(3):tmpEndSamples(3));
        tmpMeanDiffPhase23 = tmpMeanPhases(2,tmpStartSamples(2):tmpEndSamples(2)) - tmpMeanPhases(3,tmpStartSamples(3):tmpEndSamples(3));
        % 

        tmpPhases_diffToMean = NaN*ones(tmpNMonitors,tmpNPulses,tmpNSamples);
        for mon=1:tmpNMonitors
            for p=1:tmpNPulses
                tmpMeanPhase = nanmean(squeeze(tmpPhases(mon,p,tmpStartSamples(mon):tmpEndSamples(mon))));
                tmpPhases_diffToMean(mon,p,:) = tmpPhases(mon,p,:) - tmpMeanPhase;
            end
        end

        tmpCorr12_samples = NaN*ones(1,pulseWidth+1);
        tmpCorr13_samples = NaN*ones(1,pulseWidth+1);
        tmpCorr23_samples = NaN*ones(1,pulseWidth+1);
        tmpCorr12_diffToMean = NaN*ones(1,pulseWidth+1);
        tmpCorr13_diffToMean = NaN*ones(1,pulseWidth+1);
        tmpCorr23_diffToMean = NaN*ones(1,pulseWidth+1);
        tmpFit12_samples = NaN*ones(1,pulseWidth+1);
        tmpFit13_samples = NaN*ones(1,pulseWidth+1);
        tmpFit23_samples = NaN*ones(1,pulseWidth+1);
        tmpFit12_diffToMean = NaN*ones(1,pulseWidth+1);
        tmpFit13_diffToMean = NaN*ones(1,pulseWidth+1);
        tmpFit23_diffToMean = NaN*ones(1,pulseWidth+1);
        for s=1:pulseWidth+1
            sMon1 = tmpStartSamples(1)+s-1;
            sMon2 = tmpStartSamples(2)+s-1;
            sMon3 = tmpStartSamples(3)+s-1;

            tmp12 = corrcoef(squeeze(tmpPhases(1,:,sMon1)),squeeze(tmpPhases(2,:,sMon2)));
            tmp13 = corrcoef(squeeze(tmpPhases(1,:,sMon1)),squeeze(tmpPhases(3,:,sMon3)));
            tmp23 = corrcoef(squeeze(tmpPhases(2,:,sMon2)),squeeze(tmpPhases(3,:,sMon3)));
            tmpCorr12_samples(s) = tmp12(1,2);
            tmpCorr13_samples(s) = tmp13(1,2);
            tmpCorr23_samples(s) = tmp23(1,2);

            tmp12 = polyfit(squeeze(tmpPhases(1,:,sMon1)),squeeze(tmpPhases(2,:,sMon2)),1);
            tmp13 = polyfit(squeeze(tmpPhases(1,:,sMon1)),squeeze(tmpPhases(3,:,sMon3)),1);
            tmp23 = polyfit(squeeze(tmpPhases(2,:,sMon2)),squeeze(tmpPhases(3,:,sMon3)),1);
            tmpFit12_samples(s) = tmp12(1);
            tmpFit13_samples(s) = tmp13(1);
            tmpFit23_samples(s) = tmp23(1);

            tmp12 = corrcoef(squeeze(tmpPhases_diffToMean(1,:,sMon1)),squeeze(tmpPhases_diffToMean(2,:,sMon2)));
            tmp13 = corrcoef(squeeze(tmpPhases_diffToMean(1,:,sMon1)),squeeze(tmpPhases_diffToMean(3,:,sMon3)));
            tmp23 = corrcoef(squeeze(tmpPhases_diffToMean(2,:,sMon2)),squeeze(tmpPhases_diffToMean(3,:,sMon3)));
            tmpCorr12_diffToMean(s) = tmp12(1,2);
            tmpCorr13_diffToMean(s) = tmp13(1,2);
            tmpCorr23_diffToMean(s) = tmp23(1,2);

            tmp12 = polyfit(squeeze(tmpPhases_diffToMean(1,:,sMon1)),squeeze(tmpPhases_diffToMean(2,:,sMon2)),1);
            tmp13 = polyfit(squeeze(tmpPhases_diffToMean(1,:,sMon1)),squeeze(tmpPhases_diffToMean(3,:,sMon3)),1);
            tmp23 = polyfit(squeeze(tmpPhases_diffToMean(2,:,sMon2)),squeeze(tmpPhases_diffToMean(3,:,sMon3)),1);
            tmpFit12_diffToMean(s) = tmp12(1);
            tmpFit13_diffToMean(s) = tmp13(1);
            tmpFit23_diffToMean(s) = tmp23(1);

        end

    end
    %% PETS
    if (~tmpSkipPhase)
        tmpPets_meanPhases = nanmean(tmpPets_Phase); % mean per sample
        tmpSubtractMean = nanmean(tmpPets_meanPhases(tmpPetsStartSample:tmpPetsEndSample));
        tmpPets_meanPhases = tmpPets_meanPhases; %- tmpSubtractMean;
        for p=1:tmpNPulses
            tmpPets_Phase(p,:) = tmpPets_Phase(p,:);% - tmpSubtractMean;
        end
        tmpPets_meanPulsePhases = nanmean(tmpPets_Phase(:,tmpPetsStartSample:tmpPetsEndSample),2); % mean per pulse
        tmpPets_stdPhases = nanstd(tmpPets_Phase); % std per sample
        tmpPets_meanStdPhase = nanmean(tmpPets_stdPhases(tmpPetsStartSample:tmpPetsEndSample)); % mean sample std
        tmpPets_stdPulsePhase = std(tmpPets_meanPulsePhases); 

        tmpCorr_PetsFrascati3_meanPulsePhase = corrcoef(tmpMeanPulsePhases(3,:),tmpPets_meanPulsePhases);
        tmpCorr_PetsFrascati3_meanPulsePhase = tmpCorr_PetsFrascati3_meanPulsePhase(1,2);
        tmpFit_PetsFrascati3_meanPulsePhase = polyfit(tmpMeanPulsePhases(3,:)',tmpPets_meanPulsePhases,1);
        tmpFit_PetsFrascati3_meanPulsePhase = tmpFit_PetsFrascati3_meanPulsePhase(1);

        tmpCorr_PetsFrascati1_meanPulsePhase = corrcoef(tmpMeanPulsePhases(1,:),tmpPets_meanPulsePhases);
        tmpCorr_PetsFrascati1_meanPulsePhase = tmpCorr_PetsFrascati1_meanPulsePhase(1,2);
        tmpFit_PetsFrascati1_meanPulsePhase = polyfit(tmpMeanPulsePhases(1,:)',tmpPets_meanPulsePhases,1);
        tmpFit_PetsFrascati1_meanPulsePhase = tmpFit_PetsFrascati1_meanPulsePhase(1);

        tmpCorr_PetsFrascati2_meanPulsePhase = corrcoef(tmpMeanPulsePhases(2,:),tmpPets_meanPulsePhases);
        tmpCorr_PetsFrascati2_meanPulsePhase = tmpCorr_PetsFrascati2_meanPulsePhase(1,2);
        tmpFit_PetsFrascati2_meanPulsePhase = polyfit(tmpMeanPulsePhases(2,:)',tmpPets_meanPulsePhases,1);
        tmpFit_PetsFrascati2_meanPulsePhase = tmpFit_PetsFrascati2_meanPulsePhase(1);
    
    end
    %% BPM
    if (~tmpSkipBPM)
        tmpMeanBPMCC930_H = nanmean(tmpBPMCC930_H);
    end
    
    %% Save new data to arrays
    skipPhase{dataSetIndex} = tmpSkipPhase;
    skipBPM{dataSetIndex} = tmpSkipBPM;
    
    if (~tmpSkipPhase)
        diodes{dataSetIndex} = tmpDiodes;
        mixers{dataSetIndex} = tmpMixers;
        phases{dataSetIndex} = tmpPhases;
        pets_Phase{dataSetIndex} = tmpPets_Phase;
        pets_I{dataSetIndex} = tmpPets_I;
        pets_Q{dataSetIndex} = tmpPets_Q;
        nMonitors{dataSetIndex} = tmpNMonitors;
        nPulses{dataSetIndex} = tmpNPulses;
        nSamples{dataSetIndex} = tmpNSamples;
        petsNSamples{dataSetIndex} = tmpPetsNSamples;
        frascatiTimePerSample{dataSetIndex} = tmpFrascatiTimePerSample;
        petsTimePerSample{dataSetIndex} = tmpPetsTimePerSample;
        startSamples{dataSetIndex} = tmpStartSamples;
        endSamples{dataSetIndex} = tmpEndSamples;
        frascatiTimeAxes{dataSetIndex} = tmpFrascatiTimeAxes;
        petsStartSample{dataSetIndex} = tmpPetsStartSample;
        petsEndSample{dataSetIndex} = tmpPetsEndSample;
        petsTimeAxis{dataSetIndex} = tmpPetsTimeAxis;
        meanDiodes{dataSetIndex} = tmpMeanDiodes;
        meanMixers{dataSetIndex} = tmpMeanMixers;
        meanPhases{dataSetIndex} = tmpMeanPhases;
        stdPhases{dataSetIndex} = tmpStdPhases;
        meanStdPhase{dataSetIndex} = tmpMeanStdPhase;
        diffPhases12{dataSetIndex} = tmpDiffPhases12;
        diffPhases13{dataSetIndex} = tmpDiffPhases13;
        diffPhases23{dataSetIndex} = tmpDiffPhases23;
        resolution12{dataSetIndex} = tmpResolution12;
        resolution13{dataSetIndex} = tmpResolution13;
        resolution23{dataSetIndex} = tmpResolution23;
        meanResolution12{dataSetIndex} = tmpMeanResolution12;
        meanResolution13{dataSetIndex} = tmpMeanResolution13;
        meanResolution23{dataSetIndex} = tmpMeanResolution23;
        meanPulseDiodes{dataSetIndex} = tmpMeanPulseDiodes;
        meanPulseMixers{dataSetIndex} = tmpMeanPulseMixers;
        meanPulsePhases{dataSetIndex} = tmpMeanPulsePhases;
        stdPulseDiodes{dataSetIndex} = tmpStdPulseDiodes;
        stdPulseMixers{dataSetIndex} = tmpStdPulseMixers;
        stdPulsePhases{dataSetIndex} = tmpStdPulsePhases;
        corr12_meanPulsePhase{dataSetIndex} = tmpCorr12_meanPulsePhase;
        corr13_meanPulsePhase{dataSetIndex} = tmpCorr13_meanPulsePhase;
        corr23_meanPulsePhase{dataSetIndex} = tmpCorr23_meanPulsePhase;
        fit12_meanPulsePhase{dataSetIndex} = tmpFit12_meanPulsePhase;
        fit13_meanPulsePhase{dataSetIndex} = tmpFit13_meanPulsePhase;
        fit23_meanPulsePhase{dataSetIndex} = tmpFit23_meanPulsePhase;
        meanDiffPhase12{dataSetIndex} = tmpMeanDiffPhase12;
        meanDiffPhase13{dataSetIndex} = tmpMeanDiffPhase13;
        meanDiffPhase23{dataSetIndex} = tmpMeanDiffPhase23;
        phases_diffToMean{dataSetIndex} = tmpPhases_diffToMean;
        corr12_samples{dataSetIndex} = tmpCorr12_samples;
        corr13_samples{dataSetIndex} = tmpCorr13_samples;
        corr23_samples{dataSetIndex} = tmpCorr23_samples;
        corr12_diffToMean{dataSetIndex} = tmpCorr12_diffToMean;
        corr13_diffToMean{dataSetIndex} = tmpCorr13_diffToMean;
        corr23_diffToMean{dataSetIndex} = tmpCorr23_diffToMean;
        fit12_samples{dataSetIndex} = tmpFit12_samples;
        fit13_samples{dataSetIndex} = tmpFit13_samples;
        fit23_samples{dataSetIndex} = tmpFit23_samples;
        fit12_diffToMean{dataSetIndex} = tmpFit12_diffToMean;
        fit13_diffToMean{dataSetIndex} = tmpFit13_diffToMean;
        fit23_diffToMean{dataSetIndex} = tmpFit23_diffToMean;
        pets_meanPhases{dataSetIndex} = tmpPets_meanPhases;
        pets_meanPulsePhases{dataSetIndex} = tmpPets_meanPulsePhases;
        pets_stdPhases{dataSetIndex} = tmpPets_stdPhases;
        pets_meanStdPhase{dataSetIndex} = tmpPets_meanStdPhase;
        pets_stdPulsePhase{dataSetIndex} = tmpPets_stdPulsePhase;
        corr_PetsFrascati3_meanPulsePhase{dataSetIndex} = tmpCorr_PetsFrascati3_meanPulsePhase;
        fit_PetsFrascati3_meanPulsePhase{dataSetIndex} = tmpFit_PetsFrascati3_meanPulsePhase;
        corr_PetsFrascati1_meanPulsePhase{dataSetIndex} = tmpCorr_PetsFrascati1_meanPulsePhase;
        fit_PetsFrascati1_meanPulsePhase{dataSetIndex} = tmpFit_PetsFrascati1_meanPulsePhase;
        corr_PetsFrascati2_meanPulsePhase{dataSetIndex} = tmpCorr_PetsFrascati2_meanPulsePhase;
        fit_PetsFrascati2_meanPulsePhase{dataSetIndex} = tmpFit_PetsFrascati2_meanPulsePhase;
    end
    if (~tmpSkipBPM)
        bpmCC930_H{dataSetIndex} = tmpBPMCC930_H;
        meanBPMCC930_H{dataSetIndex} = tmpMeanBPMCC930_H;
    end
    %% Plots
    if (~plotData)
        continue;
    end
    
    if (~tmpSkipPhase)
        figure;
        plotHandles = NaN*ones(1,4);
        for mon=1:tmpNMonitors
            plotHandles(mon) = plot(tmpFrascatiTimeAxes(mon,tmpStartSamples(mon):tmpEndSamples(mon)),squeeze(tmpMeanPhases(mon,tmpStartSamples(mon):tmpEndSamples(mon))),monColours{mon});
            hold all;
        end
        plotHandles(4) = plot(tmpPetsTimeAxis(tmpPetsStartSample:tmpPetsEndSample),tmpPets_meanPhases(tmpPetsStartSample:tmpPetsEndSample),monColours{4});
        legend(plotHandles,legendString);
        xlabel('Time [ns]');
        ylabel('Phase [degrees]');
        title('Mean Phase vs. Sample No.');
        print([saveDir 'MeanPhaseVsSampleNo.png'],'-dpng');
        hold off;
        
        %figure;
        plotHandles = NaN*ones(1,4);
        legStrs = cell(1,4);
        for mon=1:tmpNMonitors
            plotHandles(mon) = plot(tmpFrascatiTimeAxes(mon,tmpStartSamples(mon):tmpEndSamples(mon)),squeeze(tmpStdPhases(mon,tmpStartSamples(mon):tmpEndSamples(mon))),monColours{mon});
            hold all;
            legStrs{mon} = sprintf('Mon%d (mean: %.2f^o)',mon,mean(tmpStdPhases(mon,tmpStartSamples(mon):tmpEndSamples(mon))));
        end
        plotHandles(4) = plot(tmpPetsTimeAxis(tmpPetsStartSample:tmpPetsEndSample),tmpPets_stdPhases(tmpPetsStartSample:tmpPetsEndSample),monColours{4});
        legStrs{4} = sprintf('PETS (mean: %.2f^o)',mean(tmpPets_stdPhases(tmpPetsStartSample:tmpPetsEndSample)));
        legend(plotHandles,legStrs);
        xlabel('Time [ns]');
        ylabel('Std Phase [degrees]');
        title('Std Phase vs. Sample No.');
        print([saveDir 'StdPhaseVsSampleNo.png'],'-dpng');

        %figure;
        plot(tmpResolution12,'b');
        hold all;
        plot(tmpResolution13,'r');
        plot(tmpResolution23,'g');
        title('Resolution vs. Sample No.');
        leg1Str = sprintf('Mon1-Mon2 (%.2f^o)',tmpMeanResolution12);
        leg2Str = sprintf('Mon1-Mon3 (%.2f^o)',tmpMeanResolution13);
        leg3Str = sprintf('Mon2-Mon3 (%.2f^o)',tmpMeanResolution23);
        legend(leg1Str,leg2Str,leg3Str);
        xlabel('Sample No.');
        ylabel('Resolution [12GHz Degrees]');
        print([saveDir 'Resolution.png'],'-dpng');
        hold off;
        
        %figure;
        plot(tmpMeanPulsePhases(1,:),tmpMeanPulsePhases(2,:),'o');
        xlabel('Mon1 Phase [12GHz Degrees]');
        ylabel('Mon2 Phase [12GHz Degrees]');
        title(sprintf('Mon1 - Mon2: Corr: %.2f, Grad: %.2f',tmpCorr12_meanPulsePhase,tmpFit12_meanPulsePhase));
        print([saveDir 'CorrelationMeanMon1Mon2.png'],'-dpng');
        hold off;
        
        %figure;
        plot(tmpMeanPulsePhases(1,:),tmpMeanPulsePhases(3,:),'o');
        hold all
        xlabel('Mon1 Phase [12GHz Degrees]');
        ylabel('Mon3 Phase [12GHz Degrees]');
        title(sprintf('Mon1 - Mon3: Corr: %.2f, Grad: %.2f',tmpCorr13_meanPulsePhase,tmpFit13_meanPulsePhase));
        print([saveDir 'CorrelationMeanMon1Mon3.png'],'-dpng');
        hold off;
        
        %figure;
        plot(tmpMeanPulsePhases(2,:),tmpMeanPulsePhases(3,:),'o');
        hold all
        xlabel('Mon2 Phase [12GHz Degrees]');
        ylabel('Mon3 Phase [12GHz Degrees]');
        title(sprintf('Mon2 - Mon3: Corr: %.2f, Grad: %.2f',tmpCorr23_meanPulsePhase,tmpFit23_meanPulsePhase));
        print([saveDir 'CorrelationMeanMon2Mon3.png'],'-dpng');
        hold off;
        
        %figure;
        plot(tmpMeanPulsePhases(1,:),tmpPets_meanPulsePhases,'o');
        xlabel('Mon3 Phase [12GHz Degrees]');
        ylabel('PETS Phase [12GHz Degrees]');
        title(sprintf('Mon1 - PETS: Corr: %.2f, Grad: %.2f',tmpCorr_PetsFrascati1_meanPulsePhase,tmpFit_PetsFrascati1_meanPulsePhase));
        print([saveDir 'CorrelationMeanMon1PETS.png'],'-dpng');
        hold off;
        
        %figure;
        plot(tmpMeanPulsePhases(2,:),tmpPets_meanPulsePhases,'o');
        xlabel('Mon3 Phase [12GHz Degrees]');
        ylabel('PETS Phase [12GHz Degrees]');
        title(sprintf('Mon2 - PETS: Corr: %.2f, Grad: %.2f',tmpCorr_PetsFrascati2_meanPulsePhase,tmpFit_PetsFrascati2_meanPulsePhase));
        print([saveDir 'CorrelationMeanMon2PETS.png'],'-dpng');
        hold off;
        
        %figure;
        plot(tmpMeanPulsePhases(3,:),tmpPets_meanPulsePhases,'o');
        xlabel('Mon3 Phase [12GHz Degrees]');
        ylabel('PETS Phase [12GHz Degrees]');
        title(sprintf('Mon3 - PETS: Corr: %.2f, Grad: %.2f',tmpCorr_PetsFrascati3_meanPulsePhase,tmpFit_PetsFrascati3_meanPulsePhase));
        print([saveDir 'CorrelationMeanMon3PETS.png'],'-dpng');
        hold off;
        
        %figure;
        plot(tmpCorr12_samples,'b');
        hold all;
        plot(tmpCorr13_samples,'r');
        plot(tmpCorr23_samples,'g');
        legend('Mon1-Mon2','Mon1-Mon3','Mon2-Mon3','Location','best');
        xlabel('Sample No.');
        ylabel('Corr Coeff');
        title('Correlation vs. Sample No.');
        print([saveDir 'CorrelationVsSampleNo.png'],'-dpng');
        hold off;
        
        %figure;
        plot(tmpCorr12_diffToMean,'b');
        hold all;
        plot(tmpCorr13_diffToMean,'r');
        plot(tmpCorr23_diffToMean,'g');
        legend('Mon1-Mon2','Mon1-Mon3','Mon2-Mon3','Location','best');
        xlabel('Sample No.');
        ylabel('Corr Coeff');
        title('Correlation in Difference to Mean Phase vs. Sample No.');
        print([saveDir 'CorrelationDiffToMeanVsSampleNo.png'],'-dpng');
        hold off;
        
        %figure;
        plot(tmpMeanPulsePhases(1,:),'b');
        hold all;
        plot(tmpMeanPulsePhases(2,:),'r')
        plot(tmpMeanPulsePhases(3,:),'g')
        plot(tmpPets_meanPulsePhases,'m');
        title('Mean Phase vs. Time');
        leg1 = sprintf('Mon1 (std = %.2f^o)',tmpStdPulsePhases(1)); 
        leg2 = sprintf('Mon2 (std = %.2f^o)',tmpStdPulsePhases(2)); 
        leg3 = sprintf('Mon3 (std = %.2f^o)',tmpStdPulsePhases(3)); 
        leg4 = sprintf('PETS (std = %.2f^o)',tmpPets_stdPulsePhase); 
        legend(leg1,leg2,leg3,leg4);
        xlabel('Time [Pulse No.]');
        ylabel('Phase [12GHz Degrees]');
        print([saveDir 'MeanPhaseVsTime.png'],'-dpng');
        hold off;
        
        %figure;
        plot(tmpMeanPulsePhases(2,:)-tmpMeanPulsePhases(1,:),'b');
        title(sprintf('Mon2-Mon1 Mean Phase vs. Time (std: %.2f^o)',std(tmpMeanPulsePhases(2,:)-tmpMeanPulsePhases(1,:))));
        xlabel('Time [Pulse No.]');
        ylabel('Phase [degrees]');
        ylim([-5 5]);
        print([saveDir 'DiffVsTime_Mon1Mon2.png'],'-dpng');
        hold off;
        
        %figure;
        plot(tmpMeanPulsePhases(3,:)-tmpMeanPulsePhases(1,:),'b');
        title(sprintf('Mon3-Mon1 Mean Phase vs. Time (std: %.2f^o)',std(tmpMeanPulsePhases(3,:)-tmpMeanPulsePhases(1,:))));
        xlabel('Time [Pulse No.]');
        ylabel('Phase [degrees]');
        ylim([-5 5]);
        print([saveDir 'DiffVsTime_Mon1Mon3.png'],'-dpng');
        hold off;
        
        %figure;
        plot(tmpMeanPulsePhases(3,:)-tmpMeanPulsePhases(2,:),'b');
        title(sprintf('Mon3-Mon2 Mean Phase vs. Time (std: %.2f^o)',std(tmpMeanPulsePhases(3,:)-tmpMeanPulsePhases(2,:))));
        xlabel('Time [Pulse No.]');
        ylabel('Phase [degrees]');
        ylim([-5 5]);
        print([saveDir 'DiffVsTime_Mon2Mon3.png'],'-dpng');
        hold off;
        
        %figure;
        plot(tmpPets_meanPulsePhases'-tmpMeanPulsePhases(3,:),'b');
        title(sprintf('PETS-Mon3 Mean Phase vs. Time (std: %.2f^o)',std(tmpPets_meanPulsePhases'-tmpMeanPulsePhases(3,:))));
        xlabel('Time [Pulse No.]');
        ylabel('Phase [degrees]');
        ylim([-5 5]);
        print([saveDir 'DiffVsTime_PETSMon3.png'],'-dpng');
        hold off;
        
        %figure;
        plot(tmpPets_meanPulsePhases'-tmpMeanPulsePhases(1,:),'b');
        title(sprintf('PETS-Mon1 Mean Phase vs. Time (std: %.2f^o)',std(tmpPets_meanPulsePhases'-tmpMeanPulsePhases(1,:))));
        xlabel('Time [Pulse No.]');
        ylabel('Phase [degrees]');
        ylim([-5 5]);
        print([saveDir 'DiffVsTime_PETSMon1.png'],'-dpng');
        hold off;
        
        %figure;
        plot(tmpPets_meanPulsePhases'-tmpMeanPulsePhases(2,:),'b');
        title(sprintf('PETS-Mon2 Mean Phase vs. Time (std: %.2f^o)',std(tmpPets_meanPulsePhases'-tmpMeanPulsePhases(2,:))));
        xlabel('Time [Pulse No.]');
        ylabel('Phase [degrees]');
        ylim([-5 5]);
        print([saveDir 'DiffVsTime_PETSMon2.png'],'-dpng');
        hold off;
    end

end % end of dataset loop

%% Data set comparison plots

% Mean Phase along the pulse
% subplot(3,1,1:2);
figure;
toPlot = cell(1,nDataSets);
legstr = {};
for i=1:nDataSets
    if (~skipPhase{i})
        toPlot{i} = squeeze(meanPhases{i}(3,startSamples{i}(3):endSamples{i}(3)));
        %toPlot{i} = toPlot{i} - mean(toPlot{i}([1:20 70:90]));
        plot(frascatiTimeAxes{i}(3,startSamples{i}(3):endSamples{i}(3)), toPlot{i},'Color',dataSetColours{i},'LineWidth',2,'LineStyle',dataSetStyles{i})
        hold all
        legstr{length(legstr)+1} = dataSetLabels{i};
    end
end
legend(legstr);
ylabel('Phase [degrees]');
xlabel('Time [ns]');
grid on;
% set(gca,'Ytick',-6:2:10);
% ylim([-6 10])
title('Mon3 Pulse Phase (mean)');

figure;
toPlot = cell(1,nDataSets);
legstr = {};
for i=1:nDataSets
    if (~skipPhase{i})
        toPlot{i} = squeeze(meanPhases{i}(1,startSamples{i}(1):endSamples{i}(1)));
        %toPlot{i} = toPlot{i} - mean(toPlot{i}([1:20 70:90]));
        plot(frascatiTimeAxes{i}(1,startSamples{i}(1):endSamples{i}(1)), toPlot{i},'Color',dataSetColours{i},'LineWidth',2,'LineStyle',dataSetStyles{i})
        hold all
        legstr{length(legstr)+1} = dataSetLabels{i};
    end
end
legend(legstr);
ylabel('Phase [degrees]');
xlabel('Time [ns]');
grid on;
% set(gca,'Ytick',-6:2:10);
% ylim([-6 10])
title('Mon1 Pulse Phase (mean)');


% subplot(3,1,3);
% plot(frascatiTimeAxes{2}(3,startSamples{i}(3):endSamples{i}(3)), toPlot{2}-toPlot{1},'Color',dataSetColours{2},'LineWidth',2);
% hold all;
% plot(frascatiTimeAxes{3}(3,startSamples{i}(3):endSamples{i}(3)), toPlot{3}-toPlot{1},dataSetColours{3},'LineWidth',3);
% xlabel('Time [ns]')
% ylabel('Phase [degrees]');
% legend(dataSetLabels{2:3});
% grid on;
% title('Difference to Nominal Phase');


figure;
for i=1:nDataSets
    plot(meanBPMCC930_H{i},'Color',dataSetColours{i},'LineWidth',2,'LineStyle',dataSetStyles{i});
    hold all;
end
legend(dataSetLabels);
xlabel('Sample No.');
ylabel('Position [mm]');
title('Feedforward Applied to One Kicker Only')

%%
% bpmSampleRange = 420:500;
% bpmSampleRate = 7.7;
% phaseSampleRange = 84:237;%106:210;
% bpmSubtractTime = 3638;
% phaseSubtractTime = 730;
% matchSign = [1,1,-1,1,-1,1];
% 
% phaseTrace = squeeze(meanPhases{1}(1,:));
% phaseTrace = phaseTrace./max(phaseTrace);
% phaseTime = squeeze(frascatiTimeAxes{1}(1,:))-phaseSubtractTime;
% figure;
% plotColours = {[],'b',[0,0.5,0],'r',[0.8,0.2,0.8],[]};
% for i=2:nDataSets
%     bpmTrace = meanBPMCC930_H{i}-meanBPMCC930_H{1};
%     bpmTrace = matchSign(i)*bpmTrace;
%     bpmTrace = bpmTrace./max(bpmTrace);
%     bpmTime = bpmSampleRate*(1:length(bpmTrace))-bpmSubtractTime;
%     plot(bpmTime(bpmSampleRange),bpmTrace(bpmSampleRange),'-','Color',plotColours{i},'LineWidth',2)
%     hold all;
% end
% plot(phaseTime(phaseSampleRange),phaseTrace(phaseSampleRange),'Color','k','LineWidth',2,'LineStyle','-');
% legend(dataSetLabels{2:nDataSets})
% xlabel('Time [ns]');
% ylabel('Output [a.u.]');
% ylim([-1.8 1.2])
% grid on;

%%
figure;
for i=1:nDataSets
    plot(squeeze(stdPhases{i}(3,:)),'Color',dataSetColours{i},'LineWidth',2,'LineStyle',dataSetStyles{i});
    hold all;
end
legend(dataSetLabels)
title('Mon3');

figure;
for i=1:nDataSets
    plot(pets_stdPhases{i},'Color',dataSetColours{i},'LineWidth',2,'LineStyle',dataSetStyles{i});
    hold all;
end
legend(dataSetLabels)
title('PETS')

figure;
for i=1:nDataSets
    plot(corr23_samples{i},'Color',dataSetColours{i},'LineWidth',2,'LineStyle',dataSetStyles{i});
    hold all;
end
legend(dataSetLabels)
title('Corr Mon2-Mon3');

% %%
% figure;
% toPlot = cell(1,nDataSets);
% legstr = {};
% for i=1:nDataSets
%     if (~skipPhase{i})
%         toPlot{i} = squeeze(meanPhases{i}(3,startSamples{i}(3):endSamples{i}(3)));
%         %toPlot{i} = toPlot{i} - mean(toPlot{i}([1:20 70:90]));
%         plot(toPlot{i},'Color',dataSetColours{i},'LineWidth',2,'LineStyle',dataSetStyles{i})
%         hold all
%         legstr{length(legstr)+1} = dataSetLabels{i};
%     end
% end
% legend(legstr);
% ylabel('Phase [degrees]');
% xlabel('Sample No.');
% grid on;
% % set(gca,'Ytick',-6:2:10);
% % ylim([-6 10])
% title('Mon3 Pulse Phase (mean)');