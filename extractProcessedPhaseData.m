function [processedPhaseData] = extractProcessedPhaseData( CTFData, calibrationConstants )

    monitorNames = {'Mon1','Mon2','Mon3'};

    calibrationAmplitudes = calibrationConstants(:,1);
    calibrationOffsets = calibrationConstants(:,4);

    phaseSignalSkipPercent = 10; % skip this percentage of the found start and end of the pulse.

    nMonitors = 3;

    petsStartSample = 550;
    petsEndSample = 800;

    removeBadPulses = true; % if true remove bad pulses
    removeThreshold = 0.01; % remove pulses that are removeThreshold % away from max on diode
    removeMonIndices = [1, 2]; % look at diode for these monitors when deciding which pulses to remove (require good transmission for all monitors with give indices)

    %% Load data, remove bad pulses
    [mixers,diodes] = extractMixerDiode(CTFData);
    [~,nPulses,nSamples] = size(mixers);

    meanDiodes = squeeze(nanmean(diodes,2));
    meanMixers = squeeze(nanmean(mixers,2));

    % find pulse range using diode
    startSamples = NaN*ones(1,nMonitors);
    endSamples = NaN*ones(1,nMonitors);
    for mon=1:nMonitors
        [startSamples(mon),endSamples(mon)] = getSignalRange(meanDiodes(mon,:),phaseSignalSkipPercent);
    end

    % to make a quick phase difference etc. calculation, need the number of
    % samples in each range to be the same. Fix it here.
    sampleRanges = (endSamples-startSamples)+1;
    minSampleRange = min(sampleRanges);
    diffSampleRange = sampleRanges - minSampleRange;
    for mon=1:nMonitors
        if (diffSampleRange(mon) > 0)
            if (mod(diffSampleRange(mon),2) == 0)
                skipSamples = diffSampleRange(mon)./2;
                startSamples(mon) = startSamples(mon) + skipSamples;
                endSamples(mon) = endSamples(mon) - skipSamples;
            else
                skipSamples = floor(diffSampleRange(mon)./2);
                startSamples(mon) = startSamples(mon) + skipSamples + 1;
                endSamples(mon) = endSamples(mon) - skipSamples;
            end
        end
    end

    if (removeBadPulses==true && nPulses>=20) % require at least 20 pulses just so we actually have something to take a decent mean from
        transmission = cell(1,nMonitors);
        isGoodPulseMon = cell(1,nMonitors);
        for mon = 1:nMonitors
            transmission{mon} = squeeze(nanmean(diodes(mon,:,startSamples(mon):endSamples(mon)),3));
            isGoodPulseMon{mon} = (transmission{mon}-min(transmission{mon})) < removeThreshold;
        end

        isGoodPulseAll = ones(1,nPulses);
        for mon = removeMonIndices
            isGoodPulseAll = isGoodPulseAll & isGoodPulseMon{mon};
        end

        CTFData = CTFData(isGoodPulseAll);
        mixers = mixers(:,isGoodPulseAll,:);
        diodes = diodes(:,isGoodPulseAll,:);
        nPulses = length(CTFData);

    end

    %% Process Frascati data

    phases = NaN*ones(nMonitors,nPulses,nSamples);
    for mon = 1:nMonitors
        phases(mon,:,:) = getPhaseMixerDiode(mixers(mon,:,:), diodes(mon,:,:), calibrationAmplitudes(mon), calibrationOffsets(mon)); 
    end


    meanPhases = squeeze(nanmean(phases,2)); % mean per sample
    stdPhases = squeeze(nanstd(phases,0,2)); % std per sample
    meanStdPhase = NaN*ones(1,nMonitors);
    for mon=1:nMonitors
        meanStdPhase(mon) = squeeze(nanmean(stdPhases(mon,startSamples(mon):endSamples(mon)),2)); % mean sample std
    end

    nAvg = 10;
    meanStdPhase_nAvgPulses = NaN*ones(nMonitors,ceil(nPulses/nAvg));
    for p=1:10:nPulses
        if p+nAvg-1 > nPulses
            pulseRange = p:nPulses;
        else
            pulseRange = p:p+nAvg-1;
        end

        for mon=1:nMonitors
            tmpStdPhase = squeeze(nanstd(phases(mon,pulseRange,:),0,2));
            meanStdPhase_nAvgPulses(mon,ceil(p/nAvg)) = nanmean(tmpStdPhase(startSamples(mon):endSamples(mon)));
        end
    end

    diffPhases12 = phases(1,:,startSamples(1):endSamples(1)) - phases(2,:,startSamples(2):endSamples(2));
    diffPhases13 = phases(1,:,startSamples(1):endSamples(1)) - phases(3,:,startSamples(3):endSamples(3));
    diffPhases23 = phases(2,:,startSamples(2):endSamples(2)) - phases(3,:,startSamples(3):endSamples(3));

    resolution12 = squeeze(nanstd(diffPhases12,0,3))./sqrt(2);
    resolution13 = squeeze(nanstd(diffPhases13,0,3))./sqrt(2);
    resolution23 = squeeze(nanstd(diffPhases23,0,3))./sqrt(2);

    meanPulseDiodes = NaN*ones(nMonitors,nPulses);
    meanPulseMixers = NaN*ones(nMonitors,nPulses);
    meanPulsePhases = NaN*ones(nMonitors,nPulses); % mean per pulse
    for mon = 1:nMonitors
        meanPulseDiodes(mon,:) = squeeze(nanmean(diodes(mon,:,startSamples(mon):endSamples(mon)),3));
        meanPulseMixers(mon,:) = squeeze(nanmean(mixers(mon,:,startSamples(mon):endSamples(mon)),3));
        meanPulsePhases(mon,:) = squeeze(nanmean(phases(mon,:,startSamples(mon):endSamples(mon)),3));
    end
    stdPulseDiodes =std(meanPulseDiodes,0,2); 
    stdPulseMixers =std(meanPulseMixers,0,2);
    stdPulsePhases =std(meanPulsePhases,0,2); % error bars on mean pulse phase

    corr12_meanPulsePhase = corrcoef(meanPulsePhases(1,:),meanPulsePhases(2,:));
    corr13_meanPulsePhase = corrcoef(meanPulsePhases(1,:),meanPulsePhases(3,:));
    corr23_meanPulsePhase = corrcoef(meanPulsePhases(2,:),meanPulsePhases(3,:));
    corr12_meanPulsePhase = corr12_meanPulsePhase(1,2);
    corr13_meanPulsePhase = corr13_meanPulsePhase(1,2);
    corr23_meanPulsePhase = corr23_meanPulsePhase(1,2);

    meanPhases_offsetSubtract = NaN*ones(nMonitors,nSamples);
    for mon=1:nMonitors
        meanPhases_offsetSubtract(mon,:) = meanPhases(mon,:)-nanmean(meanPhases(mon,startSamples(mon):endSamples(mon)));
    end
    meanDiffPhase12 = meanPhases_offsetSubtract(1,startSamples(1):endSamples(1)) - meanPhases_offsetSubtract(2,startSamples(2):endSamples(2));
    meanDiffPhase13 = meanPhases_offsetSubtract(1,startSamples(1):endSamples(1)) - meanPhases_offsetSubtract(3,startSamples(3):endSamples(3));
    meanDiffPhase23 = meanPhases_offsetSubtract(2,startSamples(2):endSamples(2)) - meanPhases_offsetSubtract(3,startSamples(3):endSamples(3));

    %% Process PETS data
    
    pets_I = extractCTFSignalFromMergedData('CE.SCOPE03.CH01.Acquisition.value.value',CTFData);
    pets_Q = extractCTFSignalFromMergedData('CE.SCOPE03.CH02.Acquisition.value.value',CTFData);
    
    pets_I_Sensitivity = extractCTFSignalFromMergedData('CE.SCOPE03.CH01.Acquisition.sensitivity.value',CTFData);
    pets_Q_Sensitivity = extractCTFSignalFromMergedData('CE.SCOPE03.CH02.Acquisition.sensitivity.value',CTFData);
    
    pets_I = pets_I.*pets_I_Sensitivity(1);
    pets_Q = pets_Q.*pets_Q_Sensitivity(1);
    
    [petsStartSample_I,petsEndSample_I] = getSignalRange(pets_I,phaseSignalSkipPercent);
    [petsStartSample_Q,petsEndSample_Q] = getSignalRange(pets_Q,phaseSignalSkipPercent);
    petsStartSample = max([petsStartSample_I petsStartSample_Q]);
    petsEndSample = min([petsEndSample_I petsEndSample_Q]);
    
    [pets_Phase,pets_Power] = getPhaseIQ(pets_I, pets_Q);

    pets_meanPulsePhases = nanmean(pets_Phase(:,petsStartSample:petsEndSample),2); % mean per pulse
    pets_meanPhases = nanmean(pets_Phase); % mean per sample
    pets_stdPhases = nanstd(pets_Phase); % std per sample
    pets_meanStdPhase = nanmean(pets_stdPhases(petsStartSample:petsEndSample)); % mean sample std

    corr_PetsFrascati3_meanPulsePhase = corrcoef(meanPulsePhases(3,:),pets_meanPulsePhases);
    corr_PetsFrascati3_meanPulsePhase = corr_PetsFrascati3_meanPulsePhase(1,2);
    
    %% Process BPR data
    
    %% Save data to struct.
    % At the moment not in a clever way, if anything is added above that
    % should be output it needs to be added below.
    processedPhaseData = struct();
    processedPhaseData.calibrationConstants = calibrationConstants;
    processedPhaseData.corr12_meanPulsePhase = corr12_meanPulsePhase;
    processedPhaseData.corr13_meanPulsePhase = corr13_meanPulsePhase;
    processedPhaseData.corr23_meanPulsePhase = corr23_meanPulsePhase;
    processedPhaseData.corr_PetsFrascati3_meanPulsePhase = corr_PetsFrascati3_meanPulsePhase;
    processedPhaseData.diffPhases12 = diffPhases12;
    processedPhaseData.diffPhases13 = diffPhases13;
    processedPhaseData.diffPhases23 = diffPhases23;
    processedPhaseData.diffSampleRange = diffSampleRange;
    processedPhaseData.diodes = diodes;
    processedPhaseData.endSamples = endSamples;
    processedPhaseData.isGoodPulseAll = isGoodPulseAll;
    processedPhaseData.isGoodPulseMon = isGoodPulseMon;
    processedPhaseData.meanDiffPhase12 = meanDiffPhase12;
    processedPhaseData.meanDiffPhase13 = meanDiffPhase13;
    processedPhaseData.meanDiffPhase23 = meanDiffPhase23;
    processedPhaseData.meanDiodes = meanDiodes;
    processedPhaseData.meanMixers = meanMixers;
    processedPhaseData.meanPhases = meanPhases;
    processedPhaseData.meanPhases_offsetSubtract = meanPhases_offsetSubtract;
    processedPhaseData.meanPulseDiodes = meanPulseDiodes;
    processedPhaseData.meanPulseMixers = meanPulseMixers;
    processedPhaseData.meanPulsePhases = meanPulsePhases;
    processedPhaseData.meanStdPhase = meanStdPhase;
    processedPhaseData.meanStdPhase_nAvgPulses = meanStdPhase_nAvgPulses;
    processedPhaseData.minSampleRange = minSampleRange;
    processedPhaseData.mixers = mixers;
    processedPhaseData.monitorNames = monitorNames;
    processedPhaseData.nAvg = nAvg;
    processedPhaseData.nMonitors = nMonitors;
    processedPhaseData.nPulses = nPulses;
    processedPhaseData.nSamples = nSamples;
    processedPhaseData.petsEndSample = petsEndSample;
    processedPhaseData.petsEndSample_I = petsEndSample_I;
    processedPhaseData.petsEndSample_Q = petsEndSample_Q;
    processedPhaseData.petsStartSample = petsStartSample;
    processedPhaseData.petsStartSample_I = petsStartSample_I;
    processedPhaseData.petsStartSample_Q = petsStartSample_Q;
    processedPhaseData.pets_I = pets_I;
    processedPhaseData.pets_Phase = pets_Phase;
    processedPhaseData.pets_Power = pets_Power;
    processedPhaseData.pets_Q = pets_Q;
    processedPhaseData.pets_meanPhases = pets_meanPhases;
    processedPhaseData.pets_meanPulsePhases = pets_meanPulsePhases;
    processedPhaseData.pets_meanStdPhase = pets_meanStdPhase;
    processedPhaseData.pets_stdPhases = pets_stdPhases;
    processedPhaseData.phaseSignalSkipPercent = phaseSignalSkipPercent;
    processedPhaseData.phases = phases;
    processedPhaseData.pulseRange = pulseRange;
    processedPhaseData.removeBadPulses = removeBadPulses;
    processedPhaseData.removeMonIndices = removeMonIndices;
    processedPhaseData.removeThreshold = removeThreshold;
    processedPhaseData.resolution12 = resolution12;
    processedPhaseData.resolution13 = resolution13;
    processedPhaseData.resolution23 = resolution23;
    processedPhaseData.startSamples = startSamples;
    processedPhaseData.stdPhases = stdPhases;
    processedPhaseData.stdPulseDiodes = stdPulseDiodes;
    processedPhaseData.stdPulseMixers = stdPulseMixers;
    processedPhaseData.stdPulsePhases = stdPulsePhases;
    processedPhaseData.tmpStdPhase = tmpStdPhase;
    processedPhaseData.transmission = transmission;
     
    try 
        processedPhaseData.dataDir = dataDir;
        processedPhaseData.dataSetComment = CTFData(1).comment;
    catch 
    end

end