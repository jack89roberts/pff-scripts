%% processData_Frascati.m
% All signal process, calculations related to the Frascati monitors (only).
% This script should only be called from processData.m, not used by itself.
% June 2015, Jack Roberts

% Required inputs (all from processData.m):
    % sampleRange 
        % Used to calculate means etc. If this is left empty, will prompt to ask for sample range
    % CTFData or FONTData
        % Loaded merged data.
    % calibrationConstants
    % useMixerOverSqrtDiode
        % Calibration information.
        
%% signal processing: frascati monitors

if (isFONTData)
    [ mixers, diodes ] = extractFONTMixerDiode(FONTData.ADCs, 1, mixerADCs, diodeIsPresent, diodeADCs);
else
    [ mixers, diodes ] = extractMixerDiode(CTFData);
end

[nMons,nPulses,nSamples] = size(mixers);
baselineStartSamples = 1:round(nSamples/10);
baselineEndSamples = (nSamples-round(nSamples/10)):nSamples;

% calculate phases
phases = NaN(nMons,nPulses,nSamples);
for mon=1:nMons
    CalibrationConstant = calibrationConstants(mon,1);
    CalibrationOffset = calibrationConstants(mon,4);
    Mixer =  squeeze(mixers(mon,:,:));
    if (diodeIsPresent)
        Diode = squeeze(diodes(mon,:,:));
    else
        Diode = [];
    end
    phases(mon,:,:) = getPhaseMixerDiode( Mixer, Diode, CalibrationConstant, CalibrationOffset, useMixerOverSqrtDiode );

end

% remove duplicate pulses (helps to eradicate any acquisition issues - e.g.
% too slow acquisition, acquisition all zeroes etc.)
if (stripDuplicates)
    for mon=1:nMons
        if (diodeIsPresent)
            Diode = squeeze(diodes(mon,:,:));
            diodes(mon,:,:) = removeDuplicatePulses(Diode);
        end
        
        Mixer = squeeze(mixers(mon,:,:));
        Phase = squeeze(phases(mon,:,:));
        
        mixers(mon,:,:) = removeDuplicatePulses(Mixer);
        phases(mon,:,:) = removeDuplicatePulses(Phase);
    end
end

% extract baseline for noise check before aligning
% signals which strips it
if (diodeIsPresent)
    diodeBaselineStart = diodes(:,:,baselineStartSamples);
    diodeBaselineEnd = diodes(:,:,baselineEndSamples);
end
mixerBaselineStart = mixers(:,:,baselineStartSamples);
mixerBaselineEnd = mixers(:,:,baselineEndSamples);
phaseBaselineStart = phases(:,:,baselineStartSamples);
phaseBaselineEnd = phases(:,:,baselineEndSamples);

% align mixers, diodes, phases
pulseSampleRange = cell(1,nMons);

for mon=1:nMons
    if (diodeIsPresent)
        Diode = squeeze(diodes(mon,:,:));
    end
    Mixer = squeeze(mixers(mon,:,:));
    Phase = squeeze(phases(mon,:,:));
    
    if (~useStaticFrascatiDelay)
        %[Diode,others] = getAligned2(Diode,20,nSamples,1,{Mixer,Phase});
        [Diode,others,pulseSampleRange{mon}] = getAlignedXCorr(Diode,alignTo,{Mixer,Phase}); % align start/end of pulses for each monitor to be at the sample place
        diodes(mon,:,:) = Diode;
        mixers(mon,:,:) = others{1};
        phases(mon,:,:) = others{2};
    elseif (diodeIsPresent)
        [Diode,others,pulseSampleRange{mon}] = getAlignedXCorr(Diode,'first',{Mixer,Phase}); % align all pulses to the first pulse for each monitor (but don't move pulse to start at sample x)
        diodes(mon,:,:) = Diode;
        mixers(mon,:,:) = others{1};
        phases(mon,:,:) = others{2};      
    else
        pulseSampleRange = [];
    end
    
end

if (useStaticFrascatiDelay) % use static delays to align signals now
    if (strcmp(staticFrascatiDelayType,'calculate') && diodeIsPresent)
        staticFrascatiDelay = calculateDiodeDelay(diodes);
    elseif (strcmp(staticFrascatiDelayType,'plot'))
        % PLOT AND PROMPT TO CHOOSE DELAY (AND PICK SAMPLE RANGE)
        [ pulseSampleRange, staticFrascatiDelay ] = pickPulseRange( mixers,diodes,phases );        
    end
    
    for mon=1:nMons
        if (diodeIsPresent)
            diodes(mon,:,:) = delaySignal(squeeze(diodes(mon,:,:)),staticFrascatiDelay(mon));
        end
        mixers(mon,:,:) = delaySignal(squeeze(mixers(mon,:,:)),staticFrascatiDelay(mon));
        phases(mon,:,:) = delaySignal(squeeze(phases(mon,:,:)),staticFrascatiDelay(mon));
        pulseSampleRange{mon} = pulseSampleRange{mon} + staticFrascatiDelay(mon);
    end
end

% For PFF On FONTData, remove downstream drifts caused by changes made to
% channel offset.
if (isFONTData && stripChanOffset)
    FFStatus = FONTData.FFStatus;
    ChanOffset = FONTData.ChanOffset;
    isFFOn = FFStatus==1;
    if (sum(isFFOn)>0)
        tmpChanOffset = ChanOffset(isFFOn);
        tmpMon3Phase = squeeze(phases(3,isFFOn,:));
        tmpMon3Phase = subtractChanOffset(tmpMon3Phase,...
                                          tmpChanOffset,...
                                          calibrationConstants(1,1),...
                                          calibrationConstants(1,4));

        phases(3,isFFOn,:) = tmpMon3Phase;

    else
        fprintf('WARNING: Requested to strip channel offset, but found no FF on data so did nothing.\n');
    end
    
end

% Make a time axis based on sample frequency, with zero at found pulse
% start for Mon3 (if present) or Mon1/2 (if not).
if (isFONTData)
    if (isempty(sampInterval))
        sampInterval = 1000/357; % FONT clock frequency
    end
else
    try
        sampInterval = extractCTFSignalFromMergedData('CT_SCOPE01_CH01.Acquisition.sampleInterval.value',CTFData);
    catch
        sampInterval = extractCTFSignalFromMergedData('CT.SCOPE01.CH01.Acquisition.sampleInterval.value',CTFData);        
    end
    sampInterval = sampInterval(~isnan(sampInterval));
    sampInterval = sampInterval(1);
end
timeAxis = sampInterval.*(0:(nSamples-1));    
try
    timeRefSample = pulseSampleRange{3}(1);
catch
    timeRefSample = NaN;
end
if (isnan(timeRefSample))
    timeRefSample = max([pulseSampleRange{1}(1) pulseSampleRange{2}(1)]);
end
timeAxis = timeAxis-timeAxis(timeRefSample);

% length of pulse in ns
pulseLength = NaN(1,nMons);
for mon=1:nMons
    pulseLength(mon) = sampInterval.*(length(pulseSampleRange{mon})-1);
end

% prompt for sample range to use if not given
if (isempty(sampleRange))
    a = figure();
    
    for mon=1:nMons
        if (diodeIsPresent)
            subplot(1,2,1)
        end
        plot(squeeze(mixers(mon,1,:)));
        hold all;
        if (diodeIsPresent)
            subplot(1,2,2)
            plot(squeeze(diodes(mon,1,:)));
            hold all;
        end
    end
    
    if (diodeIsPresent)
        subplot(1,2,1);
    end
    title('MIXERS');
    if (diodeIsPresent)
        subplot(1,2,2);
        title('DIODES')
    end
    
    startSample = input('Start sample: ');
    endSample = input('End sample: ');
    sampleRange = startSample:endSample;
    
    try
        close(a);
    catch
    end
end

% delay between start/end pulse and sample range
startDelay = NaN(1,nMons);
endDelay = NaN(1,nMons);
for mon=1:nMons
    try
        startDelay(mon) = sampInterval.*(sampleRange(1)-pulseSampleRange{mon}(1));
        endDelay(mon) = sampInterval.*(pulseSampleRange{mon}(end)-sampleRange(end));
    catch   
    end
end

startDelayCT = nanmean(startDelay(1:2));
startDelayCB = startDelay(3);
endDelayCT = nanmean(endDelay(1:2));
endDelayCB = endDelay(3);
if (isnan(startDelayCT)&&isnan(startDelayCB))
    error('Problem defining pulse length');
elseif (isnan(startDelayCT))
    startDelayCT = startDelayCB;
    endDelayCT = endDelayCB;
elseif(isnan(startDelayCB))
    startDelayCB = startDelayCT;
    endDelayCB = endDelayCT;
end

% remove bad pulses: outside 3 sigma on diode (will remove pulses with bad transmission)
if (diodeIsPresent)
    for mon=1:nMons
        Diode = squeeze(diodes(mon,:,:));
        Mixer = squeeze(mixers(mon,:,:));
        Phase = squeeze(phases(mon,:,:));
        dioBaseSt = squeeze(diodeBaselineStart(mon,:,:));
        dioBaseEn = squeeze(diodeBaselineEnd(mon,:,:));
        mixBaseSt = squeeze(mixerBaselineStart(mon,:,:));
        mixBaseEn = squeeze(mixerBaselineEnd(mon,:,:));
        phBaseSt= squeeze(phaseBaselineStart(mon,:,:));
        phBaseEn= squeeze(phaseBaselineEnd(mon,:,:));

        [Diode,others] = removeBadPulses(Diode,sampleRange,{Mixer,dioBaseSt,dioBaseEn,mixBaseSt,mixBaseEn,phBaseSt,phBaseEn,Phase});
        diodes(mon,:,:) = Diode;
        mixers(mon,:,:) = others{1};
        diodeBaselineStart(mon,:,:) = others{2};
        diodeBaselineEnd(mon,:,:) = others{3};
        mixerBaselineStart(mon,:,:) = others{4};
        mixerBaselineEnd(mon,:,:) = others{5};
        phaseBaselineStart(mon,:,:) = others{6};
        phaseBaselineEnd(mon,:,:) = others{7};
        phases(mon,:,:) = others{8};
    end
end

% remove drift (only phases)
if (stripDrift)
    for mon=1:nMons
        Phase = squeeze(phases(mon,:,:));
        phases(mon,:,:) = removeDrift(Phase,sampleRange,driftNAvg);
    end    
end

% remove bad pulses: outside 3 sigma on phase
if (stripOutliers)
    for mon=1:nMons
        phases(mon,:,:) = removeBadPulses(squeeze(phases(mon,:,:)),sampleRange);
    end
end
%% baseline noise: mean, std, resolution

% diode
if (diodeIsPresent)
    [meanSampleDioBaseSt,stdSampleDioBaseSt,meanSampleDioBaseSt_err,stdSampleDioBaseSt_err] = nanMeanStdErr(diodeBaselineStart,2);
    [meanSampleDioBaseEn,stdSampleDioBaseEn,meanSampleDioBaseEn_err,stdSampleDioBaseEn_err] = nanMeanStdErr(diodeBaselineEnd,2);

    meanBaselineDiodeStart = NaN(1,nMons);
    stdBaselineDiodeStart = NaN(1,nMons);
    meanBaselineDiodeEnd = NaN(1,nMons);
    stdBaselineDiodeEnd = NaN(1,nMons);
    meanBaselineDiodeStart_err = NaN(1,nMons);
    stdBaselineDiodeStart_err = NaN(1,nMons);
    meanBaselineDiodeEnd_err = NaN(1,nMons);
    stdBaselineDiodeEnd_err = NaN(1,nMons);
    for mon=1:nMons
        tmpSt = squeeze(diodeBaselineStart(mon,:,:));
        tmpEn = squeeze(diodeBaselineEnd(mon,:,:));

        [meanBaselineDiodeStart(mon),...
         stdBaselineDiodeStart(mon),...
         meanBaselineDiodeStart_err(mon),...
         stdBaselineDiodeStart_err(mon)] = nanMeanStdErr(tmpSt(:));

        [meanBaselineDiodeEnd(mon),...
         stdBaselineDiodeEnd(mon),...
         meanBaselineDiodeEnd_err(mon),...
         stdBaselineDiodeEnd_err(mon)] = nanMeanStdErr(tmpEn(:));
    end

    res12DioBaseSt = squeeze(diodeBaselineStart(1,:,:)-diodeBaselineStart(2,:,:));
    res13DioBaseSt = squeeze(diodeBaselineStart(1,:,:)-diodeBaselineStart(3,:,:));
    res23DioBaseSt = squeeze(diodeBaselineStart(2,:,:)-diodeBaselineStart(3,:,:));
    [~,res12DioBaseSt,~,res12DioBaseSt_err] = nanMeanStdErr(res12DioBaseSt./sqrt(2),1);
    [~,res13DioBaseSt,~,res13DioBaseSt_err] = nanMeanStdErr(res13DioBaseSt./sqrt(2),1);
    [~,res23DioBaseSt,~,res23DioBaseSt_err] = nanMeanStdErr(res23DioBaseSt./sqrt(2),1);

    res12DioBaseEn = squeeze(diodeBaselineEnd(1,:,:)-diodeBaselineEnd(2,:,:));
    res13DioBaseEn = squeeze(diodeBaselineEnd(1,:,:)-diodeBaselineEnd(3,:,:));
    res23DioBaseEn = squeeze(diodeBaselineEnd(2,:,:)-diodeBaselineEnd(3,:,:));
    [~,res12DioBaseEn,~,res12DioBaseEn_err] = nanMeanStdErr(res12DioBaseEn./sqrt(2),1);
    [~,res13DioBaseEn,~,res13DioBaseEn_err] = nanMeanStdErr(res13DioBaseEn./sqrt(2),1);
    [~,res23DioBaseEn,~,res23DioBaseEn_err] = nanMeanStdErr(res23DioBaseEn./sqrt(2),1);

    [meanRes12DioBaseSt,~,meanRes12DioBaseSt_err] = nanMeanStdErr(res12DioBaseSt);
    [meanRes13DioBaseSt,~,meanRes13DioBaseSt_err] = nanMeanStdErr(res13DioBaseSt);
    [meanRes23DioBaseSt,~,meanRes23DioBaseSt_err] = nanMeanStdErr(res23DioBaseSt);
    [meanRes12DioBaseEn,~,meanRes12DioBaseEn_err] = nanMeanStdErr(res12DioBaseEn);
    [meanRes13DioBaseEn,~,meanRes13DioBaseEn_err] = nanMeanStdErr(res13DioBaseEn);
    [meanRes23DioBaseEn,~,meanRes23DioBaseEn_err] = nanMeanStdErr(res23DioBaseEn);
end

% mixer
[meanSampleMixBaseSt,stdSampleMixBaseSt,meanSampleMixBaseSt_err,stdSampleMixBaseSt_err] = nanMeanStdErr(mixerBaselineStart,2);
[meanSampleMixBaseEn,stdSampleMixBaseEn,meanSampleMixBaseEn_err,stdSampleMixBaseEn_err] = nanMeanStdErr(mixerBaselineEnd,2);

meanBaselineMixerStart = NaN(1,nMons);
stdBaselineMixerStart = NaN(1,nMons);
meanBaselineMixerEnd = NaN(1,nMons);
stdBaselineMixerEnd = NaN(1,nMons);
meanBaselineMixerStart_err = NaN(1,nMons);
stdBaselineMixerStart_err = NaN(1,nMons);
meanBaselineMixerEnd_err = NaN(1,nMons);
stdBaselineMixerEnd_err = NaN(1,nMons);
for mon=1:nMons
    tmpSt = squeeze(mixerBaselineStart(mon,:,:));
    tmpEn = squeeze(mixerBaselineEnd(mon,:,:));
    
    [meanBaselineMixerStart(mon),...
     stdBaselineMixerStart(mon),...
     meanBaselineMixerStart_err(mon),...
     stdBaselineMixerStart_err(mon)] = nanMeanStdErr(tmpSt(:));
    
    [meanBaselineMixerEnd(mon),...
     stdBaselineMixerEnd(mon),...
     meanBaselineMixerEnd_err(mon),...
     stdBaselineMixerEnd_err(mon)] = nanMeanStdErr(tmpEn(:));
end

res12MixBaseSt = squeeze(mixerBaselineStart(1,:,:)-mixerBaselineStart(2,:,:));
res13MixBaseSt = squeeze(mixerBaselineStart(1,:,:)-mixerBaselineStart(3,:,:));
res23MixBaseSt = squeeze(mixerBaselineStart(2,:,:)-mixerBaselineStart(3,:,:));
[~,res12MixBaseSt,~,res12MixBaseSt_err] = nanMeanStdErr(res12MixBaseSt./sqrt(2),1);
[~,res13MixBaseSt,~,res13MixBaseSt_err] = nanMeanStdErr(res13MixBaseSt./sqrt(2),1);
[~,res23MixBaseSt,~,res23MixBaseSt_err] = nanMeanStdErr(res23MixBaseSt./sqrt(2),1);

res12MixBaseEn = squeeze(mixerBaselineEnd(1,:,:)-mixerBaselineEnd(2,:,:));
res13MixBaseEn = squeeze(mixerBaselineEnd(1,:,:)-mixerBaselineEnd(3,:,:));
res23MixBaseEn = squeeze(mixerBaselineEnd(2,:,:)-mixerBaselineEnd(3,:,:));
[~,res12MixBaseEn,~,res12MixBaseEn_err] = nanMeanStdErr(res12MixBaseEn./sqrt(2),1);
[~,res13MixBaseEn,~,res13MixBaseEn_err] = nanMeanStdErr(res13MixBaseEn./sqrt(2),1);
[~,res23MixBaseEn,~,res23MixBaseEn_err] = nanMeanStdErr(res23MixBaseEn./sqrt(2),1);

[meanRes12MixBaseSt,~,meanRes12MixBaseSt_err] = nanMeanStdErr(res12MixBaseSt);
[meanRes13MixBaseSt,~,meanRes13MixBaseSt_err] = nanMeanStdErr(res13MixBaseSt);
[meanRes23MixBaseSt,~,meanRes23MixBaseSt_err] = nanMeanStdErr(res23MixBaseSt);
[meanRes12MixBaseEn,~,meanRes12MixBaseEn_err] = nanMeanStdErr(res12MixBaseEn);
[meanRes13MixBaseEn,~,meanRes13MixBaseEn_err] = nanMeanStdErr(res13MixBaseEn);
[meanRes23MixBaseEn,~,meanRes23MixBaseEn_err] = nanMeanStdErr(res23MixBaseEn);

% phase
[meanSamplePhBaseSt,stdSamplePhBaseSt,meanSamplePhBaseSt_err,stdSamplePhBaseSt_err] = nanMeanStdErr(phaseBaselineStart,2);
[meanSamplePhBaseEn,stdSamplePhBaseEn,meanSamplePhBaseEn_err,stdSamplePhBaseEn_err] = nanMeanStdErr(phaseBaselineEnd,2);

meanBaselinePhaseStart = NaN(1,nMons);
stdBaselinePhaseStart = NaN(1,nMons);
meanBaselinePhaseEnd = NaN(1,nMons);
stdBaselinePhaseEnd = NaN(1,nMons);
meanBaselinePhaseStart_err = NaN(1,nMons);
stdBaselinePhaseStart_err = NaN(1,nMons);
meanBaselinePhaseEnd_err = NaN(1,nMons);
stdBaselinePhaseEnd_err = NaN(1,nMons);
for mon=1:nMons
    tmpSt = squeeze(phaseBaselineStart(mon,:,:));
    tmpEn = squeeze(phaseBaselineEnd(mon,:,:));
    
    [meanBaselinePhaseStart(mon),...
     stdBaselinePhaseStart(mon),...
     meanBaselinePhaseStart_err(mon),...
     stdBaselinePhaseStart_err(mon)] = nanMeanStdErr(tmpSt(:));
    
    [meanBaselinePhaseEnd(mon),...
     stdBaselinePhaseEnd(mon),...
     meanBaselinePhaseEnd_err(mon),...
     stdBaselinePhaseEnd_err(mon)] = nanMeanStdErr(tmpEn(:));
end

res12PhBaseSt = squeeze(phaseBaselineStart(1,:,:)-phaseBaselineStart(2,:,:));
res13PhBaseSt = squeeze(phaseBaselineStart(1,:,:)-phaseBaselineStart(3,:,:));
res23PhBaseSt = squeeze(phaseBaselineStart(2,:,:)-phaseBaselineStart(3,:,:));
[~,res12PhBaseSt,~,res12PhBaseSt_err] = nanMeanStdErr(res12PhBaseSt./sqrt(2),1);
[~,res13PhBaseSt,~,res13PhBaseSt_err] = nanMeanStdErr(res13PhBaseSt./sqrt(2),1);
[~,res23PhBaseSt,~,res23PhBaseSt_err] = nanMeanStdErr(res23PhBaseSt./sqrt(2),1);

res12PhBaseEn = squeeze(phaseBaselineEnd(1,:,:)-phaseBaselineEnd(2,:,:));
res13PhBaseEn = squeeze(phaseBaselineEnd(1,:,:)-phaseBaselineEnd(3,:,:));
res23PhBaseEn = squeeze(phaseBaselineEnd(2,:,:)-phaseBaselineEnd(3,:,:));
[~,res12PhBaseEn,~,res12PhBaseEn_err] = nanMeanStdErr(res12PhBaseEn./sqrt(2),1);
[~,res13PhBaseEn,~,res13PhBaseEn_err] = nanMeanStdErr(res13PhBaseEn./sqrt(2),1);
[~,res23PhBaseEn,~,res23PhBaseEn_err] = nanMeanStdErr(res23PhBaseEn./sqrt(2),1);

[meanRes12PhBaseSt,~,meanRes12PhBaseSt_err] = nanMeanStdErr(res12PhBaseSt);
[meanRes13PhBaseSt,~,meanRes13PhBaseSt_err] = nanMeanStdErr(res13PhBaseSt);
[meanRes23PhBaseSt,~,meanRes23PhBaseSt_err] = nanMeanStdErr(res23PhBaseSt);
[meanRes12PhBaseEn,~,meanRes12PhBaseEn_err] = nanMeanStdErr(res12PhBaseEn);
[meanRes13PhBaseEn,~,meanRes13PhBaseEn_err] = nanMeanStdErr(res13PhBaseEn);
[meanRes23PhBaseEn,~,meanRes23PhBaseEn_err] = nanMeanStdErr(res23PhBaseEn);


%% mean/std phases

% means/jitters per sample
[meanPhaseAlongPulse,...
 stdPhaseAlongPulse,...
 meanPhaseAlongPulse_err,...
 stdPhaseAlongPulse_err] = nanMeanStdErr(phases, 2);

[meanStdPhaseAlongPulse,~,meanStdPhaseAlongPulse_err,~] = nanMeanStdErr(stdPhaseAlongPulse(:,sampleRange),2); % mean sample jitter
[flatnessMeanPhaseAlongPulse,~,flatnessMeanPhaseAlongPulse_err,~] = nanMeanStdErr(meanPhaseAlongPulse(:,sampleRange),2); % flatness of mean phase along pulse

if (diodeIsPresent)
    [meanDiodeAlongPulse,...
     stdDiodeAlongPulse,...
     meanDiodeAlongPulse_err,...
     stdDiodeAlongPulse_err] = nanMeanStdErr(diodes, 2);

    [meanStdDiodeAlongPulse,~,meanStdDiodeAlongPulse_err,~] = nanMeanStdErr(stdDiodeAlongPulse(:,sampleRange),2); % mean sample jitter
end

[meanMixerAlongPulse,...
 stdMixerAlongPulse,...
 meanMixerAlongPulse_err,...
 stdMixerAlongPulse_err] = nanMeanStdErr(mixers, 2);

[meanStdMixerAlongPulse,~,meanStdMixerAlongPulse_err,~] = nanMeanStdErr(stdMixerAlongPulse(:,sampleRange),2); % mean sample jitter


% mean/jitters of pulse
[meanPulsePhase,...
 pulsePhaseFlatness,...
 meanPulsePhase_err,...
 pulsePhaseFlatness_err] = nanMeanStdErr(phases(:,:,sampleRange), 3);

[~,stdMeanPulsePhase,~,stdMeanPulsePhase_err] = nanMeanStdErr(meanPulsePhase,2); % jitter of mean phase
[meanPulsePhaseFlatness, stdPulsePhaseFlatness, meanPulsePhaseFlatness_err, stdPulsePhaseFlatness_err] = nanMeanStdErr(pulsePhaseFlatness,2); % mean/stds of all pulse flatness

if (diodeIsPresent)
    [meanDiode,...
     stdDiode,...
     meanDiode_err,...
     stdDiode_err] = nanMeanStdErr(diodes(:,:,sampleRange), 3);

    [~,stdMeanDiode,~,stdMeanDiode_err] = nanMeanStdErr(meanDiode,2); % jitter of mean diode
end

% Phase subtraction
if(strcmp(subtractType,'mean'))
    subtractPhase = nanmean(meanPulsePhase,2);
elseif(strcmp(subtractType,'none'))
    subtractPhase = zeros(1,nMons);
elseif(strcmp(subtractType,'value'))
    % keep the user given values
elseif(strcmp(subtractType,'file'))
    load(subtractFile,'subtractPhase');
end
for mon=1:nMons
    meanPhaseAlongPulse(mon,:) = meanPhaseAlongPulse(mon,:)-subtractPhase(mon);
    meanPulsePhase(mon,:) = meanPulsePhase(mon,:)-subtractPhase(mon);
end

% Calculate ratios of jitter between monitors (can be useful for comparing
% data sets etc.)
ratio12StdPhaseAlongPulse = stdPhaseAlongPulse(1,:)./stdPhaseAlongPulse(2,:);
ratio13StdPhaseAlongPulse = stdPhaseAlongPulse(1,:)./stdPhaseAlongPulse(3,:);
ratio23StdPhaseAlongPulse = stdPhaseAlongPulse(2,:)./stdPhaseAlongPulse(3,:);

ratio12StdMeanPulsePhase = stdMeanPulsePhase(1)./stdMeanPulsePhase(2);
ratio13StdMeanPulsePhase = stdMeanPulsePhase(1)./stdMeanPulsePhase(3);
ratio23StdMeanPulsePhase = stdMeanPulsePhase(2)./stdMeanPulsePhase(3);

% Flatness


% flag for NaN pulses (used for BPM processing)
isNaNPulse = isnan(meanPulsePhase);

%% correlations

% mean
[corrMeanMix1Mix2,corrMeanMix1Mix2_err] = nancorrcoef(meanPulsePhase(1,:),meanPulsePhase(2,:));
[corrMeanMix1Mix3,corrMeanMix1Mix3_err] = nancorrcoef(meanPulsePhase(1,:),meanPulsePhase(3,:));
[corrMeanMix2Mix3,corrMeanMix2Mix3_err] = nancorrcoef(meanPulsePhase(2,:),meanPulsePhase(3,:));

[fitMeanMix1Mix2,fitMeanMix1Mix2_rsq,fitMeanMix1Mix2_conf] = nanpolyfit(meanPulsePhase(1,:),meanPulsePhase(2,:),1);
[fitMeanMix1Mix3,fitMeanMix1Mix3_rsq,fitMeanMix1Mix3_conf] = nanpolyfit(meanPulsePhase(1,:),meanPulsePhase(3,:),1);
[fitMeanMix2Mix3,fitMeanMix2Mix3_rsq,fitMeanMix2Mix3_conf] = nanpolyfit(meanPulsePhase(2,:),meanPulsePhase(3,:),1);

% pulse shape (samples for one trigger)
corrShapeMix1Mix2 = NaN(1,nPulses);
corrShapeMix1Mix3 = NaN(1,nPulses);
corrShapeMix2Mix3 = NaN(1,nPulses);
corrShapeMix1Mix2_err = NaN(1,nPulses);
corrShapeMix1Mix3_err = NaN(1,nPulses);
corrShapeMix2Mix3_err = NaN(1,nPulses);

fitShapeMix1Mix2 = NaN(nPulses,2);
fitShapeMix1Mix3 = NaN(nPulses,2);
fitShapeMix2Mix3 = NaN(nPulses,2);
fitShapeMix1Mix2_rsq = NaN(1,nPulses);
fitShapeMix1Mix3_rsq = NaN(1,nPulses);
fitShapeMix2Mix3_rsq = NaN(nPulses,2,2);
fitShapeMix1Mix2_conf = NaN(nPulses,2,2);
fitShapeMix1Mix3_conf = NaN(nPulses,2,2);
fitShapeMix2Mix3_conf = NaN(nPulses,2,2);

for p=1:nPulses
    tmp1 = squeeze(phases(1,p,sampleRange));
    tmp2 = squeeze(phases(2,p,sampleRange));
    tmp3 = squeeze(phases(3,p,sampleRange));
    
    [corrShapeMix1Mix2(p),corrShapeMix1Mix2_err(p)] =...
        nancorrcoef(tmp1,tmp2);
       
    [corrShapeMix1Mix3(p),corrShapeMix1Mix3_err(p)] =...
        nancorrcoef(tmp1,tmp3);
    
    [corrShapeMix2Mix3(p),corrShapeMix2Mix3_err(p)] =...
        nancorrcoef(tmp2,tmp3);
    
    [fitShapeMix1Mix2(p,:),fitShapeMix1Mix2_rsq(p),fitShapeMix1Mix2_conf(p,:,:)] =...
        nanpolyfit(tmp1,tmp2,1);
    
    [fitShapeMix1Mix3(p,:),fitShapeMix1Mix3_rsq(p),fitShapeMix1Mix3_conf(p,:,:)] =...
        nanpolyfit(tmp1,tmp3,1);
    
    [fitShapeMix2Mix3(p,:),fitShapeMix2Mix3_rsq(p),fitShapeMix2Mix3_conf(p,:,:)] =...
        nanpolyfit(tmp2,tmp3,1);
    
end
[meanCorrShapeMix1Mix2,~,meanCorrShapeMix1Mix2_err,~] = nanMeanStdErr(corrShapeMix1Mix2);
[meanCorrShapeMix1Mix3,~,meanCorrShapeMix1Mix3_err,~] = nanMeanStdErr(corrShapeMix1Mix3);
[meanCorrShapeMix2Mix3,~,meanCorrShapeMix2Mix3_err,~] = nanMeanStdErr(corrShapeMix2Mix3);

[meanFitShapeMix1Mix2,~,meanFitShapeMix1Mix2_err,~] = nanMeanStdErr(fitShapeMix1Mix2(:,1));
[meanFitShapeMix1Mix3,~,meanFitShapeMix1Mix3_err,~] = nanMeanStdErr(fitShapeMix1Mix3(:,1));
[meanFitShapeMix2Mix3,~,meanFitShapeMix2Mix3_err,~] = nanMeanStdErr(fitShapeMix2Mix3(:,1));

% sample by sample (samples vs. time)
corrSamplesMix1Mix2 = NaN(1,nSamples);
corrSamplesMix1Mix3 = NaN(1,nSamples);
corrSamplesMix2Mix3 = NaN(1,nSamples);
corrSamplesMix1Mix2_err = NaN(1,nSamples);
corrSamplesMix1Mix3_err = NaN(1,nSamples);
corrSamplesMix2Mix3_err = NaN(1,nSamples);
for s=1:nSamples
    [corrSamplesMix1Mix2(s), corrSamplesMix1Mix2_err(s)] =...
        nancorrcoef(squeeze(phases(1,:,s)),squeeze(phases(2,:,s)));
    
    [corrSamplesMix1Mix3(s),corrSamplesMix1Mix3_err(s)] =...
        nancorrcoef(squeeze(phases(1,:,s)),squeeze(phases(3,:,s)));
    
    [corrSamplesMix2Mix3(s),corrSamplesMix2Mix3_err(s)] =...
        nancorrcoef(squeeze(phases(2,:,s)),squeeze(phases(3,:,s)));
end
[meanCorrSamplesMix1Mix2,~,meanCorrSamplesMix1Mix2_err,~] = nanMeanStdErr(corrSamplesMix1Mix2(sampleRange));
[meanCorrSamplesMix1Mix3,~,meanCorrSamplesMix1Mix3_err,~] = nanMeanStdErr(corrSamplesMix1Mix3(sampleRange));
[meanCorrSamplesMix2Mix3,~,meanCorrSamplesMix2Mix3_err,~] = nanMeanStdErr(corrSamplesMix2Mix3(sampleRange));

%% resolutions
diffPhase12 = squeeze(phases(1,:,:)-phases(2,:,:));
diffPhase13 = squeeze(phases(1,:,:)-phases(3,:,:));
diffPhase23 = squeeze(phases(2,:,:)-phases(3,:,:));

[~,resolution12,~,resolution12_err] = nanMeanStdErr(diffPhase12./sqrt(2),1);
[~,resolution13,~,resolution13_err] = nanMeanStdErr(diffPhase13./sqrt(2),1);
[~,resolution23,~,resolution23_err] = nanMeanStdErr(diffPhase23./sqrt(2),1);

[meanResolution12,~,meanResolution12_err,~] = nanMeanStdErr(resolution12(sampleRange));
[meanResolution13,~,meanResolution13_err,~] = nanMeanStdErr(resolution13(sampleRange));
[meanResolution23,~,meanResolution23_err,~] = nanMeanStdErr(resolution23(sampleRange));
