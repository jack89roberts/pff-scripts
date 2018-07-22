%% PETS
% Currently uses I as the "diode like" signal.

%% signal processing

% extract
try
    petsI = double(extractCTFSignalFromMergedData('CE_SCOPE03_CH01.Acquisition.value.value',CTFData));
    petsQ = double(extractCTFSignalFromMergedData('CE_SCOPE03_CH02.Acquisition.value.value',CTFData));
    petsISensitivity = double(extractCTFSignalFromMergedData('CE_SCOPE03_CH01.Acquisition.sensitivity.value',CTFData));
    petsQSensitivity = double(extractCTFSignalFromMergedData('CE_SCOPE03_CH02.Acquisition.sensitivity.value',CTFData));
    petsISensitivity = petsISensitivity(1);
    petsQSensitivity = petsQSensitivity(1);
    petsI = petsI.*petsISensitivity;
    petsQ = petsQ.*petsQSensitivity;
    
catch
    fprintf('Error extracting PETS data! PETS analysis not completed.\n');
    return;
end

petsSampInterval = extractCTFSignalFromMergedData('CE_SCOPE03_CH01.Acquisition.sampleInterval.value',CTFData);
petsSampInterval = petsSampInterval(1);

% Some data had scope setup error with different number of points in I and
% Q. Use cubic spline to resample the data (should check effect on
% resolution but by eye results look sensible.
[~,nSamplesI] = size(petsI);
[~,nSamplesQ] = size(petsQ);
if (nSamplesI ~= nSamplesQ)
    fprintf('Number of points in PETS I and PETS Q channels is different. Re-smapling data.\n')
    
    % resample to size of signal with fewer of points
    if (nSamplesI>nSamplesQ)
        petsI = spline(1:nSamplesI,petsI,linspace(1,nSamplesI,nSamplesQ));
        petsSampInterval = extractCTFSignalFromMergedData('CE_SCOPE03_CH02.Acquisition.sampleInterval.value',CTFData);
        petsSampInterval = petsSampInterval(1);
    else
        petsQ = spline(1:nSamplesQ,petsQ,linspace(1,nSamplesQ,nSamplesI));
        petsSampInterval = extractCTFSignalFromMergedData('CE_SCOPE03_CH01.Acquisition.sampleInterval.value',CTFData);
        petsSampInterval = petsSampInterval(1);

    end
end

% calculate phases
[petsPhase,petsPower] = getPhaseIQ(petsI, petsQ);

% remove duplicate pulses (helps to eradicate any acquisition issues - e.g.
% too slow acquisition, acquisition all zeroes etc.)
if (stripDuplicates)
    petsI = removeDuplicatePulses(petsI);
    petsQ = removeDuplicatePulses(petsQ);
    petsPower = removeDuplicatePulses(petsPower);
    petsPhase = removeDuplicatePulses(petsPhase);
end

% align
% uses petsI to align. Maybe could be improved by checking which between
% petsI and petsQ is flatter and using that.
[petsI,others,petsPulseSampleRange] = getAlignedXCorr(petsI,alignTo,{petsQ,petsPhase,petsPower});
petsQ = others{1};
petsPhase = others{2};
petsPower = others{3};

% make time axis
[~,petsNSamples] = size(petsI);
petsTimeAxis = petsSampInterval.*(0:petsNSamples-1);
petsTimeAxis = petsTimeAxis-petsTimeAxis(petsPulseSampleRange(1));
petsPulseLength = petsSampInterval.*(length(petsPulseSampleRange)-1);  % ns

% choose sample range
if (~isempty(petsPulseSampleRange))
    tmpStartSample = petsPulseSampleRange(1) + round(startDelayCB./petsSampInterval);
    tmpEndSample = petsPulseSampleRange(end) - round(endDelayCB./petsSampInterval);
else
    tmpStartSample = 1;
    tmpEndSample = petsNSamples;
end
petsSampleRange = tmpStartSample:tmpEndSample;

% remove duplicate pulses
petsI = removeDuplicatePulses(petsI);
petsQ = removeDuplicatePulses(petsQ);
petsPhase = removeDuplicatePulses(petsPhase);
petsPower  = removeDuplicatePulses(petsPower);

% remove drift
if (stripDrift)
    petsPhase = removeDrift(petsPhase,petsSampleRange,driftNAvg);
end

% remove bad pulses
if (stripOutliers)
    petsPhase = removeBadPulses(petsPhase,petsSampleRange);
end

%% mean/std phases

% means/jitters per sample
[petsMeanPhaseAlongPulse,...
 petsStdPhaseAlongPulse,...
 petsMeanPhaseAlongPulse_err,...
 petsStdPhaseAlongPulse_err] = nanMeanStdErr(petsPhase);

[petsMeanStdPhaseAlongPulse,~,petsMeanStdPhaseAlongPulse_err,~] = nanMeanStdErr(petsStdPhaseAlongPulse(petsSampleRange)); % mean sample jitter
[petsFlatnessMeanPhaseAlongPulse,~,petsFlatnessMeanPhaseAlongPulse_err,~] = nanMeanStdErr(petsMeanPhaseAlongPulse(petsSampleRange)); % flatness of mean phase along pulse

[petsMeanPowerAlongPulse,...
 petsStdPowerAlongPulse,...
 petsMeanPowerAlongPulse_err,...
 petsStdPowerAlongPulse_err] = nanMeanStdErr(petsPower);

[petsMeanStdPowerAlongPulse,~,petsMeanStdPowerAlongPulse_err,~] = nanMeanStdErr(petsStdPowerAlongPulse(petsSampleRange)); % mean sample jitter

[petsIMeanAlongPulse,...
 petsIStdAlongPulse,...
 petsIMeanAlongPulse_err,...
 petsIStdAlongPulse_err] = nanMeanStdErr(petsI);

[petsIMeanStdAlongPulse,~,petsIMeanStdAlongPulse_err,~] = nanMeanStdErr(petsIStdAlongPulse(petsSampleRange)); % mean sample jitter

[petsQMeanAlongPulse,...
 petsQStdAlongPulse,...
 petsQMeanAlongPulse_err,...
 petsQStdAlongPulse_err] = nanMeanStdErr(petsQ);

[petsQMeanStdAlongPulse,~,petsQMeanStdAlongPulse_err,~] = nanMeanStdErr(petsQStdAlongPulse(petsSampleRange)); % mean sample jitter

% mean/jitters of pulse
[petsMeanPulsePhase,...
 petsPulsePhaseFlatness,...
 petsMeanPulsePhase_err,...
 petsPulsePhaseFlatness_err] = nanMeanStdErr(petsPhase(:,petsSampleRange),2);

[~,petsStdMeanPulsePhase,~,petsStdMeanPulsePhase_err] = nanMeanStdErr(petsMeanPulsePhase); % jitter of mean phase
[petsMeanPulsePhaseFlatness, petsStdPulsePhaseFlatness, petsMeanPulsePhaseFlatness_err, petsStdPulsePhaseFlatness_err] = nanMeanStdErr(petsPulsePhaseFlatness); % mean/stds of all pulse flatness

[petsMeanPower,...
 petsPowerFlatness,...
 petsMeanPower_err,...
 petsStdPower_err] = nanMeanStdErr(petsPower(:,petsSampleRange), 2);

[~,petsStdMeanDiode,~,petsStdMeanDiode_err] = nanMeanStdErr(petsMeanPower); % jitter of mean phase
[petsMeanPowerFlatness, petsStdPowerFlatness, petsMeanPowerFlatness_err, petsStdPowerFlatness_err] = nanMeanStdErr(petsPowerFlatness); % mean/stds of all pulse flatness

% Phase subtraction
if(strcmp(subtractType,'mean'))
    petsSubtractPhase = nanmean(petsMeanPulsePhase);
elseif(strcmp(subtractType,'none'))
    petsSubtractPhase = 0;
elseif(strcmp(subtractType,'value'))
    % keep the user given values
elseif(strcmp(subtractType,'file'))
    load(subtractFile,'petsSubtractPhase');
end
petsMeanPhaseAlongPulse = petsMeanPhaseAlongPulse-petsSubtractPhase;
petsMeanPulsePhase = petsMeanPulsePhase-petsSubtractPhase;

% flag for NaN pulses (used for BPM processing)
petsIsNaNPulse = isnan(petsMeanPulsePhase);

% correlation on mean
[corrMeanMix3PETS,corrMeanMix3PETS_err] = nancorrcoef(meanPulsePhase(3,:),petsMeanPulsePhase);
[fitMeanMix3PETS,fitMeanMix3PETS_rsq,fitMeanMix3PETS_conf] = nanpolyfit(meanPulsePhase(3,:),petsMeanPulsePhase,1);

[corrMeanMix2PETS,corrMeanMix2PETS_err] = nancorrcoef(meanPulsePhase(2,:),petsMeanPulsePhase);
[fitMeanMix2PETS,fitMeanMix2PETS_rsq,fitMeanMix2PETS_conf] = nanpolyfit(meanPulsePhase(2,:),petsMeanPulsePhase,1);

[corrMeanMix1PETS,corrMeanMix1PETS_err] = nancorrcoef(meanPulsePhase(1,:),petsMeanPulsePhase);
[fitMeanMix1PETS,fitMeanMix1PETS_rsq,fitMeanMix1PETS_conf] = nanpolyfit(meanPulsePhase(1,:),petsMeanPulsePhase,1);

% To consider implementing - correlations vs. sample number.
% Difficult due to different sampling rates.
% pulse shape (samples for one trigger)
% sample by sample (samples vs. time)
% resolution
