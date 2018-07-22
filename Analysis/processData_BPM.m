%% processData_Frascati.m
% All signal processing, calculations related to the BPMs.
% This script should only be called from processData.m, not used by itself.
% It requires all inputs from processData.m and all outputs from
% processData_Frascati.m to work.

% June 2015, Jack Roberts
% TO ADD: Fit gradients etc.
%% extract signals

% some parameters, shouldn't need to change them
bpmCCTimePerSample = 5.208;
bpmCTTimePerSample = 10.417;
devsPath = '../devices/';
bpmDevFiles = {'bpmCL.devs', 'bpmCT.devs', 'bpmCR.devs', 'bpmCC.devs', 'bpmCB.devs'};

bpmNames = {};
for i=1:length(bpmDevFiles)
    tmpBPMFile = fopen([devsPath bpmDevFiles{i}]);
    tmpBPMNames = textscan(tmpBPMFile,'%s');
    fclose(tmpBPMFile);

    for j=1:length(tmpBPMNames)
        bpmNames = [bpmNames; tmpBPMNames{j}];
    end
end
bpmNames = strrep(bpmNames,'.','_');
nBPMs = length(bpmNames);

bpmS = cell(1,nBPMs);
bpmH = cell(1,nBPMs);
bpmV = cell(1,nBPMs);
useBPM = zeros(1,nBPMs);
bpmNSamples = NaN(1,nBPMs);
useVertical = 0;
for i=1:nBPMs
    try
        [bpmH{i},bpmS{i},bpmV{i}] = extractBPMFromCTFData(bpmNames{i},CTFData);
        [~,bpmNSamples(i)] = size(bpmH{i});

        if (bpmNSamples(i)>100)
            useBPM(i) = 1;
        else
            continue;
        end

        if (sum(isnan(bpmV{i})) ~= bpmNSamples(i)) % if at least 1 vertical bpm signal has non nan data, include vertical analysis
            useVertical = 1;
        end
        
    catch
        % some error in extracting this BPM, so don't use it.
    end
end

refBPM1Index = find(ismember(bpmNames,refBPM1)); 
refBPM2Index = find(ismember(bpmNames,refBPM2));

if (isempty(refBPM1Index))
    useRefBPM1 = 0;
    warning('refBPM1 %s not found',refBPM1);
else
    useRefBPM1 = 1;
end

if (isempty(refBPM2Index))
    useRefBPM2 = 0;
    warning('refBPM2 %s not found',refBPM2);
else
    useRefBPM2 = 1;
end

%% remove (NaN) any duplicate pulses

if (stripDuplicates)
    for i=1:nBPMs
        if (useBPM(i))
            bpmS{i} = removeDuplicatePulses(bpmS{i});
            bpmH{i} = removeDuplicatePulses(bpmH{i});
            if (useVertical)
                bpmV{i} = removeDuplicatePulses(bpmV{i});
            end
        end
    end
end

%% align signals
bpmPulseRange = cell(1,nBPMs);
bpmTimePerSample = NaN(1,nBPMs);
for i=1:nBPMs
    if (useBPM(i))
        tmpBPMS = bpmS{i};
        tmpBPMH = bpmH{i};
        tmpBPMV = bpmV{i};

        if (bpmNSamples(i) > 400)
            tmpTimePerSample = bpmCCTimePerSample;
        else
            tmpTimePerSample = bpmCTTimePerSample;
        end

        if (useVertical)
            [tmpBPMS,tmpOthers,tmpPulseRange,errMessage] = getAlignedXCorr(tmpBPMS,alignTo,{tmpBPMH, tmpBPMV});
            tmpBPMH = tmpOthers{1};
            tmpBPMV = tmpOthers{2};
        else
            [tmpBPMS,tmpOthers,tmpPulseRange,errMessage] = getAlignedXCorr(tmpBPMS,alignTo,{tmpBPMH});
            tmpBPMH = tmpOthers{1};
        end

        bpmS{i} = tmpBPMS;
        bpmH{i} = tmpBPMH;
        if (useVertical)
            bpmV{i} = tmpBPMV;
        end

        bpmPulseRange{i} = tmpPulseRange;
        bpmTimePerSample(i) = tmpTimePerSample;

        if (~isempty(errMessage))
            fprintf('%s: %s\n',bpmNames{i},errMessage);
        end
    end
end
%% find sample ranges to use to calculate means etc. (pulse starts at nRMS)

bpmSampleRange = cell(1,nBPMs);
for i=1:nBPMs
    if (useBPM(i))
        if (strcmp(bpmNames{i}(1:2),'CC') || strcmp(bpmNames{i}(1:2),'CB')) % assume length of pulse in TBL frascati monitor for bpms in TL2 and TBL
            tmpStartDelay = startDelayCB;
            tmpEndDelay = endDelayCB;
        else % assume length of upstream frascati phase monitors for bpms in and before combiner ring
            tmpStartDelay = startDelayCT;
            tmpEndDelay = endDelayCT;
        end

        if (~isempty(bpmPulseRange{i}))
            tmpStartSample = bpmPulseRange{i}(1) + round(tmpStartDelay./bpmTimePerSample(i));
            tmpEndSample = bpmPulseRange{i}(end) - round(tmpEndDelay./bpmTimePerSample(i));
        else
            tmpStartSample = 1;
            tmpEndSample = bpmNSamples(i);
        end

        bpmSampleRange{i} = tmpStartSample:tmpEndSample;
    end
end

%% remove bad pulses - BPM S
if (stripOutliers) 
    for i=1:nBPMs
        if (useBPM(i))
            tmpRange = bpmSampleRange{i};

            % remove pulses with bad transmission (S) in H and V as well
            if (useVertical)
                [bpmS{i},others] = removeBadPulses(bpmS{i},tmpRange,{bpmH{i},bpmV{i}});
                bpmH{i} = others{1};
                bpmV{i} = others{2};
            else
                [bpmS{i},others] = removeBadPulses(bpmS{i},tmpRange,{bpmH{i},bpmV{i}});
                bpmH{i} = others{1};
            end

            % then remove outliers separately in H and V
            bpmH{i} = removeBadPulses(bpmH{i},tmpRange);
            if (useVertical)
                bpmV{i} = removeBadPulses(bpmV{i},tmpRange);
            end
        end
    end
end

%% remove drifts in BPM H and BPMV 
if (stripDrift) 
    for i=1:nBPMs
        if (useBPM(i))
            tmpRange = bpmSampleRange{i};

            bpmH{i} = removeDrift(bpmH{i},tmpRange,driftNAvg);
            if (useVertical)
                bpmV{i} = removeDrift(bpmV{i},tmpRange,driftNAvg);
            end
        end
    end
end

%% remove bad pulses - BPM H and BPM V
if (stripOutliers)
    for i=1:nBPMs
        if (useBPM(i))
            tmpRange = bpmSampleRange{i};

            bpmH{i} = removeBadPulses(bpmH{i},tmpRange);
            if (useVertical)
                bpmV{i} = removeBadPulses(bpmV{i},tmpRange);
            end
        end
    end
end

%% calculate means

% Pre-define arrays
meanBPMHAlongPulse = cell(1,nBPMs);
stdBPMHAlongPulse = cell(1,nBPMs);
meanBPMHAlongPulse_err = cell(1,nBPMs);
stdBPMHAlongPulse_err = cell(1,nBPMs);
meanStdBPMHAlongPulse = cell(1,nBPMs);
meanStdBPMHAlongPulse_err = cell(1,nBPMs);
meanBPMH = cell(1,nBPMs);
stdPulseBPMH = cell(1,nBPMs);
meanBPMH_err = cell(1,nBPMs);
stdPulseBPMH_err = cell(1,nBPMs);
stdMeanBPMH = cell(1,nBPMs);
stdMeanBPMH_err = cell(1,nBPMs);
meanBPMSAlongPulse = cell(1,nBPMs);
stdBPMSAlongPulse = cell(1,nBPMs);
meanBPMSAlongPulse_err = cell(1,nBPMs);
stdBPMSAlongPulse_err = cell(1,nBPMs);
meanStdBPMSAlongPulse = cell(1,nBPMs);
meanStdBPMSAlongPulse_err = cell(1,nBPMs);
meanBPMS = cell(1,nBPMs);
stdPulseBPMS = cell(1,nBPMs);
meanBPMS_err = cell(1,nBPMs);
stdPulseBPMS_err = cell(1,nBPMs);
stdMeanBPMS = cell(1,nBPMs);
stdMeanBPMS_err = cell(1,nBPMs);
if (useVertical)
    meanBPMVAlongPulse = cell(1,nBPMs);
    stdBPMVAlongPulse = cell(1,nBPMs);
    meanBPMVAlongPulse_err = cell(1,nBPMs);
    stdBPMVAlongPulse_err = cell(1,nBPMs);
    meanStdBPMVAlongPulse = cell(1,nBPMs);
    meanStdBPMVAlongPulse_err = cell(1,nBPMs);
    meanBPMV = cell(1,nBPMs);
    stdPulseBPMV = cell(1,nBPMs);
    meanBPMV_err = cell(1,nBPMs);
    stdPulseBPMV_err = cell(1,nBPMs);
    stdMeanBPMV = cell(1,nBPMs);
    stdMeanBPMV_err = cell(1,nBPMs);
end

% calculate means, errors
for i=1:nBPMs
    if (useBPM(i))
        [meanBPMHAlongPulse{i},stdBPMHAlongPulse{i},meanBPMHAlongPulse_err{i},stdBPMHAlongPulse_err{i}] = nanMeanStdErr(bpmH{i});
        [meanStdBPMHAlongPulse{i},~,meanStdBPMHAlongPulse_err{i},~] = nanMeanStdErr(stdBPMHAlongPulse{i}(bpmSampleRange{i}));
        [meanBPMH{i},stdPulseBPMH{i},meanBPMH_err{i},stdPulseBPMH_err{i}] = nanMeanStdErr(bpmH{i}(:,bpmSampleRange{i}),2);
        [~,stdMeanBPMH{i},~,stdMeanBPMH_err{i}] = nanMeanStdErr(meanBPMH{i});

        [meanBPMSAlongPulse{i},stdBPMSAlongPulse{i},meanBPMSAlongPulse_err{i},stdBPMSAlongPulse_err{i}] = nanMeanStdErr(bpmS{i});
        [meanStdBPMSAlongPulse{i},~,meanStdBPMSAlongPulse_err{i},~] = nanMeanStdErr(stdBPMSAlongPulse{i}(bpmSampleRange{i}));
        [meanBPMS{i},stdPulseBPMS{i},meanBPMS_err{i},stdPulseBPMS_err{i}] = nanMeanStdErr(bpmS{i}(:,bpmSampleRange{i}),2);
        [~,stdMeanBPMS{i},~,stdMeanBPMS_err{i}] = nanMeanStdErr(meanBPMS{i});

        if (useVertical)
            [meanBPMVAlongPulse{i},stdBPMVAlongPulse{i},meanBPMVAlongPulse_err{i},stdBPMVAlongPulse_err{i}] = nanMeanStdErr(bpmV{i});
            [meanStdBPMVAlongPulse{i},~,meanStdBPMVAlongPulse_err{i},~] = nanMeanStdErr(stdBPMVAlongPulse{i}(bpmSampleRange{i}));
            [meanBPMV{i},stdPulseBPMV{i},meanBPMV_err{i},stdPulseBPMV_err{i}] = nanMeanStdErr(bpmV{i}(:,bpmSampleRange{i}),2);
            [~,stdMeanBPMV{i},~,stdMeanBPMV_err{i}] = nanMeanStdErr(meanBPMV{i});
        end
    end
end

if (useRefBPM1)
    meanRefBPM1_S = meanBPMS{refBPM1Index};
    meanRefBPM1_H = meanBPMH{refBPM1Index};
    if (useVertical)
        meanRefBPM1_V = meanBPMV{refBPM1Index};
    end
end
if (useRefBPM2)
    meanRefBPM2_S = meanBPMS{refBPM2Index};
    meanRefBPM2_H = meanBPMH{refBPM2Index};
    if (useVertical)
        meanRefBPM2_V = meanBPMV{refBPM2Index};
    end
end


%% calculate correlations

corrMon3_BPMS = NaN(1,nBPMs);
corrMon3_BPMH = NaN(1,nBPMs);
corrMon3_BPMV = NaN(1,nBPMs);
corrMon2_BPMS = NaN(1,nBPMs);
corrMon2_BPMH = NaN(1,nBPMs);
corrMon2_BPMV = NaN(1,nBPMs);
corrMon1_BPMS = NaN(1,nBPMs);
corrMon1_BPMH = NaN(1,nBPMs);
corrMon1_BPMV = NaN(1,nBPMs);
corrMon3_BPMS_err = NaN(1,nBPMs);
corrMon3_BPMH_err = NaN(1,nBPMs);
corrMon3_BPMV_err = NaN(1,nBPMs);
corrMon2_BPMS_err = NaN(1,nBPMs);
corrMon2_BPMH_err = NaN(1,nBPMs);
corrMon2_BPMV_err = NaN(1,nBPMs);
corrMon1_BPMS_err = NaN(1,nBPMs);
corrMon1_BPMH_err = NaN(1,nBPMs);
corrMon1_BPMV_err = NaN(1,nBPMs);

corrRefBPM1H_BPMH = NaN(1,nBPMs);
corrRefBPM1H_BPMS = NaN(1,nBPMs);
corrRefBPM1H_BPMV = NaN(1,nBPMs);
corrRefBPM1V_BPMH = NaN(1,nBPMs);
corrRefBPM1V_BPMS = NaN(1,nBPMs);
corrRefBPM1V_BPMV = NaN(1,nBPMs);
corrRefBPM1S_BPMH = NaN(1,nBPMs);
corrRefBPM1S_BPMS = NaN(1,nBPMs);
corrRefBPM1S_BPMV = NaN(1,nBPMs);
corrRefBPM2H_BPMH = NaN(1,nBPMs);
corrRefBPM2H_BPMS = NaN(1,nBPMs);
corrRefBPM2H_BPMV = NaN(1,nBPMs);
corrRefBPM2V_BPMH = NaN(1,nBPMs);
corrRefBPM2V_BPMS = NaN(1,nBPMs);
corrRefBPM2V_BPMV = NaN(1,nBPMs);
corrRefBPM2S_BPMH = NaN(1,nBPMs);
corrRefBPM2S_BPMS = NaN(1,nBPMs);
corrRefBPM2S_BPMV = NaN(1,nBPMs);
corrRefBPM1H_BPMH_err = NaN(1,nBPMs);
corrRefBPM1H_BPMS_err = NaN(1,nBPMs);
corrRefBPM1H_BPMV_err = NaN(1,nBPMs);
corrRefBPM1V_BPMH_err = NaN(1,nBPMs);
corrRefBPM1V_BPMS_err = NaN(1,nBPMs);
corrRefBPM1V_BPMV_err = NaN(1,nBPMs);
corrRefBPM1S_BPMH_err = NaN(1,nBPMs);
corrRefBPM1S_BPMS_err = NaN(1,nBPMs);
corrRefBPM1S_BPMV_err = NaN(1,nBPMs);
corrRefBPM2H_BPMH_err = NaN(1,nBPMs);
corrRefBPM2H_BPMS_err = NaN(1,nBPMs);
corrRefBPM2H_BPMV_err = NaN(1,nBPMs);
corrRefBPM2V_BPMH_err = NaN(1,nBPMs);
corrRefBPM2V_BPMS_err = NaN(1,nBPMs);
corrRefBPM2V_BPMV_err = NaN(1,nBPMs);
corrRefBPM2S_BPMH_err = NaN(1,nBPMs);
corrRefBPM2S_BPMS_err = NaN(1,nBPMs);
corrRefBPM2S_BPMV_err = NaN(1,nBPMs);

for i=1:nBPMs
    if (useBPM(i))
        %%%%%%%%%%%%%%%%%CORRELATIONS BETWEEN PHASE MONITORS AND ALL BPMS%%%%%%%%%%%%%%%%%%%%%%%%
        % Mon3 - BPMH
        [corrMon3_BPMH(i),corrMon3_BPMH_err(i)] = nancorrcoef(meanPulsePhase(3,:),meanBPMH{i});

        % Mon3 - BPMS
        [corrMon3_BPMS(i),corrMon3_BPMS_err(i)] = nancorrcoef(meanPulsePhase(3,:),meanBPMS{i});

        % Mon2 - BPMH
        [corrMon2_BPMH(i),corrMon2_BPMH_err(i)] = nancorrcoef(meanPulsePhase(2,:),meanBPMH{i});

        % Mon2 - BPMS
        [corrMon2_BPMS(i),corrMon2_BPMS_err(i)] = nancorrcoef(meanPulsePhase(2,:),meanBPMS{i});

        % Mon1 - BPMH
        [corrMon1_BPMH(i),corrMon1_BPMH_err(i)] = nancorrcoef(meanPulsePhase(1,:),meanBPMH{i});

        % Mon1 - BPMS
        [corrMon1_BPMS(i),corrMon1_BPMS_err(i)] = nancorrcoef(meanPulsePhase(1,:),meanBPMS{i});

        if (useVertical)

            % Mon3 - BPMV
            [corrMon3_BPMV(i),corrMon3_BPMV_err(i)] = nancorrcoef(meanPulsePhase(3,:),meanBPMV{i});

            % Mon2 - BPMV
            [corrMon2_BPMV(i),corrMon2_BPMV_err(i)] = nancorrcoef(meanPulsePhase(2,:),meanBPMV{i});

            % Mon1 - BPMV
            [corrMon1_BPMV(i),corrMon1_BPMV_err(i)] = nancorrcoef(meanPulsePhase(1,:),meanBPMV{i});
        end

        %%%%%%%%%%%%%%%%%CORRELATIONS BETWEEN REFBPM2 AND ALL OTHER BPMS%%%%%%%%%%%%%%%%%%%%%%%%
        if (useRefBPM2)
            % RefBPM2H - BPMH
            [corrRefBPM2H_BPMH(i),corrRefBPM2H_BPMH_err(i)] = nancorrcoef(meanRefBPM2_H,meanBPMH{i});

            % RefBPM2H - BPMS
            [corrRefBPM2H_BPMS(i),corrRefBPM2H_BPMS_err(i)] = nancorrcoef(meanRefBPM2_H,meanBPMS{i});

             % RefBPM2S - BPMH
            [corrRefBPM2S_BPMH(i),corrRefBPM2S_BPMH_err(i)] = nancorrcoef(meanRefBPM2_S,meanBPMH{i});

            % RefBPM2S - BPMS
            [corrRefBPM2S_BPMS(i),corrRefBPM2S_BPMS_err(i)] = nancorrcoef(meanRefBPM2_S,meanBPMS{i});


            if (useVertical)

                % RefBPM2H - BPMV
                [corrRefBPM2H_BPMV(i),corrRefBPM2H_BPMV_err(i)] = nancorrcoef(meanRefBPM2_H,meanBPMV{i});

                % RefBPM2V - BPMH
                [corrRefBPM2V_BPMH(i),corrRefBPM2V_BPMH_err(i)] = nancorrcoef(meanRefBPM2_V,meanBPMH{i});

                % RefBPM2V - BPMS
                [corrRefBPM2V_BPMS(i),corrRefBPM2V_BPMS_err(i)] = nancorrcoef(meanRefBPM2_V,meanBPMS{i});

                % RefBPM2V - BPMV
                [corrRefBPM2V_BPMV(i),corrRefBPM2V_BPMV_err(i)] = nancorrcoef(meanRefBPM2_V,meanBPMV{i});

                % RefBPM2S - BPMV
                [corrRefBPM2S_BPMV(i),corrRefBPM2S_BPMV_err(i)] = nancorrcoef(meanRefBPM2_S,meanBPMV{i});

            end

        end

        %%%%%%%%%%%%%%%%%CORRELATIONS BETWEEN REFBPM1 AND ALL OTHER BPMS%%%%%%%%%%%%%%%%%%%%%%%%
        if (useRefBPM1)
            % RefBPM1H - BPMH
            [corrRefBPM1H_BPMH(i),corrRefBPM1H_BPMH_err(i)] = nancorrcoef(meanRefBPM1_H,meanBPMH{i});

            % RefBPM1H - BPMS
            [corrRefBPM1H_BPMS(i),corrRefBPM1H_BPMS_err(i)] = nancorrcoef(meanRefBPM1_H,meanBPMS{i});

            % RefBPM1S - BPMH
            [corrRefBPM1S_BPMH(i),corrRefBPM1S_BPMH_err(i)] = nancorrcoef(meanRefBPM1_S,meanBPMH{i});

            % RefBPM1S - BPMS
            [corrRefBPM1S_BPMS(i),corrRefBPM1S_BPMS_err(i)] = nancorrcoef(meanRefBPM1_S,meanBPMS{i});

            if (useVertical)
                % RefBPM1H - BPMV
                [corrRefBPM1H_BPMV(i),corrRefBPM1H_BPMV_err(i)] = nancorrcoef(meanRefBPM1_H,meanBPMV{i});

                % RefBPM1V - BPMH
                [corrRefBPM1V_BPMH(i),corrRefBPM1V_BPMH_err(i)] = nancorrcoef(meanRefBPM1_V,meanBPMH{i});

                % RefBPM1V - BPMS
                [corrRefBPM1V_BPMS(i),corrRefBPM1V_BPMS_err(i)] = nancorrcoef(meanRefBPM1_V,meanBPMS{i});

                % RefBPM1V - BPMV
                [corrRefBPM1V_BPMV(i),corrRefBPM1V_BPMV_err(i)] = nancorrcoef(meanRefBPM1_V,meanBPMV{i});

                % RefBPM1S - BPMV
                [corrRefBPM1S_BPMV(i),corrRefBPM1S_BPMV_err(i)] = nancorrcoef(meanRefBPM1_S,meanBPMV{i});

            end

        end
    end
end

%% Correlation between phase correlation and bpm correlation

if (useRefBPM1)
    
    [bpmH_CORR_Mon3_RefBPM1H,bpmH_CORR_Mon3_RefBPM1H_err] = nancorrcoef(corrMon3_BPMH,corrRefBPM1H_BPMH);
    [bpmS_CORR_Mon3_RefBPM1H,bpmS_CORR_Mon3_RefBPM1H_err] = nancorrcoef(corrMon3_BPMS,corrRefBPM1H_BPMS);
    [bpmH_CORR_Mon3_RefBPM1S,bpmH_CORR_Mon3_RefBPM1S_err] = nancorrcoef(corrMon3_BPMH,corrRefBPM1S_BPMH);
    [bpmS_CORR_Mon3_RefBPM1S,bpmS_CORR_Mon3_RefBPM1S_err] = nancorrcoef(corrMon3_BPMS,corrRefBPM1S_BPMS);
    
    [bpmH_CORR_Mon2_RefBPM1H,bpmH_CORR_Mon2_RefBPM1H_err] = nancorrcoef(corrMon2_BPMH,corrRefBPM1H_BPMH);
    [bpmS_CORR_Mon2_RefBPM1H,bpmS_CORR_Mon2_RefBPM1H_err] = nancorrcoef(corrMon2_BPMS,corrRefBPM1H_BPMS);
    [bpmH_CORR_Mon2_RefBPM1S,bpmH_CORR_Mon2_RefBPM1S_err] = nancorrcoef(corrMon2_BPMH,corrRefBPM1S_BPMH);
    [bpmS_CORR_Mon2_RefBPM1S,bpmS_CORR_Mon2_RefBPM1S_err] = nancorrcoef(corrMon2_BPMS,corrRefBPM1S_BPMS);
    
    [bpmH_CORR_Mon1_RefBPM1H,bpmH_CORR_Mon1_RefBPM1H_err] = nancorrcoef(corrMon1_BPMH,corrRefBPM1H_BPMH);
    [bpmS_CORR_Mon1_RefBPM1H,bpmS_CORR_Mon1_RefBPM1H_err] = nancorrcoef(corrMon1_BPMS,corrRefBPM1H_BPMS);
    [bpmH_CORR_Mon1_RefBPM1S,bpmH_CORR_Mon1_RefBPM1S_err] = nancorrcoef(corrMon1_BPMH,corrRefBPM1S_BPMH);
    [bpmS_CORR_Mon1_RefBPM1S,bpmS_CORR_Mon1_RefBPM1S_err] = nancorrcoef(corrMon1_BPMS,corrRefBPM1S_BPMS);
    
    
    if (useVertical)
        [bpmH_CORR_Mon3_RefBPM1V,bpmH_CORR_Mon3_RefBPM1V_err] = nancorrcoef(corrMon3_BPMH,corrRefBPM1V_BPMH);
        [bpmS_CORR_Mon3_RefBPM1V,bpmS_CORR_Mon3_RefBPM1V_err] = nancorrcoef(corrMon3_BPMS,corrRefBPM1V_BPMS);
        [bpmV_CORR_Mon3_RefBPM1V,bpmV_CORR_Mon3_RefBPM1V_err] = nancorrcoef(corrMon3_BPMV,corrRefBPM1V_BPMV);
        [bpmV_CORR_Mon3_RefBPM1H,bpmV_CORR_Mon3_RefBPM1H_err] = nancorrcoef(corrMon3_BPMV,corrRefBPM1H_BPMV);
        [bpmV_CORR_Mon3_RefBPM1S,bpmV_CORR_Mon3_RefBPM1S_err] = nancorrcoef(corrMon3_BPMV,corrRefBPM1S_BPMV);
        
        [bpmH_CORR_Mon2_RefBPM1V,bpmH_CORR_Mon2_RefBPM1V_err] = nancorrcoef(corrMon2_BPMH,corrRefBPM1V_BPMH);
        [bpmS_CORR_Mon2_RefBPM1V,bpmS_CORR_Mon2_RefBPM1V_err] = nancorrcoef(corrMon2_BPMS,corrRefBPM1V_BPMS);
        [bpmV_CORR_Mon2_RefBPM1V,bpmV_CORR_Mon2_RefBPM1V_err] = nancorrcoef(corrMon2_BPMV,corrRefBPM1V_BPMV);
        [bpmV_CORR_Mon2_RefBPM1H,bpmV_CORR_Mon2_RefBPM1H_err] = nancorrcoef(corrMon2_BPMV,corrRefBPM1H_BPMV);
        [bpmV_CORR_Mon2_RefBPM1S,bpmV_CORR_Mon2_RefBPM1S_err] = nancorrcoef(corrMon2_BPMV,corrRefBPM1S_BPMV);
        
        [bpmH_CORR_Mon1_RefBPM1V,bpmH_CORR_Mon1_RefBPM1V_err] = nancorrcoef(corrMon1_BPMH,corrRefBPM1V_BPMH);
        [bpmS_CORR_Mon1_RefBPM1V,bpmS_CORR_Mon1_RefBPM1V_err] = nancorrcoef(corrMon1_BPMS,corrRefBPM1V_BPMS);
        [bpmV_CORR_Mon1_RefBPM1V,bpmV_CORR_Mon1_RefBPM1V_err] = nancorrcoef(corrMon1_BPMV,corrRefBPM1V_BPMV);
        [bpmV_CORR_Mon1_RefBPM1H,bpmV_CORR_Mon1_RefBPM1H_err] = nancorrcoef(corrMon1_BPMV,corrRefBPM1H_BPMV);
        [bpmV_CORR_Mon1_RefBPM1S,bpmV_CORR_Mon1_RefBPM1S_err] = nancorrcoef(corrMon1_BPMV,corrRefBPM1S_BPMV);
        
    end
    
end

if (useRefBPM2)
    
    [bpmH_CORR_Mon3_RefBPM2H,bpmH_CORR_Mon3_RefBPM2H_err] = nancorrcoef(corrMon3_BPMH,corrRefBPM2H_BPMH);
    [bpmS_CORR_Mon3_RefBPM2H,bpmS_CORR_Mon3_RefBPM2H_err] = nancorrcoef(corrMon3_BPMS,corrRefBPM2H_BPMS);
    [bpmH_CORR_Mon3_RefBPM2S,bpmH_CORR_Mon3_RefBPM2S_err] = nancorrcoef(corrMon3_BPMH,corrRefBPM2S_BPMH);
    [bpmS_CORR_Mon3_RefBPM2S,bpmS_CORR_Mon3_RefBPM2S_err] = nancorrcoef(corrMon3_BPMS,corrRefBPM2S_BPMS);
    
    [bpmH_CORR_Mon2_RefBPM2H,bpmH_CORR_Mon2_RefBPM2H_err] = nancorrcoef(corrMon2_BPMH,corrRefBPM2H_BPMH);
    [bpmS_CORR_Mon2_RefBPM2H,bpmS_CORR_Mon2_RefBPM2H_err] = nancorrcoef(corrMon2_BPMS,corrRefBPM2H_BPMS);
    [bpmH_CORR_Mon2_RefBPM2S,bpmH_CORR_Mon2_RefBPM2S_err] = nancorrcoef(corrMon2_BPMH,corrRefBPM2S_BPMH);
    [bpmS_CORR_Mon2_RefBPM2S,bpmS_CORR_Mon2_RefBPM2S_err] = nancorrcoef(corrMon2_BPMS,corrRefBPM2S_BPMS);

    [bpmH_CORR_Mon1_RefBPM2H,bpmH_CORR_Mon1_RefBPM2H_err] = nancorrcoef(corrMon1_BPMH,corrRefBPM2H_BPMH);
    [bpmS_CORR_Mon1_RefBPM2H,bpmS_CORR_Mon1_RefBPM2H_err] = nancorrcoef(corrMon1_BPMS,corrRefBPM2H_BPMS);
    [bpmH_CORR_Mon1_RefBPM2S,bpmH_CORR_Mon1_RefBPM2S_err] = nancorrcoef(corrMon1_BPMH,corrRefBPM2S_BPMH);
    [bpmS_CORR_Mon1_RefBPM2S,bpmS_CORR_Mon1_RefBPM2S_err] = nancorrcoef(corrMon1_BPMS,corrRefBPM2S_BPMS);

    if (useVertical)
        [bpmV_CORR_Mon3_RefBPM2H,bpmV_CORR_Mon3_RefBPM2H_err] = nancorrcoef(corrMon3_BPMV,corrRefBPM2H_BPMV);
        [bpmH_CORR_Mon3_RefBPM2V,bpmH_CORR_Mon3_RefBPM2V_err] = nancorrcoef(corrMon3_BPMH,corrRefBPM2V_BPMH);
        [bpmS_CORR_Mon3_RefBPM2V,bpmS_CORR_Mon3_RefBPM2V_err] = nancorrcoef(corrMon3_BPMS,corrRefBPM2V_BPMS);
        [bpmV_CORR_Mon3_RefBPM2V,bpmV_CORR_Mon3_RefBPM2V_err] = nancorrcoef(corrMon3_BPMV,corrRefBPM2V_BPMV);
        [bpmV_CORR_Mon3_RefBPM2S,bpmV_CORR_Mon3_RefBPM2S_err] = nancorrcoef(corrMon3_BPMV,corrRefBPM2S_BPMV);
        
        [bpmV_CORR_Mon2_RefBPM2H,bpmV_CORR_Mon2_RefBPM2H_err] = nancorrcoef(corrMon2_BPMV,corrRefBPM2H_BPMV);
        [bpmH_CORR_Mon2_RefBPM2V,bpmH_CORR_Mon2_RefBPM2V_err] = nancorrcoef(corrMon2_BPMH,corrRefBPM2V_BPMH);
        [bpmS_CORR_Mon2_RefBPM2V,bpmS_CORR_Mon2_RefBPM2V_err] = nancorrcoef(corrMon2_BPMS,corrRefBPM2V_BPMS);
        [bpmV_CORR_Mon2_RefBPM2V,bpmV_CORR_Mon2_RefBPM2V_err] = nancorrcoef(corrMon2_BPMV,corrRefBPM2V_BPMV);
        [bpmV_CORR_Mon2_RefBPM2S,bpmV_CORR_Mon2_RefBPM2S_err] = nancorrcoef(corrMon2_BPMV,corrRefBPM2S_BPMV);
        
        [bpmV_CORR_Mon1_RefBPM2H,bpmV_CORR_Mon1_RefBPM2H_err] = nancorrcoef(corrMon1_BPMV,corrRefBPM2H_BPMV);
        [bpmH_CORR_Mon1_RefBPM2V,bpmH_CORR_Mon1_RefBPM2V_err] = nancorrcoef(corrMon1_BPMH,corrRefBPM2V_BPMH);
        [bpmS_CORR_Mon1_RefBPM2V,bpmS_CORR_Mon1_RefBPM2V_err] = nancorrcoef(corrMon1_BPMS,corrRefBPM2V_BPMS);
        [bpmV_CORR_Mon1_RefBPM2V,bpmV_CORR_Mon1_RefBPM2V_err] = nancorrcoef(corrMon1_BPMV,corrRefBPM2V_BPMV);
        [bpmV_CORR_Mon1_RefBPM2S,bpmV_CORR_Mon1_RefBPM2S_err] = nancorrcoef(corrMon1_BPMV,corrRefBPM2S_BPMV);

    end
    
end