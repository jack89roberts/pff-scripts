function scanResults = phase_processScanResults(corr465Offsets,calFile,scanDiodes,scanMixers )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    % define some variables, load calibration constants.
    [nCorrOffsets,nAvg,nMons,nSamples] = size(scanDiodes);
    [calibrationConstants, useMixerOverSqrtDiode] = loadFrascatiCalibrationConstants(calFile);
    
    % calculate phases
    scanPhases = NaN(nCorrOffsets,nAvg,nMons,nSamples);
    for offset=1:nCorrOffsets
        for mon=1:nMons
            tmpMixer = squeeze(scanMixers(offset,:,mon,:));
            tmpDiode = squeeze(scanDiodes(offset,:,mon,:));
            if (useMixerOverSqrtDiode)
                scanPhases(offset,:,mon,:) = getPhaseMixerDiode(tmpMixer,tmpDiode,calibrationConstants(mon,1),calibrationConstants(mon,4));
            else
                scanPhases(offset,:,mon,:) = getPhaseMixerDiode(tmpMixer,[],calibrationConstants(mon,1),calibrationConstants(mon,4));
            end
        end
    end
    
    % align signals
    for offset=1:nCorrOffsets
        for mon=1:nMons
            tmpMixer = squeeze(scanMixers(offset,:,mon,:));
            tmpDiode = squeeze(scanDiodes(offset,:,mon,:));
            tmpPhase = squeeze(scanPhases(offset,:,mon,:));
            
            [tmpDiode,tmpOther] = getAlignedXCorr(tmpDiode,'end',{tmpMixer,tmpPhase});
            tmpMixer = tmpOther{1};
            tmpPhase = tmpOther{2};
            
            scanDiodes(offset,:,mon,:) = tmpDiode;
            scanMixers(offset,:,mon,:) = tmpMixer;
            scanPhases(offset,:,mon,:) = tmpPhase;
        end
    end

    % Get sample range
    sampleRange = getFrascatiSampleRange(squeeze(scanMixers(1,1,:,:)),squeeze(scanDiodes(1,1,:,:)));
    
    % Calculate mean phases
    meanPhases = nanmean(nanmean(scanPhases(:,:,:,sampleRange),4),2);
    stdPhases = nanstd(nanmean(scanPhases(:,:,:,sampleRange),4),0,2);
    stdPhases = squeeze(stdPhases);
    
    % fit data
    [fitCoeffs,fitRsquare,fitConfInts] = nanpolyfit(corr465Offsets, meanPhases(:,3), 1);
    phaseShiftPerAmp = fitCoeffs(1);
    fprintf('Phase shift per amp: %.2f degrees\n',phaseShiftPerAmp);
    fprintf('Fit R square: %.4f\n',fitRsquare);
    
    % plot data
    figure;
    fitVals = phaseShiftPerAmp.*linspace(min(corr465Offsets),max(corr465Offsets),100) + fitCoeffs(2);
    plot(linspace(min(corr465Offsets),max(corr465Offsets),100),fitVals);
    hold on
    errorbar(corr465Offsets,meanPhases(:,3),stdPhases(:,3)./sqrt(nAvg),'o');
    xlabel('Corr 465 Offset [A]');
    ylabel('Phase [degrees]');
    grid on;
    
    % save results to data struct
    scanResults = struct();
    scanResults.calibrationConstants=calibrationConstants;
    scanResults.useMixerOverSqrtDiode=useMixerOverSqrtDiode;
    scanResults.scanPhases=scanPhases;
    scanResults.scanDiodes=scanDiodes;
    scanResults.scanMixers=scanMixers;
    scanResults.sampleRange=sampleRange;
    scanResults.meanPhases=meanPhases;
    scanResults.stdPhases=stdPhases;
    scanResults.fitCoeffs=fitCoeffs;
    scanResults.fitRsquare=fitRsquare;
    scanResults.fitConfInts=fitConfInts;
    scanResults.phaseShiftPerAmp=phaseShiftPerAmp;
    scanResults.fitVals=fitVals;
    
end

