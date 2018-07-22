function [ phaseData, dataIsGood ] = extractPhaseData(dataStruct, frascatiCalStruct)
    % Have removed:
    % Initial phase subtraction
    % Time axis calculation (based on pulse start sample)
    frascatiCalibrationFactors = frascatiCalStruct.calibrationFactors;
    useMixerOverSqrtDiode = frascatiCalStruct.useMixerOverSqrtDiode;
    
    try
        [frascatiMixers,frascatiDiodes] = extractMixerDiode(dataStruct); % JACK - Changed extractMixerDiode function to try to deal with '.' being replaced with '_'
        frascatiMixers = squeeze(frascatiMixers);
        frascatiDiodes = squeeze(frascatiDiodes);

        [~,frascatiNSamples] = size(frascatiMixers);
        [frascatiNMonitors,~] = size(frascatiCalibrationFactors);
        frascatiPhases = NaN*ones(frascatiNMonitors,frascatiNSamples);

        for mon = 1:frascatiNMonitors
            if (useMixerOverSqrtDiode)
                frascatiPhases(mon,:,:) = getPhaseMixerDiode(frascatiMixers(mon,:,:),...
                                                     frascatiDiodes(mon,:,:),...
                                                     frascatiCalibrationFactors(mon,1),...
                                                     frascatiCalibrationFactors(mon,4)...
                                                    );
            else
                frascatiPhases(mon,:,:) = getPhaseMixerDiode(frascatiMixers(mon,:,:),...
                                                     [],...
                                                     frascatiCalibrationFactors(mon,1),...
                                                     frascatiCalibrationFactors(mon,4)...
                                                    );
            end
        end
        frascatiPhases = squeeze(frascatiPhases);

        try
            petsI = extractCTFSignalFromMergedData('CE_SCOPE03_CH01.Acquisition.value.value',dataStruct);  % JACK - Changed to try to deal with '.' being replaced with '_'
            petsQ = extractCTFSignalFromMergedData('CE_SCOPE03_CH02.Acquisition.value.value',dataStruct);
            petsISensitivity = extractCTFSignalFromMergedData('CE_SCOPE03_CH01.Acquisition.sensitivity.value',dataStruct);
            petsQSensitivity = extractCTFSignalFromMergedData('CE_SCOPE03_CH02.Acquisition.sensitivity.value',dataStruct);
        catch
            petsI = extractCTFSignalFromMergedData('CE.SCOPE03.CH01.Acquisition.value.value',dataStruct);  % JACK - Changed to try to deal with '.' being replaced with '_'
            petsQ = extractCTFSignalFromMergedData('CE.SCOPE03.CH02.Acquisition.value.value',dataStruct);
            petsISensitivity = extractCTFSignalFromMergedData('CE.SCOPE03.CH01.Acquisition.sensitivity.value',dataStruct);
            petsQSensitivity = extractCTFSignalFromMergedData('CE.SCOPE03.CH02.Acquisition.sensitivity.value',dataStruct);
        end

        petsNSamples = length(petsI);
        petsI = double(petsI).*petsISensitivity;
        petsQ = double(petsQ).*petsQSensitivity;
        [petsPhase,petsPower] = getPhaseIQ(petsI, petsQ);
        petsPhase = squeeze(petsPhase);

        try
            frascatiTimePerSample = extractCTFSignalFromMergedData('CT_SCOPE01_CH02.Acquisition.sampleInterval.value',dataStruct); % JACK - Changed to try to deal with '.' being replaced with '_'
            petsTimePerSample = extractCTFSignalFromMergedData('CE_SCOPE03_CH01.Acquisition.sampleInterval.value',dataStruct);
        catch  
            frascatiTimePerSample = extractCTFSignalFromMergedData('CT.SCOPE01.CH02.Acquisition.sampleInterval.value',dataStruct); % JACK - Changed to try to deal with '.' being replaced with '_'
            petsTimePerSample = extractCTFSignalFromMergedData('CE.SCOPE03.CH01.Acquisition.sampleInterval.value',dataStruct);
        end
        
        phaseData = struct();
        phaseData.frascatiMixers = frascatiMixers;
        phaseData.frascatiDiodes = frascatiDiodes;
        phaseData.frascatiPhases = frascatiPhases;
        phaseData.petsI = petsI;
        phaseData.petsQ = petsQ;
        phaseData.petsPower = petsPower;
        phaseData.petsPhase = petsPhase;
        phaseData.frascatiNSamples = frascatiNSamples;
        phaseData.frascatiNMonitors = frascatiNMonitors;
        phaseData.petsNSamples = petsNSamples;
        phaseData.frascatiTimePerSample = frascatiTimePerSample;
        phaseData.petsTimePerSample = petsTimePerSample;

        if (isempty(frascatiMixers) || isnan(min(frascatiMixers(:))) || isempty(petsI) || isnan(min(petsI(:))))
            dataIsGood = 0;
        else
            dataIsGood = 1;
        end
        
    catch
        phaseData = struct();
        dataIsGood = 0;
    end
end