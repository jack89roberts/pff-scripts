function [ treatedPhaseData, dataIsGood ] = treatPhaseData(dataStruct, frascatiCalibrationFactors)
    % Have removed:
    % Initial phase subtraction
    % Time axis calculation (based on pulse start sample)
    try
        [frascatiMixers,frascatiDiodes] = extractMixerDiode(dataStruct); % JACK - Changed extractMixerDiode function to try to deal with '.' being replaced with '_'
        frascatiMixers = squeeze(frascatiMixers);
        frascatiDiodes = squeeze(frascatiDiodes);

        [~,frascatiNSamples] = size(frascatiMixers);
        [frascatiNMonitors,~] = size(frascatiCalibrationFactors);
        frascatiPhases = NaN*ones(frascatiNMonitors,frascatiNSamples);

        for mon = 1:frascatiNMonitors
            frascatiPhases(mon,:,:) = getPhaseMixerDiode(frascatiMixers(mon,:,:),...
                                                 frascatiDiodes(mon,:,:),...
                                                 frascatiCalibrationFactors(mon,1),...
                                                 frascatiCalibrationFactors(mon,4)...
                                                );
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
        
        treatedPhaseData = struct();
        treatedPhaseData.frascatiMixers = frascatiMixers;
        treatedPhaseData.frascatiDiodes = frascatiDiodes;
        treatedPhaseData.frascatiPhases = frascatiPhases;
        treatedPhaseData.petsI = petsI;
        treatedPhaseData.petsQ = petsQ;
        treatedPhaseData.petsPower = petsPower;
        treatedPhaseData.petsPhase = petsPhase;
        treatedPhaseData.frascatiNSamples = frascatiNSamples;
        treatedPhaseData.frascatiNMonitors = frascatiNMonitors;
        treatedPhaseData.petsNSamples = petsNSamples;
        treatedPhaseData.frascatiTimePerSample = frascatiTimePerSample;
        treatedPhaseData.petsTimePerSample = petsTimePerSample;

        if (isempty(frascatiMixers) || isnan(min(frascatiMixers(:))) || isempty(petsI) || isnan(min(petsI(:))))
            dataIsGood = 0;
        else
            dataIsGood = 1;
        end
        
    catch
        treatedPhaseData = struct();
        dataIsGood = 0;
    end
end