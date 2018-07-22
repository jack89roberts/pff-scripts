function [ phaseStatistics ] = calculatePhaseStatistics( processedPhaseData)
%calculatePhaseStatistics Takes processedPhaseData struct and uses it to
%calculate some phase statistics.
%   Detailed explanation goes here

    frascatiNMonitors = processedPhaseData.frascatiNMonitors;
    [~,nPulses,~] = size(processedPhaseData.frascatiMixers);
    frascatiStartSamples = processedPhaseData.frascatiStartSamples;
    frascatiEndSamples = processedPhaseData.frascatiEndSamples;
    petsStartSample = processedPhaseData.petsStartSample;
    petsEndSample = processedPhaseData.petsEndSample;
    
    % Means and standard deviations for the Frascati monitors
    frascatiMeanSamplePhases = squeeze(nanmean(processedPhaseData.frascatiPhases,2));
    frascatiStdSamplePhases = squeeze(nanstd(processedPhaseData.frascatiPhases,0,2));
    frascatiMeanSampleDiodes = squeeze(nanmean(processedPhaseData.frascatiDiodes,2));
    frascatiStdSampleDiodes = squeeze(nanstd(processedPhaseData.frascatiDiodes,0,2));
    frascatiMeanSampleMixers = squeeze(nanmean(processedPhaseData.frascatiMixers,2));
    frascatiStdSampleMixers = squeeze(nanstd(processedPhaseData.frascatiMixers,0,2));

    frascatiMeanStdSamplePhases = NaN*ones(1,processedPhaseData.frascatiNMonitors);
    frascatiMeanStdSampleDiodes = NaN*ones(1,processedPhaseData.frascatiNMonitors);
    frascatiMeanStdSampleMixers = NaN*ones(1,processedPhaseData.frascatiNMonitors);

    frascatiMeanPulsePhases = NaN*ones(frascatiNMonitors,nPulses);
    frascatiMeanPulseDiodes = NaN*ones(frascatiNMonitors,nPulses);
    frascatiMeanPulseMixers = NaN*ones(frascatiNMonitors,nPulses);

    for mon=1:frascatiNMonitors
        frascatiMeanStdSamplePhases(mon) = squeeze(nanmean(frascatiStdSamplePhases(mon,frascatiStartSamples(mon):frascatiEndSamples(mon)),2));
        frascatiMeanStdSampleDiodes(mon) = squeeze(nanmean(frascatiStdSampleDiodes(mon,frascatiStartSamples(mon):frascatiEndSamples(mon)),2));
        frascatiMeanStdSampleMixers(mon) = squeeze(nanmean(frascatiStdSampleMixers(mon,frascatiStartSamples(mon):frascatiEndSamples(mon)),2));

        frascatiMeanPulsePhases(mon,:) = squeeze(nanmean(processedPhaseData.frascatiPhases(mon,:,frascatiStartSamples(mon):frascatiEndSamples(mon)),3));
        frascatiMeanPulseDiodes(mon,:) = squeeze(nanmean(processedPhaseData.frascatiDiodes(mon,:,frascatiStartSamples(mon):frascatiEndSamples(mon)),3));
        frascatiMeanPulseMixers(mon,:) = squeeze(nanmean(processedPhaseData.frascatiMixers(mon,:,frascatiStartSamples(mon):frascatiEndSamples(mon)),3));


    end

    frascatiStdPulsePhases = nanstd(frascatiMeanPulsePhases,0,2); 
    frascatiStdPulseDiodes = nanstd(frascatiMeanPulseDiodes,0,2); 
    frascatiStdPulseMixers = nanstd(frascatiMeanPulseMixers,0,2); 


    % Means and standard deviations for the PETS
    petsMeanSamplePhase = nanmean(processedPhaseData.petsPhase);
    petsMeanPulsePhase = nanmean(processedPhaseData.petsPhase(:,petsStartSample:petsEndSample),2);
    petsStdSamplePhase = nanstd(processedPhaseData.petsPhase);
    petsMeanStdSamplePhase = nanmean(petsStdSamplePhase(petsStartSample:petsEndSample));
    petsStdPulsePhase = nanstd(petsMeanPulsePhase);


    % Correlations
    corr12_meanPulsePhase = corrcoef(frascatiMeanPulsePhases(1,:),frascatiMeanPulsePhases(2,:));           
    corr13_meanPulsePhase = corrcoef(frascatiMeanPulsePhases(1,:),frascatiMeanPulsePhases(3,:));
    corr23_meanPulsePhase = corrcoef(frascatiMeanPulsePhases(2,:),frascatiMeanPulsePhases(3,:));
    corrPetsFrascati3_meanPulsePhase = corrcoef(frascatiMeanPulsePhases(3,:),petsMeanPulsePhase);

    corr12_meanPulsePhase = corr12_meanPulsePhase(1,2);
    corr13_meanPulsePhase = corr13_meanPulsePhase(1,2);
    corr23_meanPulsePhase = corr23_meanPulsePhase(1,2);
    corrPetsFrascati3_meanPulsePhase = corrPetsFrascati3_meanPulsePhase(1,2);

    % Fits
    fit12_meanPulsePhase = polyfit(frascatiMeanPulsePhases(1,:),frascatiMeanPulsePhases(2,:),1);
    fit13_meanPulsePhase = polyfit(frascatiMeanPulsePhases(1,:),frascatiMeanPulsePhases(3,:),1);
    fit23_meanPulsePhase = polyfit(frascatiMeanPulsePhases(2,:),frascatiMeanPulsePhases(3,:),1);
    fitPetsFrascati3_meanPulsePhase = polyfit(frascatiMeanPulsePhases(3,:),petsMeanPulsePhase',1);

    % save everything to the phaseStatistics struct
    phaseStatistics = struct();
    phaseStatistics.frascatiMeanSamplePhases = frascatiMeanSamplePhases;
    phaseStatistics.frascatiStdSamplePhases = frascatiStdSamplePhases;
    phaseStatistics.frascatiMeanStdSamplePhases = frascatiMeanStdSamplePhases;
    phaseStatistics.frascatiMeanPulsePhases = frascatiMeanPulsePhases;
    phaseStatistics.frascatiStdPulsePhases = frascatiStdPulsePhases;

    phaseStatistics.frascatiMeanSampleDiodes = frascatiMeanSampleDiodes;
    phaseStatistics.frascatiStdSampleDiodes = frascatiStdSampleDiodes;
    phaseStatistics.frascatiMeanStdSampleDiodes = frascatiMeanStdSampleDiodes;
    phaseStatistics.frascatiMeanPulseDiodes = frascatiMeanPulseDiodes;
    phaseStatistics.frascatiStdPulseDiodes = frascatiStdPulseDiodes;

    phaseStatistics.frascatiMeanSampleMixers = frascatiMeanSampleMixers;
    phaseStatistics.frascatiStdSampleMixers = frascatiStdSampleMixers;
    phaseStatistics.frascatiMeanStdSampleMixers = frascatiMeanStdSampleMixers;
    phaseStatistics.frascatiMeanPulseMixers = frascatiMeanPulseMixers;
    phaseStatistics.frascatiStdPulseMixers = frascatiStdPulseMixers;

    phaseStatistics.petsMeanSamplePhase = petsMeanSamplePhase;
    phaseStatistics.petsMeanPulsePhase = petsMeanPulsePhase;
    phaseStatistics.petsMeanStdSamplePhase = petsMeanStdSamplePhase;
    phaseStatistics.petsStdSamplePhase = petsStdSamplePhase;
    phaseStatistics.petsStdPulsePhase = petsStdPulsePhase;

    phaseStatistics.corr12_meanPulsePhase = corr12_meanPulsePhase;
    phaseStatistics.corr13_meanPulsePhase = corr13_meanPulsePhase;
    phaseStatistics.corr23_meanPulsePhase = corr23_meanPulsePhase;
    phaseStatistics.corrPetsFrascati3_meanPulsePhase = corrPetsFrascati3_meanPulsePhase;           

    phaseStatistics.fit12_meanPulsePhase = fit12_meanPulsePhase;
    phaseStatistics.fit13_meanPulsePhase = fit13_meanPulsePhase;
    phaseStatistics.fit23_meanPulsePhase = fit23_meanPulsePhase;
    phaseStatistics.fitPetsFrascati3_meanPulsePhase = fitPetsFrascati3_meanPulsePhase;    

end

