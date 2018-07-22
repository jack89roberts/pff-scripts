function [ ampSigs ] = extractAmpMon( CTFData,useOffset,useSensitivity,convertToAmpVolts )
%extractAmpMon Extracts amplifier monitoring signals from CTFData.
% ampSigs: 1=LA, 2=LB, 3=RA, 4=RB
    if (nargin < 4)
        convertToAmpVolts = 1; % convert ADC voltages in to actual amplifier output voltage.
    end
    if (nargin < 3)
        useSensitivity = 1;
    end
    if (nargin < 2)
        useOffset = 1; % 2=mean subtraction, 1=using acquired offset, 0 = no subtraction
    end
    
    
    % Define channel names for all the phase monitors (should be saved
    % somewhere else really).
    sigNames = {'CT_SCOPE02_CH05', 'CT_SCOPE02_CH06', 'CT_SCOPE02_CH07', 'CT_SCOPE02_CH08'};
    acquisitionName = 'Acquisition.value.value';
    sensitivityName = 'Acquisition.sensitivity.value';
    offsetName = 'Acquisition.offset.value';
    nSignals = length(sigNames);
    
    % get number of pulses and number of samples
    nPulses = length(CTFData);
    
    try
        nSamples = eval(sprintf('length(CTFData(end).%s.%s)',sigNames{1},acquisitionName));
    catch
        error('Amplifier monitoring signals not present');
    end
    offsetSamples = 1:round(nSamples./10); % subtract offset using average of 1st 10% of samples (pulse must not arrive in 1st 10%!!!)
    
    % Initialise arrays
    signalSensitivities = NaN(nSignals,nPulses);
    signalOffsets = NaN(nSignals,nPulses);
    ampSigs = NaN(nSignals,nPulses,nSamples);
    
    % Extract all signals, multiply by sensitivities and subtract offsets.
    for s=1:nSignals
        % Extract signals
        tmpSensName = [sigNames{s} '.' sensitivityName];
        signalSensitivities(s,:) = extractCTFSignalFromMergedData(tmpSensName,CTFData);
                
        tmpSigName = [sigNames{s} '.' acquisitionName];
        ampSigs(s,:,:) = extractCTFSignalFromMergedData(tmpSigName,CTFData);
        ampSigs = double(ampSigs);
        
        % Sensitivities/offsets
        if (useSensitivity == 1)
            for p=1:nPulses
                ampSigs(s,p,:) = ampSigs(s,p,:).*signalSensitivities(s,p);                
            end            
        end
        if (useOffset==2) % mean subtraction
            for p=1:nPulses                
                ampSigs(s,p,:) = ampSigs(s,p,:) - nanmean(ampSigs(s,p,offsetSamples));
            end
        elseif (useOffset==1) % using acquired offset
            tmpMixOffsetName = [sigNames{s} '.' offsetName];
            signalOffsets(s,:) = extractCTFSignalFromMergedData(tmpMixOffsetName,CTFData);
            for p=1:nPulses         
                ampSigs(s,p,:) = ampSigs(s,p,:) + signalOffsets(s,p);
            end
        end
    end
    
    if (convertToAmpVolts)
        ampSigs = ampSigs*(115*4); % amp monitoring outputs factor 115, 12dB before SiS
    end
    
    if (nPulses == 1)
        ampSigs = squeeze(ampSigs);
    end

end

