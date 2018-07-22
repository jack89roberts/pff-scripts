function [ mixers, diodes ] = extractMixerDiode( CTFData, useOffset,useSensitivity )
%extractMixerDiode Extracts phase monitor mixer and diode signals (in
%Volts, offset subtracted) from saved data.
%   CTFData: Data struct from mergeMatMonData function.
%
%   useOffset: If 2 subtracts mean of first 10% of points from mixer/diode.
%   If 1 adds the acquired OASIS offset to the mixer/diode. If 0 does no
%   subtraction. Default: 2.
%
%   useSensitivity: If 1, multiplies mixer/diode values by sensitivity.
%   Default: 1
%
%   Always extracts data for all three monitors.
    if (nargin < 3)
        useSensitivity = 1;
    end
    if (nargin < 2)
        useOffset = 2; % 2=mean subtraction, 1=using acquired offset, 0 = no subtraction
    end
    
    
    % Define channel names for all the phase monitors (should be saved
    % somewhere else really).
    mixerNames = {'CT.SCOPE01.CH02', 'CT.SCOPE01.CH04', 'CT.SCOPE01.CH06'};
    diodeNames = {'CT.SCOPE01.CH01', 'CT.SCOPE01.CH03', 'CT.SCOPE01.CH05'};
    acquisitionName = 'Acquisition.value.value';
    sensitivityName = 'Acquisition.sensitivity.value';
    offsetName = 'Acquisition.offset.value';
    nMonitors = length(mixerNames);
    
    % get number of pulses and number of samples
    nPulses = length(CTFData);
    
    
    try 
        nSamples = eval(sprintf('length(CTFData(1).%s.%s)',mixerNames{1},acquisitionName));
    catch
        mixerNames = {'CT_SCOPE01_CH02', 'CT_SCOPE01_CH04', 'CT_SCOPE01_CH06'};
        diodeNames = {'CT_SCOPE01_CH01', 'CT_SCOPE01_CH03', 'CT_SCOPE01_CH05'};
        nSamples = eval(sprintf('length(CTFData(end).%s.%s)',mixerNames{1},acquisitionName)); % Jack - changed to using last pulse instead of first pulse as less likely to be hit by acquisition issues
    end
    
    offsetSamples = 1:round(nSamples./10); % subtract offset using average of 1st 10% of samples (pulse must not arrive in 1st 10%!!!)
    
    % Initialise arrays
    mixerSensitivities = NaN(nMonitors,nPulses);
    diodeSensitivities = NaN(nMonitors,nPulses);
    mixerOffsets = NaN(nMonitors,nPulses);
    diodeOffsets = NaN(nMonitors,nPulses);
    mixers = NaN(nMonitors,nPulses,nSamples);
    diodes = NaN(nMonitors,nPulses,nSamples);
    
    % Extract all the requested monitor signals, multiply by sensitivities
    % and subtract offsets.
    for m=1:nMonitors
        % Extract signals
        tmpMixSensName = [mixerNames{m} '.' sensitivityName];
        mixerSensitivities(m,:) = extractCTFSignalFromMergedData(tmpMixSensName,CTFData);
        
        tmpDioSensName = [diodeNames{m} '.' sensitivityName];
        diodeSensitivities(m,:) = extractCTFSignalFromMergedData(tmpDioSensName,CTFData);
        
        tmpMixName = [mixerNames{m} '.' acquisitionName];
        mixers(m,:,:) = extractCTFSignalFromMergedData(tmpMixName,CTFData);
        mixers = double(mixers);
        
        tmpDioName = [diodeNames{m} '.' acquisitionName];
        diodes(m,:,:) = extractCTFSignalFromMergedData(tmpDioName,CTFData);
        diodes = double(diodes);
        
        % Sensitivities/offsets
        if (useSensitivity == 1)
            for p=1:nPulses
                mixers(m,p,:) = mixers(m,p,:).*mixerSensitivities(m,p);
                diodes(m,p,:) = diodes(m,p,:).*diodeSensitivities(m,p);
                
                mixers(m,p,:) = mixers(m,p,:);
                diodes(m,p,:) = diodes(m,p,:);
            end            
        end
        if (useOffset==2) % mean subtraction
            for p=1:nPulses                
                mixers(m,p,:) = mixers(m,p,:) - nanmean(mixers(m,p,offsetSamples));
                diodes(m,p,:) = diodes(m,p,:) - nanmean(diodes(m,p,offsetSamples));
            end
        elseif (useOffset==1) % using acquired offset
            tmpMixOffsetName = [mixerNames{m} '.' offsetName];
            mixerOffsets(m,:) = extractCTFSignalFromMergedData(tmpMixOffsetName,CTFData);
            tmpDioOffsetName = [diodeNames{m} '.' offsetName];
            diodeOffsets(m,:) = extractCTFSignalFromMergedData(tmpDioOffsetName,CTFData);

            for p=1:nPulses         
                mixers(m,p,:) = mixers(m,p,:) + mixerOffsets(m,p);
                diodes(m,p,:) = diodes(m,p,:) + diodeOffsets(m,p);
            end
        end
    end
    
    if (nPulses == 1)
        mixers = squeeze(mixers);
        diodes = squeeze(diodes);
    end

end

