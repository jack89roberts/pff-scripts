function [ posH, posV ] = extractMonPos( CTFData, useOffset,useSensitivity )
%extractMonPos Extracts phase monitor position signals (in
%Volts, offset subtracted) from saved data.
%   CTFData: Data struct from mergeMatMonData function.
%
%   useOffset: If 2 subtracts mean of first 10% of points from signals.
%   If 1 adds the acquired OASIS offset to the mixer/diode. If 0 does no
%   subtraction. Default: 2.
%
%   useSensitivity: If 1, multiplies signal values by sensitivity.
%   Default: 1
%
    if (nargin < 3)
        useSensitivity = 1;
    end
    if (nargin < 2)
        useOffset = 2; % 2=mean subtraction, 1=using acquired offset, 0 = no subtraction
    end
    
    
    % Define channel names for all the phase monitors (should be saved
    % somewhere else really).
    posHNames = {'CT.SCOPE01.CH07', 'CT.SCOPE02.CH01'};
    posVNames = {'CT.SCOPE01.CH08', 'CT.SCOPE02.CH02'};
    acquisitionName = 'Acquisition.value.value';
    sensitivityName = 'Acquisition.sensitivity.value';
    offsetName = 'Acquisition.offset.value';
    nMonitors = length(posHNames);
    
    % get number of pulses and number of samples
    nPulses = length(CTFData);
    
    
    try 
        nSamples = eval(sprintf('length(CTFData(1).%s.%s)',posHNames{1},acquisitionName));
    catch
        posHNames = {'CT_SCOPE01_CH07', 'CT_SCOPE02_CH01'};
        posVNames = {'CT_SCOPE01_CH08', 'CT_SCOPE02_CH02'};
        nSamples = eval(sprintf('length(CTFData(end).%s.%s)',posHNames{1},acquisitionName)); % Jack - changed to using last pulse instead of first pulse as less likely to be hit by acquisition issues
    end
    
    offsetSamples = 1:round(nSamples./10); % subtract offset using average of 1st 10% of samples (pulse must not arrive in 1st 10%!!!)
    
    % Initialise arrays
    hSensitivities = NaN(nMonitors,nPulses);
    vSensitivities = NaN(nMonitors,nPulses);
    hOffsets = NaN(nMonitors,nPulses);
    vOffsets = NaN(nMonitors,nPulses);
    posH = NaN(nMonitors,nPulses,nSamples);
    posV = NaN(nMonitors,nPulses,nSamples);
    
    % Extract all the requested monitor signals, multiply by sensitivities
    % and subtract offsets.
    for m=1:nMonitors
        % Extract signals
        tmpHSensName = [posHNames{m} '.' sensitivityName];
        hSensitivities(m,:) = extractCTFSignalFromMergedData(tmpHSensName,CTFData);
        
        tmpVSensName = [posVNames{m} '.' sensitivityName];
        vSensitivities(m,:) = extractCTFSignalFromMergedData(tmpVSensName,CTFData);
        
        tmpHName = [posHNames{m} '.' acquisitionName];
        posH(m,:,:) = extractCTFSignalFromMergedData(tmpHName,CTFData);
        posH = double(posH);
        
        tmpVName = [posVNames{m} '.' acquisitionName];
        posV(m,:,:) = extractCTFSignalFromMergedData(tmpVName,CTFData);
        posV = double(posV);
        
        % Sensitivities/offsets
        if (useSensitivity == 1)
            for p=1:nPulses
                posH(m,p,:) = posH(m,p,:).*hSensitivities(m,p);
                posV(m,p,:) = posV(m,p,:).*vSensitivities(m,p);
                
                posH(m,p,:) = posH(m,p,:);
                posV(m,p,:) = posV(m,p,:);
            end            
        end
        if (useOffset==2) % mean subtraction
            for p=1:nPulses                
                posH(m,p,:) = posH(m,p,:) - nanmean(posH(m,p,offsetSamples));
                posV(m,p,:) = posV(m,p,:) - nanmean(posV(m,p,offsetSamples));
            end
        elseif (useOffset==1) % using acquired offset
            tmpHOffsetName = [posHNames{m} '.' offsetName];
            hOffsets(m,:) = extractCTFSignalFromMergedData(tmpHOffsetName,CTFData);
            tmpVOffsetName = [posVNames{m} '.' offsetName];
            vOffsets(m,:) = extractCTFSignalFromMergedData(tmpVOffsetName,CTFData);

            for p=1:nPulses         
                posH(m,p,:) = posH(m,p,:) + hOffsets(m,p);
                posV(m,p,:) = posV(m,p,:) + vOffsets(m,p);
            end
        end
    end
    
    if (nPulses == 1)
        posH = squeeze(posH);
        posV = squeeze(posV);
    end

end

