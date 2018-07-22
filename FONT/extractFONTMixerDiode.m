function [ mixers,diodes ] = extractFONTMixerDiode( ADCs, useOffset, mixerADCs, diodeIsPresent, diodeADCs )
%extractFONTMixerDiode Take ADCs from extracted FONT data and return arrays
%of mixers and diodes in format used for other analysis scripts.
    if nargin<2
        useOffset = 1;
    end 
    % ADC number of each mixer and diode channel.
    if (nargin < 3 || isempty(mixerADCs))
        mixerADCs = [2,4,6];
    end
    if (nargin<4)
        diodeIsPresent = 1;
    end
    if (nargin<5 || isempty(diodeADCs))
        diodeADCs = [1,3,5];  
    end
    
    % convert in to (mon,pulse,sample) dimension order as used in other
    % analysis scripts
    %ADCs = permute(ADCs,[1 3 2]); % Jack - now done by default in my
    %version of data extracting function
    
    mixers = ADCs(mixerADCs,:,:);
    if (diodeIsPresent)
        diodes = ADCs(diodeADCs,:,:);
    else
        diodes = [];
    end
    
    [nMonitors,nPulses,nSamples] = size(mixers);
    offsetSamples = 1:round(nSamples./10); % subtract offset using average of 1st 10% of samples (pulse must not arrive in 1st 10%!!!)
    
    % baseline subtraction
    for m=1:nMonitors
        if (useOffset==1) % mean subtraction
            for p=1:nPulses                
                mixers(m,p,:) = mixers(m,p,:) - nanmean(mixers(m,p,offsetSamples));
                if (diodeIsPresent)
                    diodes(m,p,:) = diodes(m,p,:) - nanmean(diodes(m,p,offsetSamples));
                end
            end
        end
    end
    
end

