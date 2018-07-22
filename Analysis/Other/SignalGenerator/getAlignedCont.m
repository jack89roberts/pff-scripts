function [ alignedMixers, alignedDiodes ] = getAlignedCont( mixers, diodes)
%getAlignedCont Version of getAlignedXCorr to use on continuous signals
%(not pulses) from signal generator phase monitor measurements.
%
% Alignment performed using mixer signals then same delays applied to diode
% (near saturation oscillations are not/barely visible on diode making
% alignment ~impossible otherwise, though in principle if oscillation is
% not visible alignment isn't having an effect).

    %% Define arrays, noise calculation
    [nPulses,nSamples] = size(mixers);
    
    containsNaNs = sum(isnan(mixers(:)))>0;
    
    alignedMixers = NaN(nPulses,nSamples); % Jack changed to NaNs from zeros so out of range output is NaN and not 0.
    alignedDiodes = NaN(nPulses,nSamples);
        
    %% Cross correlation
    % Align all pulses to a reference pulse (first pulse above noiseThreshold and, if
    % possible, without NaNs).
    
    refPulse = mixers(1,:);

    % calculate cross correlations and align
    if (nPulses>1)
        for p=1:nPulses

            if (containsNaNs)
                [c,lags] = nanxcorr(refPulse,mixers(p,:),200);
            else
                [c,lags] = xcorr(refPulse,mixers(p,:)); % c = correlation values for sample delay in lags 
            end
            
            [~,maxInd] = max(c); % index of max correlation
            bestLag = lags(maxInd); % delay of max correlation

            % align pulse (both mixer and diode) to the reference (trim points from start if
            % it is too late, add NaNs at start if it is too early).
            if (bestLag<0)
                alignedMixers(p,1:(end-abs(bestLag))) = mixers(p,(abs(bestLag)+1):end);
                alignedDiodes(p,1:(end-abs(bestLag))) = diodes(p,(abs(bestLag)+1):end);
            elseif (bestLag>0)
                alignedMixers(p,(1+abs(bestLag)):end) = mixers(p,1:(end-abs(bestLag)));
                alignedDiodes(p,(1+abs(bestLag)):end) = diodes(p,1:(end-abs(bestLag)));
            elseif (bestLag==0)
                alignedMixers(p,:) = mixers(p,:);
                alignedDiodes(p,:) = diodes(p,:);
            else
                warning('At least one alignment failed (delay from cross correlation undefined)');
            end           
            
        end
    else
        alignedMixers = mixers;
        alignedDiodes = diodes;
    end
   
end

