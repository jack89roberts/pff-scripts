function driftRemSignal = removeDrift( signal, sampleRange, nAvg )
%removeDrift Subtracts drift from signal
%   sampleRange: Drift is subtracted based on the mean across these
%   samples. If not given, uses mean across all samples.
%   nAvg: Subtracted drift is mean of this many points. If an even number
%   is given one will be added to it to create a symmetric window about
%   each sample point. Default: nPulses/10.

    [nPulses,nSamples] = size(signal);
    
    % deal with the case where signal is one dimensional
    if (nPulses==1 || nSamples==1)
        nPulses = length(signal);
        nSamples = 1;
        sampleRange = 1;
    end
    
    if (nargin<3)
        %nAvg = 5;
        nAvg = round(nPulses/10);
    end
    % use an odd number to average (value at sample is average of that
    % sample plus (nAvg-1)/2 samples before and after it
    if (mod(nAvg,2) == 0)
        nAvg = nAvg + 1;
    end
    % if nAvg became bigger than nPulses, make it the next lowest odd
    % number.
    if (nAvg > nPulses)
        if (mod(nPulses,2)==0)
            nAvg = nPulses-1;
        else
            nAvg = nPulses;
        end
    end
    nNeighbours = (nAvg-1)/2;
    
    % if no sample range given, use all the samples
    if (nargin<2 || isempty(sampleRange))
        sampleRange = 1:nSamples;
    end
    
    % calculate pulse means in signal
    if (nSamples>1)
        meanSignal = nanmean(signal(:,sampleRange),2);
    else
        meanSignal = signal;
    end
    
    driftRemSignal = signal;
    for p=1:nPulses
        
        % calculate mean to remove
        if (p < nNeighbours+1)
        % use first nAvg pulses whilst there are not nNeighbours points
        % available before p
            avgSignal = nanmean(meanSignal(1:nAvg));
            
        elseif (p > nPulses-nNeighbours)
        % use last nAvg pulses whilst there are not nNeighbours points
        % available after p
            avgSignal = nanmean(meanSignal((nPulses-nAvg+1):nPulses));
        
        else
        % use p-nNeighbours:p+nNeighbours
            avgSignal = nanmean(meanSignal((p-nNeighbours):(p+nNeighbours)));
        
        end
        
        % subtract avgSignal from signal
        if (nSamples>1)
            driftRemSignal(p,:) = driftRemSignal(p,:)-avgSignal;
        else
            driftRemSignal(p) = driftRemSignal(p)-avgSignal;
        end
    end
    
end

