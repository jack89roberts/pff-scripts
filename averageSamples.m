function [averagedSignal,midSamples] = averageSamples( signal, nAvg, typeAvg )
%averageSamples Performs a running average of nAvg samples for all the
%pulses in signal.
%averagedSignal = averageSamples( signal, nAvg )
% typeAvg: rolling or box

    [nPulses,nSamples] = size(signal);
    
    % deal with the case where signal is one dimensional
    if (nPulses==1 || nSamples==1)
        if (nPulses>nSamples)
            signal = signal';
        end
        nPulses = 1;
        nSamples = length(signal);
    end
    
    if (nargin<2)
        %nAvg = 5;
        nAvg = round(nSamples/100);
    end
    
    if (nargin<3)
        typeAvg = 'rolling';
    elseif (~(strcmp(typeAvg,'rolling') || strcmp(typeAvg,'box')))
        error('typeAvg must be rolling or box')
    end
    
    % use an odd number to average (value at sample is average of that
    % sample plus (nAvg-1)/2 samples before and after it
    if (mod(nAvg,2) == 0)
        nAvg = nAvg + 1;
    end
    % if nAvg became bigger than nSamples, make it the next lowest odd
    % number.
    if (nAvg > nSamples)
        if (mod(nSamples,2)==0)
            nAvg = nSamples-1;
        else
            nAvg = nSamples;
        end
    end
    nNeighbours = (nAvg-1)/2;
    
    averagedSignal = NaN(size(signal));

    if (strcmp(typeAvg,'rolling'))
        midSamples = 1:nSamples;
        
        for p=1:nPulses

            for s=1:nSamples

                % calculate mean to remove
                if (s < nNeighbours+1)
                % use first nAvg samples whilst there are not nNeighbours points
                % available before s
                    avgSignal = nanmean(signal(p,1:nAvg));

                elseif (s > nSamples-nNeighbours)
                % use last nAvg samples whilst there are not nNeighbours points
                % available after s
                    avgSignal = nanmean(signal(p,(nSamples-nAvg+1):nSamples));

                else
                % use s-nNeighbours:s+nNeighbours
                    avgSignal = nanmean(signal(p,(s-nNeighbours):(s+nNeighbours)));

                end

                averagedSignal(p,s) = avgSignal;

            end
        end
        
    else % box average
        midSamples = (nNeighbours+1):nAvg:(nSamples-nNeighbours);
        
        for s=midSamples              
               averagedSignal(:,s) = nanmean(signal(:,(s-nNeighbours):(s+nNeighbours)),2);
        end
        
        averagedSignal = averagedSignal(:,midSamples);
        
    end
    
end

