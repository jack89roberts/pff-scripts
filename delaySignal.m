function delayedSignal = delaySignal( signal, nDelay )
% delaySignal Shifts a signal in time by nDelay samples along second
% dimension.
%   signal(pulses,samples)
%   nDelay>0: shifts signal later, nDelay<0: shifts signal earlier.

    if (iscell(signal)); error('signal must not be a cell'); end;
    if (ndims(signal)>2); error('signal must not have more than two dimensions'); end;
    
    [nRows,nColumns] = size(signal);
    if ( (nRows==1 || nColumns==1) && (nRows>nColumns))
        signal = signal';
    end

    delayedSignal = NaN(size(signal));
    
    if (nDelay>0) % need to move pulse later
        delayedSignal(:,(1+nDelay):end) = signal(:,1:(end-nDelay));
    elseif (nDelay<0) % need to move pulse earlier
        delayedSignal(:,1:(end-abs(nDelay))) = signal(:,(abs(nDelay)+1):end);
    else
        delayedSignal = signal;
    end
    
    if ( (nRows==1 || nColumns==1) && (nRows>nColumns))
        delayedSignal = delayedSignal';
    end


end

