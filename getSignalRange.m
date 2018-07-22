function [ startSample, endSample ] = getSignalRange( signal, skipPercent )
%getSignalRange Attempts to find range of a signal
%   skipPercent - remove this percentage of samples from the beginning and
%   end of the found range.
    diffSig = diff(signal);
    [~,slope1Samp] = min(diffSig);
    [~,slope2Samp] = max(diffSig);
    
    minSamp = min([slope1Samp slope2Samp]);
    maxSamp = max([slope1Samp slope2Samp]);
    
    sigWidth = maxSamp - minSamp;
    samplesToSkip = ceil((skipPercent/100).*sigWidth);
    
    startSample = minSamp + samplesToSkip;
    endSample = maxSamp - samplesToSkip;

end

