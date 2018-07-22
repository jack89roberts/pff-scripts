function [ strippedSignal, duplicatePulses ] = removeDuplicatePulses( signal )
%removeDuplicatePulses Replaces any neighbouring duplicate pulses in signal
%with NaNs (leaves the first pulse and removes the subsequent matching
%ones)
    
    strippedSignal = signal;
    [nPulses,nSamples] = size(signal);
    
    i=1;  
    duplicatePulses = zeros(1,nPulses);
    while i < nPulses
        j = i+1;
        
        isDuplicate = (sum(signal(i,:) == signal(j,:))) == ( nSamples - ( sum(isnan(signal(i,:) == signal(j,:))) ) );

        while (isDuplicate && j<nPulses)
            strippedSignal(j,:) = NaN;
            duplicatePulses(j) = 1;
            j = j+1;
            isDuplicate = (sum(signal(i,:) == signal(j,:)) == ( nSamples - ( sum(isnan(signal(i,:) == signal(j,:))) ) ) );
        end

        i = j;
    end

end

