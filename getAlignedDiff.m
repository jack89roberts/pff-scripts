function [ aligned, foundPeak ] = getAlignedDiff( trace )
%getAlignedDiff Tries to align pulses in trace using the minimum and maximum value of its diff. 
%   Detailed explanation goes here


    nBaseline = 20;
    nStd = 5;
    
    [nPulses,nSamples] = size(trace);
    nPointsAfterBaseline = nSamples-nBaseline;
    
    diffTrace = diff(trace,1,2);
    
    %lintrace = reshape(diffTrace(:,1:nBaseline),1,nPulses*nBaseline);
    %noiseThreshold = nanstd(lintrace)*nStd;
    noiseThreshold = nStd.*nanstd(diffTrace(:));
    
    [maxVals,maxInds] = max(diffTrace,[],2);
    [minVals,minInds] = min(diffTrace,[],2);
    foundPeak = (abs(maxVals)>noiseThreshold) | (abs(minVals)>noiseThreshold);
    
    startSamples = min([maxInds minInds],[],2);
    endSamples = min([maxInds minInds],[],2);

    aligned = NaN(nPulses,nSamples);
    
    for p=1:nPulses
        
        if (foundPeak)
            if (startSamples(p)>nBaseline)
                aligned(p,1:nBaseline) = trace(p, (startSamples(p)-nBaseline):(startSamples(p)-1));
            else
                aligned(p,(nBaseline-startSamples(p)+2):nBaseline) = trace(p, 1:(startSamples(p)-1));
            end

            nPointsAfterStart = nSamples-startSamples(p)+1;

            if (nPointsAfterStart < nPointsAfterBaseline)
                aligned(p, (nBaseline+1):(nBaseline+nPointsAfterStart)) = trace(p, startSamples(p):end);
            else
                aligned(p,(nBaseline+1):end) = trace(p, startSamples(p):(startSamples(p)+nPointsAfterBaseline-1));
            end
        end
        
    end
    
end

