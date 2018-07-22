function [ aligned, otherAligned, pulseRange ] = getAlignedStd( trace, otherTraces )
%getAlignedStd Tries to align pulses in trace using a running standard
%deviation.
%   rising and falling edge cause two peaks in a running standard deviation
%   across nAvg points. Pulse range is defined as between the two tallest
%   peaks in this standard deviation.
%   VERY SLOW, e.g. at least 10 times slower than getAligned2, but may
%   produce better results in some cases. Particularly for pulse range
%   calculation (alignment itself is usually good in getAligned2).

    nBaseline = 20; % pulse start is moved to sample nBaseline + 1
    nStd = 5; % check that a peak exists with amplitude above this many standard deviations
    nAvg = 5; % odd number. Number of points to include in each standard deviation calculation.
    %%
    nNeighbours = ceil((nAvg-1)./2); % range of each std calculation is (sample-nNeighbours):(sample+nNeighbours)
    
    [nPulses,nSamples] = size(trace);
    nPointsAfterBaseline = nSamples-nBaseline; % used for constructing aligned pulse

    
    if (nargin < 2)
        otherTraces = {};
        otherAligned = {};
    else
        nOtherTraces = length(otherTraces);
        otherAligned = cell(1,nOtherTraces);
        for i=1:nOtherTraces
            [nRows,~] = size(otherTraces{i});
            if (nRows ~= nPulses)
                error('otherTraces must have same no. pulses as trace');
            end
            otherAligned{i} = NaN(nPulses,nSamples);
        end
    end
    
    
    % calculate running standard deviation
    stdTrace = NaN(nPulses,nSamples);
    for s=(1+nNeighbours):(nSamples-nNeighbours)
        stdTrace(:,s) = nanstd(trace(:,(s-2):(s+2)),0,2);
    end
    
    % look for peaks in StdTrace, define pulse as being between the two
    % largest.
    startSamples = NaN(1,nPulses);
    endSamples = NaN(1,nPulses);
    startPeakVals = NaN(1,nPulses);
    endPeakVals = NaN(1,nPulses);
    for p=1:nPulses
        [pks,locs] = findpeaks(stdTrace(p,:),'SORTSTR','descend'); % returns peaks in stdTrace from largest to smallest
        
        if (locs(1)<locs(2)) % startSample is the largest peak
            startSamples(p) = locs(1);
            endSamples(p) = locs(2);
            startPeakVals(p) = pks(1);
            endPeakVals(p) = pks(2);
        else % startSample is the second largest peak
            startSamples(p) = locs(2);
            endSamples(p) = locs(1);
            startPeakVals(p) = pks(2);
            endPeakVals(p) = pks(1);            
        end
        
    end
    
    noiseThreshold = nStd.*nanmean(stdTrace(:)); % require peak amplitudes to be larger than this
    foundPeak = (abs(startPeakVals)>noiseThreshold) | (abs(endPeakVals)>noiseThreshold); % =1 if a peak exists above threshold for that pulse.
    
    % calculate mean sample range of the pulse
    meanLength = round(nanmean(endSamples(foundPeak)-startSamples(foundPeak)));
    pulseRange = (nBaseline+1):(nBaseline+1+meanLength); % pulse starts at nBaseline+1 after alignment  
    
    aligned = NaN(nPulses,nSamples);
    for p=1:nPulses
        
        if (foundPeak) % leave any pulse that didn't go above threshold as NaNs.
            
            if (startSamples(p)>nBaseline) % if startSample>nBaseline, we can include real noise before the pulse in aligned
                aligned(p,1:nBaseline) = trace(p, (startSamples(p)-nBaseline):(startSamples(p)-1));
            else % if startSample<nBaseline, must leave some NaNs at the beginning of aligned
                aligned(p,(nBaseline-startSamples(p)+2):nBaseline) = trace(p, 1:(startSamples(p)-1));
            end

            nPointsAfterStart = nSamples-startSamples(p)+1; % number of samples left after the found start of the pulse

            if (nPointsAfterStart < nPointsAfterBaseline) %
                aligned(p, (nBaseline+1):(nBaseline+nPointsAfterStart)) = trace(p, startSamples(p):end);
            else
                aligned(p,(nBaseline+1):end) = trace(p, startSamples(p):(startSamples(p)+nPointsAfterBaseline-1));
            end
            
            if (~isempty(otherTraces))
            
                for q=1:nOtherTraces
                
                    if (startSamples(p)>nBaseline) % if startSample>nBaseline, we can include real noise before the pulse in aligned
                        otherAligned{q}(p,1:nBaseline) = otherTraces{q}(p, (startSamples(p)-nBaseline):(startSamples(p)-1));
                    else % if startSample<nBaseline, must leave some NaNs at the beginning of aligned
                        otherAligned{q}(p,(nBaseline-startSamples(p)+2):nBaseline) = otherTraces{q}(p, 1:(startSamples(p)-1));
                    end

                    nPointsAfterStart = nSamples-startSamples(p)+1; % number of samples left after the found start of the pulse

                    if (nPointsAfterStart < nPointsAfterBaseline) %
                        otherAligned{q}(p, (nBaseline+1):(nBaseline+nPointsAfterStart)) = otherTraces{q}(p, startSamples(p):end);
                    else
                        otherAligned{q}(p,(nBaseline+1):end) = otherTraces{q}(p, startSamples(p):(startSamples(p)+nPointsAfterBaseline-1));
                    end

                    
                end
                
            end
            
        end
        
    end
    

    

end

