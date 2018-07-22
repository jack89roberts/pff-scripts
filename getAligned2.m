function [aligned,otherAligned,pulseRange] = getAligned2(trace,np_rms,npoints, preserveOrder, otherTraces)
%% getAligned2: Takes an array of pulses and aligns them in time (so that the start of the pulse is always at the same sample number).
% trace(pulses,samples): 2d array of signals
% np_rms: no. points used for baseline noise calculation
% npoints: no. points saved to aligned data
% preserveOrder - Any pulse in trace that doesn't go outside threshold will return a
% NaN pulse, so that the size of aligned matches the size of trace.
% otherTraces: apply the same alignment found for trace to the signals in
% otherTraces (must be cell array and must have same no. pulses as trace).
% pulseRange: estimated sample range of pulse (startSample:endSample)
nlead = 20;% can not be smaller then 2, no. points to include before identified pulse start
nStds = 15; % start of signal where output goes above this many standard deviations % Jack - I found that 5 sigma sometimes caught unwanted noise in baseline so tried 10 instead... even 15!

[n,m] = size(trace);
alignedH2 = NaN*ones(n,npoints); % Jack changed to NaNs from zeros so out of range output is NaN and not 0.

if nargin>3 % preserveOrder given
    
    if nargin>4 % otherTraces given
        if (~iscell(otherTraces))
            error('otherTraces must be a cell array');
        else        
            nOtherSignals = length(otherTraces);
            otherAligned = cell(1,nOtherSignals);
            for i=1:nOtherSignals
                [nRows,~] = size(otherTraces{i});
                if (nRows ~= n)
                    error('otherTraces must have same no. pulses as trace');
                end
                otherAligned{i} = NaN*ones(n,npoints);
            end
        end     
    else
        otherAligned = {};
    end       
else
    preserveOrder = 0;
    otherAligned = {};
end

lintrace = reshape(trace(:,1:np_rms),1,n*np_rms); % first np_rms points for all pulses (get noise baseline)
std5 = nanstd(lintrace)*nStds;%nanstd(lintrace)*5; % std of noise 
mean5 = nanmean(lintrace); % mean of noise

% Jack - check noise at end as well just in case signal starts in first np_rms samples
% Jack - commented this as it causes issues with reflection on 3rd monitor
% diode.
% lintraceEnd = reshape(trace(:,m-np_rms+1:m),1,n*np_rms); 
% std5End = nanstd(lintraceEnd)*nStds;
% mean5End = nanmean(lintraceEnd); 
% if (std5End < 0.9*std5) % Jack - fudge factor of 0.9 to deal with reflection on diode mon3 (signal doesn't go back to 0 after end of pulse)
%     std5 = std5End;
%     mean5 = mean5End;
% end

nal = 0;
for i=1:n % loop over pulses
 % find first point where  output more than 5 std away from noise
 % ix is index of first point for aligned data
 ixP = find(trace(i,:) < mean5 - std5, 1, 'first');
 ixM = find(trace(i,:) > mean5 + std5, 1, 'first');

 if (~isempty(ixP) || ~isempty(ixM))
   nal = nal + 1; % no. pulses found with output outside 5 std
   
   if (isempty(ixP)) 
     ix = ixM; 
   elseif (isempty(ixM)) 
     ix = ixP;
   else
     ix = min(ixP, ixM);
   end
   
   ex = ix+npoints-nlead-1;
   if (ex > m) %  Jack added to deal with index exceeds dimensions error when ix too close to nSamples (m)
       ex = m;
   end
   
%    if (ix <= nlead)
%        ix = nlead+1;
%    end 
   if (ix<= nlead) % Jack - If ix is less than nlead, pad beginning of alignedH2 with NaNs so that start of pulse is always at sample nlead+1. (Rather than moving start of pulse closer to sample 1)
       nNaNs = nlead-ix+1;
       alignedH2(nal,(nNaNs+1):npoints) =  trace(i,1:ex);
   else
       alignedH2(nal,1:length(ix-nlead:ex)) =  trace(i,ix-nlead:ex); % save data in found range to aligned array
   end
       
   %fprintf(1,'%d %d %d  %f %f\n',i,ix, npoints, mean5, std5);
    
   if (~isempty(otherAligned)) % Jack - apply alignment from trace to otherTraces
       for j=1:nOtherSignals
            ex = ix+npoints-nlead-1;
            [~,nColumns] = size(otherTraces{j});
            if (ex > nColumns)
                ex = nColumns;
            end
                    
            if (ix<= nlead)
                otherAligned{j}(nal,(nNaNs+1):npoints) =  otherTraces{j}(i,1:ex);
            else
                otherAligned{j}(nal,1:length(ix-nlead:ex)) =  otherTraces{j}(i,ix-nlead:ex); % save data in found range to aligned array
            end

       end
   end
   
 else % signal never went outside 5 std
   %fprintf(1,'%d No below threshold %f %f\n',i, mean5, std5);
   
   % preserve order of pulses
   if (preserveOrder)
       nal = nal + 1; 
       
       alignedH2(nal,:) =  NaN; % save data in found range to aligned array

       if (~isempty(otherAligned))
           for j=1:nOtherSignals
                otherAligned{j}(nal,:) = NaN;
           end
       end
   end
   
 end
 
end

aligned = alignedH2(1:nal,1:npoints); % remove unfilled elements (if output of some pulses didn't go outside 5 std)

if (~isempty(otherAligned))
    for i=1:nOtherSignals
        otherAligned{i} = otherAligned{i}(1:nal,1:npoints);
    end
end

% Jack - estimate pulse range (often unreliable, particularly for end of
% the pulse in cases where the rising and falling edges are not symmetric.
pulseStarts = NaN(1,n);
pulseEnds = NaN(1,n);
for i=1:n
    tmpPulseRegion = abs(aligned(n,:)) > mean5+std5; % 1 for all values above noise threshold, 0 otherwise
    tmpEdges = diff(tmpPulseRegion); % 1 for rising edge, -1 for falling edge. 0 everywhere else.
    
    if (sum(tmpEdges==1)==1) && (sum(tmpEdges==-1)==1) % check only one region above noise threshold was found
        pulseStarts(n) = find(tmpEdges==1) + 1; % first point in pulse above noise threshold is location of max peak in diff +1
        pulseEnds(n) = find(tmpEdges==-1); % last point in pulse above noise threshold is location of min peak in diff
    end
end

meanStart = round(nanmean(pulseStarts));
meanEnd = round(nanmean(pulseEnds));
pulseRange = meanStart:meanEnd;

end
