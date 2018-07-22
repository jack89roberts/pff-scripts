function [ aligned, otherAligned, pulseRange, errMessage ] = getAlignedXCorr( trace, alignTo, otherTraces, nLead, nStds, plotResult)
%getAlignedXCorr Performs alignment of traces using an initial cross
%correlation step.
%   May produce better results than getAligned2 (less likely to have one
%   sample offsets in output)
%   alignTo:
%       start: moves all pulses to start at sample nLead+1
%       end: moves all pulses to end at nSamples-nLead
%       first: moves all pulses to be aligned with first pulse in trace
  
    %% Define arrays, noise calculation
    if nargin<3
        otherTraces = {};
    end
    if nargin<4
        nLead = 50;% can not be smaller then 2, no. points to use for baseline noise calculation. Also no. points returned before pulse (align to start) or after the pulse (align to end).
    elseif (nLead<2)
        error('nLead must be greater than 2');
    end
    if nargin<5
        nStds = 15; % pulse exists if output goes nStds standard deviations above baseline noise. Pulse start/end based on where output goes nStds/3 above baseline noise.
    end
    if nargin<6
        plotResult = 0;
    end
    
    [nPulses,nSamples] = size(trace);
    tmpAligned = NaN(nPulses,nSamples);

    if nargin>1
        if(~strcmp(alignTo,'start')) && (~strcmp(alignTo,'end') && (~strcmp(alignTo,'first')))
            error('%s is not a valid value for startOrEnd. Use "start", "end" or "first".',alignTo);
        end
    else
        alignTo = 'end';
    end
    
    if nargin>2 % otherTraces given
        if (~iscell(otherTraces))
            error('otherTraces must be a cell array');
        else        
            nOtherTraces = length(otherTraces);
            tmpOtherAligned = cell(1,nOtherTraces);
            for i=1:nOtherTraces
                [nRows,~] = size(otherTraces{i});
                if (nRows ~= nPulses)
                    error('otherTraces must have same no. pulses as trace');
                end
                tmpOtherAligned{i} = NaN(nPulses,nSamples);
            end
        end     
    else
        nOtherTraces = 0;
        tmpOtherAligned = {};
    end
    
    aligned = NaN(nPulses,nSamples);
    otherAligned = cell(1,nOtherTraces);
    if (nOtherTraces > 0)
        for q=1:nOtherTraces
            otherAligned{q} = NaN(nPulses,nSamples);
        end
    end
    pulseRange = [];
    
    errMessage = '';
    
    % calculate baseline noise, define pulse as existing if output above
    % nStds*noise exists.
    %lintrace = reshape(trace(:,1:nLead),1,nPulses*nLead); % first nLead points for all pulses (get noise baseline)
    % check beginning and end and use where jitter is lowest in case we are
    % catching the pulse.
    lintraceStart = trace(:,1:nLead);
    stdNoiseStart = nanstd(lintraceStart(:))*nStds;% std of noise 
    meanNoiseStart = nanmean(lintraceStart(:)); % mean of noise
    lintraceEnd = trace(:,(end-nLead+1):end);
    stdNoiseEnd = nanstd(lintraceEnd(:))*nStds;% std of noise 
    meanNoiseEnd = nanmean(lintraceEnd(:)); % mean of noise
    [stdNoise, lowJitLoc] = min([stdNoiseStart stdNoiseEnd]);
    if (lowJitLoc==1)
        meanNoise = meanNoiseStart;
    else
        meanNoise = meanNoiseEnd;
    end
    noiseMaxThreshold = meanNoise+stdNoise;
    noiseMinThreshold = meanNoise-stdNoise;
    
    pulseMaxOut = max(trace,[],2);
    pulseMinOut = min(trace,[],2);
    foundPulse = (pulseMaxOut > noiseMaxThreshold) | (pulseMinOut<noiseMinThreshold);
    if (sum(foundPulse)==0)
        errMessage = 'No pulses above noise threshold (alignment not attempted).';
        return;
    elseif (sum(foundPulse)<nPulses)
        errMessage = 'Some pulses do not go above noise threshold.';
    end
    
    % find pulses that have NaN sample values (have to use nancrosscorr
    % instead of xcorr if true).
    nanPulse = sum(isnan(trace),2) > 0;
    if (sum(nanPulse)>0)
        errMessage = [errMessage ' NaNs found in trace (results may be less accurate).'];
    end
    
    %% Cross correlation
    % Align all pulses to a reference pulse (first pulse above noiseThreshold and, if
    % possible, without NaNs).
    
    % Extract the reference pulse
    refPulseInd = find(foundPulse&(~nanPulse),1); % first pulse above noise threshold and with no NaNs.
    if (~isempty(refPulseInd))
        refPulse = trace(refPulseInd,:);
        refIsNaN = 0;
    else
        refPulseInd = find(foundPulse,1); % first pulse above noise threshold but with NaNs.
        refPulse = trace(refPulseInd,:);
        refIsNaN = 1;
    end
    
    % calculate cross correlations and align
    if (nPulses>1)
        sentFailedMessage = 0;
        for p=1:nPulses

            if (foundPulse(p)) % if pulse goes above noise threshold

                if (~nanPulse(p) && ~refIsNaN) % if no NaNs use xcorr, which is more robust/faster
                    [c,lags] = xcorr(refPulse,trace(p,:)); % c = correlation values for sample delay in lags 
                else % If NaNs exist use nanxcorr (but note that I've seen this function fail in certain cases. Older version of getAlignedXCorr used my nancrosscorr, which might give better results but is very slow).
                    [c,lags] = nanxcorr(refPulse,trace(p,:));
                end
                [~,maxInd] = max(c); % index of max correlation
                bestLag = lags(maxInd); % delay of max correlation


                % align pulse to the reference (trim points from start if
                % it is too late, add NaNs at start if it is too early).
                if (bestLag<0)
                    tmpAligned(p,1:(end-abs(bestLag))) = trace(p,(abs(bestLag)+1):end);
                elseif (bestLag>0)
                    tmpAligned(p,(1+abs(bestLag)):end) = trace(p,1:(end-abs(bestLag)));
                elseif (bestLag==0)
                    tmpAligned(p,:) = trace(p,:);
                else
                    if (~sentFailedMessage)
                        errMessage = [errMessage ' At least one alignment failed (delay from cross correlation undefined).'];
                        sentFailedMessage = 1;
                    end
                end

                if (nOtherTraces > 0) % perform same alignment on otherTraces
                    for q=1:nOtherTraces
                        if (bestLag<0)
                            tmpOtherAligned{q}(p,1:(end-abs(bestLag))) = otherTraces{q}(p,(abs(bestLag)+1):end);
                        elseif (bestLag>0)
                            tmpOtherAligned{q}(p,(1+abs(bestLag)):end) = otherTraces{q}(p,1:(end-abs(bestLag)));
                        elseif (bestLag==0)
                            tmpOtherAligned{q}(p,:) =  otherTraces{q}(p,:);
                        else
                            if (~sentFailedMessage)
                                errMessage = [errMessage ' At least one alignment failed (delay from cross correlation undefined).'];
                                sentFailedMessage = 1;
                            end
                        end

                    end
                end
            end
        end
    else
        tmpAligned = trace;
        tmpOtherAligned = otherTraces;
    end
    %% Shift all traces so that they start at nLead+1
    % after cross correlation all pulses are aligned but start where
    % reference pulse originally started, not at nLead+1.
    
    % Diff method: define start/end of pulse as where the min/max in diff
    % is (diff = rate of change of trace).
    if (nPulses>1)
        meanPulse = nanmean(tmpAligned);
    else
        meanPulse = trace;
    end
    rawDiffPulse = diff(meanPulse);
    diffPulse = rawDiffPulse;
    diffPulse(abs(diffPulse)<(stdNoise/3)) = NaN; % exclude all points in diff that are not above noiseThreshold
    
    [~,minPeakInd] = min(diffPulse);
    [~,maxPeakInd] = max(diffPulse);
    startPulse = min([minPeakInd maxPeakInd]);
    endPulse = max([minPeakInd maxPeakInd]);
    
    pulseRange = startPulse:endPulse;
    
    if (~strcmp(alignTo,'first'))
    
        if (strcmp(alignTo,'start'))
            shiftPulse = (nLead+1)-startPulse;
        elseif (strcmp(alignTo,'end'))
            shiftPulse = (nSamples-nLead)-endPulse;
        end

        pulseRange = pulseRange + shiftPulse;

    %     Nose threshold method
    %     meanAboveThreshold = abs(meanPulse)>noiseThreshold;
    %     startPulse = find(meanAboveThreshold,1);
    %     shiftPulse = (nLead+1)-startPulse;

        if (shiftPulse>0) % need to move pulse later
            aligned(:,(1+abs(shiftPulse)):end) = tmpAligned(:,1:(end-abs(shiftPulse)));
            if (nOtherTraces > 0)
                for q=1:nOtherTraces
                    otherAligned{q}(:,(1+abs(shiftPulse)):end) = tmpOtherAligned{q}(:,1:(end-abs(shiftPulse)));
                end
            end

        elseif (shiftPulse<0) % need to move pulse earlier
            aligned(:,1:(end-abs(shiftPulse))) = tmpAligned(:,(abs(shiftPulse)+1):end);
            if (nOtherTraces > 0)
                for q=1:nOtherTraces
                    otherAligned{q}(:,1:(end-abs(shiftPulse))) = tmpOtherAligned{q}(:,(abs(shiftPulse)+1):end);
                end
            end
        end
        
    else
        aligned = tmpAligned;
        otherAligned = tmpOtherAligned;
    end
    %% plot the result if requested
    if (plotResult)
        figure;
        subplot(2,2,1)
        plot(trace');
        title('INPUT')
        hold all;
        xLim = get(gca,'XLim');
        plot(xLim,[noiseMaxThreshold noiseMaxThreshold],'k')
        plot(xLim,[noiseMinThreshold noiseMinThreshold],'k')

        subplot(2,2,2)
        plot(tmpAligned');
        title('INITIAL ALIGNMENT')
        hold all;
        yLim = get(gca,'YLim');
        plot(xLim,[noiseMaxThreshold noiseMaxThreshold],'k')
        plot(xLim,[noiseMinThreshold noiseMinThreshold],'k')

        subplot(2,2,3)
        plot(rawDiffPulse);
        title('LOOKING FOR PULSE RANGE')
        hold all;
        xLim = get(gca,'XLim');
        yLim = get(gca,'YLim');
        plot(xLim,[(stdNoise/3) (stdNoise/3)],'k')
        plot(xLim,[-(stdNoise/3) -(stdNoise/3)],'k')
        plot([startPulse startPulse],yLim,'k');
        plot([endPulse endPulse],yLim,'k');
        ylim(yLim);
       
        subplot(2,2,4)
        plot(aligned');
        title('FINAL ALIGNMENT')
        hold all;
        yLim = get(gca,'YLim');
        plot([pulseRange(1) pulseRange(1)],yLim,'k');
        plot([pulseRange(end) pulseRange(end)],yLim,'k');
        ylim(yLim);
    end
    
   
end

