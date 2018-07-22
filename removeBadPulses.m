function [ goodSignal, goodOtherSignals ] = removeBadPulses(signal, sampleRange, otherSignals, nSigma, showResults)
%removeBadPulses Sets any pulses nSigma away from mean in signal to NaN.
%   signal: Traces to check for bad pulses, 2d array (pulses, samples).
%
%   sampleRange: calculate means/stds across these samples when deciding
%   which pulses to remove.
%
%   otherSignals: Sets same pulses to NaN as found for signal if given (otherSignals
%   MUST be a cell array and same no. pulses as signal).
%
%   nSigma: Remove (set to NaN) any pulses nSigma away from mean. Default =
%   3 if not given.
%
%   showResults: 1: print a message with no. pulses removed, 2: make a plot
% 
%   goodSignal: signal with bad pulses set to NaN
%
%   goodOtherSignals: otherSignals with bad pulses set to NaN
    [nPulses,nSamples] = size(signal);

    if (nargin <5)
        showResults = 0;
    end
    if (nargin < 4)
        nSigma = 3;  
    end
    if (nargin < 3)
        otherSignals = {};
    end
    if (nargin < 2)
        sampleRange = 1:nSamples;
    end
    
    if (nPulses == 1 || nSamples==1)
        meanSignals = signal;
        nSamples = 1;
    else
        meanSignals = nanmean(signal(:,sampleRange),2);
    end
    overallMean = nanmean(meanSignals);
    overallStd = nanstd(meanSignals);
    remAbove = overallMean + nSigma.*overallStd;
    remBelow = overallMean - nSigma.*overallStd;
    
    isBad = (meanSignals>remAbove) | (meanSignals<remBelow);
    
    goodSignal = signal;
    if (nPulses==1)
        goodSignal(isBad) = NaN;
    else
        goodSignal(isBad,:) = NaN;
    end
    if (~isempty(otherSignals))
        goodOtherSignals = otherSignals;
        for i=1:length(otherSignals)
            if (nSamples==1)
                goodOtherSignals{i}(isBad) = NaN;
            else
                goodOtherSignals{i}(isBad,:) = NaN;
            end
        end
    end
    
    if (showResults==1)
        fprintf('%d pulses out of %d removed.\n',sum(isBad),sum(~isnan(meanSignals)));
    elseif (showResults==2)
        figure;
        plot(meanSignals);
        hold all;
        xLim = [1 nPulses];
        plot(xLim,[remAbove remAbove],'k','LineWidth',3);
        plot(xLim,[remBelow remBelow],'k','LineWidth',3);
        yLim = get(gca,'YLim');
        badInds = find(isBad);
        for b=1:length(badInds)
            plot([badInds(b) badInds(b)],yLim,'r')
        end
        xlabel('Pulse No.')
        ylabel('Mean Pulse')
        title(sprintf('%d pulses out of %d removed.\n',sum(isBad),sum(~isnan(meanSignals))))
        ylim(yLim);
        xlim(xLim);
    end
end

