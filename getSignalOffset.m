function [offsetInSig1Samples, maxCrossCorrel] = getSignalOffset( sig1Data, sig1SampRate, sig2Data, sig2SampRate, plotResult )
%getSignalOffset 
%Returns the time offset that gives the best cross correlation
%between two signals which can have different sampling rates.
%   sig1Data/sig2Data : signals to cross-correlate
%   sig1SampRate/sig2SampRate: time interval between samples for each
%   signal.
    if nargin < 5
        plotResult = false;
    end

    sig1Samples = length(sig1Data);
    sig2Samples = length(sig2Data);
    
    sig1Time = sig1SampRate*(0:(sig1Samples-1));
    sig2Time = sig2SampRate*(0:(sig2Samples-1));
    
    if ( sig1SampRate ~= sig2SampRate ) 
        % if the sampling rates are different need to do some interpolation
        % to get equivalent sampling points to cross correlate
        
        sig2Data = spline(sig2Time,sig2Data, linspace(0,max(sig2Time),sig1Samples));
        sig2Time = linspace(0,max(sig2Time),sig1Samples);
        
    end

    % cross correlate the signals
    [c, lags] = xcorr(sig1Data, sig2Data);
    [maxCrossCorrel,maxCorrelIndex] = max(c);
    offsetInSig1Samples = lags(maxCorrelIndex);
    
    
    % TEST
%     diffSq = NaN*ones(1,sig1Samples); 
%     for i = 1:sig1Samples
%         sig1Shifted = [sig1Data zeros(1, i)];
%         sig2Shifted = [zeros(1, i) sig2Data];
%         diffSq(i) = sum((sig1Shifted - sig2Shifted).^2);
%     end
% 
%     [maxCrossCorrel, offsetInSig1Samples] = min(diffSq);
    
    if (plotResult)
        figure;
        plot(sig1Time,sig1Data);
        hold all;
        plot(sig2Time,sig2Data);
        plot(sig2Time+offsetInSig1Samples,sig2Data);

        if ( sig1SampRate ~= sig2SampRate ) 
            legend('Signal 1','Signal 2 (Interpolated): Original ', 'Signal 2 (Interpolated): Offset');
        else
            legend('Signal 1','Signal 2: Original', 'Signal 2: Offset');
        end
        
    end
    
end

