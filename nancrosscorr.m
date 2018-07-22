function [ xcf, lags ] = nancrosscorr(y1, y2, lags)
%nancrosscorr Attempt to make a cross correlation function that can cope
%with NaN values.
%   y1, y2: Input signals to check cross correlation. 2 dimenional:
%   nPulses, nSamples.
%
%   xcf: cross correlation values for each lag (delay)
%
%   lags: lag (y2 wrt y1) used for each index of xcf (same as input
%   argument if given, or -10:10 if not given).
%
%   Normalisation may not give coefficient values between -1 and 1 (to be
%   checked), but magnitude of xcf should at least give a good indication
%   of best delay between y1 and y2.
    
    % initialise some arrays and values
    %tic
    if (nargin < 3)
        lags = -10:10;
    end
    
    [nPulses,nSamples] = size(y1);
    
    xcf = zeros(1,length(lags));
    nGood = zeros(1,length(lags));
    
    % calculate means and normalisation
    if (nPulses==1)
        % calculate cross correlation between samples instead
        meanY1 = nanmean(y1);
        meanY2 = nanmean(y2);
        nPulses = nSamples;
    else
        meanY1 = nanmean(nanmean(y1,2));
        meanY2 = nanmean(nanmean(y2,2));
    end
    
    normY1 = 0;
    normY2 = 0;
    nGoodY1 = 0;
    nGoodY2 = 0;
    
    for p=1:nPulses 
        pulseY1 = y1(p);          
        isGoodY1 = ~(sum(isnan(pulseY1)) == length(pulseY1));            
        if (isGoodY1)
            nGoodY1 = nGoodY1 + 1;
            normY1 = normY1 + (nanmean(pulseY1)-meanY1).^2;       
        end
        
        pulseY2 = y2(p);
        isGoodY2 = ~(sum(isnan(pulseY2)) == length(pulseY2));
        if (isGoodY2)
            nGoodY2 = nGoodY2 + 1;
            normY2 = normY2 + (nanmean(pulseY2)-meanY2).^2;
        end
        
    end
    normY1 = normY1./nGoodY1;
    normY2 = normY2./nGoodY2;
    norm = sqrt(normY1.*normY2);
    
    % calculate cross correlation for each lag
    for l = 1:length(lags)
        %l
        lagVal = lags(l);
        
        if(lagVal < 0)            
            y2Range = 1:(nPulses+lagVal);
            y1Range = (1-lagVal):nPulses;
        else
            y2Range = (1+lagVal):nPulses;
            y1Range = 1:(nPulses-lagVal);
        end
        for p = 1:length(y1Range)
            
            pulseY1 = y1(y1Range(p)); 
            pulseY2 = y2(y2Range(p));
            
            % only calculate when both y1 and y2 have non-NaN values
            isGoodY1 = ~(sum(isnan(pulseY1)) == length(pulseY1));            
            isGoodY2 = ~(sum(isnan(pulseY2)) == length(pulseY2));
            
            if (isGoodY1 && isGoodY2)
                nGood(l) = nGood(l) + 1;
                
                meanPulseY1 = nanmean(pulseY1);
                meanPulseY2 = nanmean(pulseY2);  
                
                xcf(l) = xcf(l) + (meanPulseY1-meanY1).*(meanPulseY2-meanY2);
            end
            
        end
    end
    
    xcf = xcf./(nGood.*norm);
    %toc
end

