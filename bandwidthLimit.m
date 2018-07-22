function [ filtData ] = bandwidthLimit(data,maxFreq,sampleFreq,order )
%bandwidthLimit Uses a butterworth filter to remove the frequency
%components of data above maxFreq.
%   data: data to be filtered. If 1D treats as 1 pulse, if 2D treats first
%   dimension as pulse index and second dimension as sample no.
%   maxFreq: cut-off frequency for low-pass filter (i.e. bandwidth limit)
%   sampleFreq: sample frequency of data
%   order: order of filter to be used (optional, default=6)
%   filtData: the filtered data
    
    if nargin<4
        order = 6;
    end
    
    normFreq = maxFreq/(sampleFreq/2);
    [b,a] = butter(order,normFreq);

    [nPulses,nSamples] = size(data);
    
    if (nPulses==1 || nSamples==1) % just 1 pulse
        filtData = filtfilt(b,a,data);
    
    else % loop over pulses
        filtData = NaN(nPulses,nSamples);
        
        for p=1:nPulses
            filtData(p,:) = filtfilt(b,a,data(p,:));
        end
    end
    
end

