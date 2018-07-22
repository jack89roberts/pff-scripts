function [ subFONTData ] = extractSubsetFONTData( FONTData, pulseRange )
%extractSubsetFONTData Extracts pulses in pulseRange from FONTData
%   Detailed explanation goes here


    subFONTData = FONTData;

    subFONTData.ADCs = subFONTData.ADCs(:,pulseRange,:);
    subFONTData.DACs = subFONTData.DACs(:,pulseRange,:);
    subFONTData.BoardStamp = subFONTData.BoardStamp(:,pulseRange);
    subFONTData.rbs = subFONTData.rbs(pulseRange,:);
    subFONTData.Time = subFONTData.Time(pulseRange,:);
    subFONTData.ChanOffset = subFONTData.ChanOffset(pulseRange,:);
    subFONTData.FFStatus = subFONTData.FFStatus(pulseRange);


end

