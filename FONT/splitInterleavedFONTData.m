function [ FONTDataFFOn, FONTDataFFOff ] = splitInterleavedFONTData( FONTData )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%     FONTDataOdd = FONTData;
%     FONTDataEven = FONTData;
%     
%     FONTDataOdd.ADCs = FONTDataOdd.ADCs(:,1:2:end,:);
%     FONTDataOdd.DACs = FONTDataOdd.DACs(1:2:end,:,:);
%     FONTDataOdd.BoardStamp = FONTDataOdd.BoardStamp(:,1:2:end);
%     FONTDataOdd.rbs = FONTDataOdd.rbs(1:2:end,:);
%     FONTDataOdd.Time = FONTDataOdd.Time(1:2:end,:);
%     
%     FONTDataEven.ADCs = FONTDataEven.ADCs(:,2:2:end,:);
%     FONTDataEven.DACs = FONTDataEven.DACs(2:2:end,:,:);
%     FONTDataEven.BoardStamp = FONTDataEven.BoardStamp(:,2:2:end);
%     FONTDataEven.rbs = FONTDataEven.rbs(2:2:end,:);
%     FONTDataEven.Time = FONTDataEven.Time(2:2:end,:);
    
    try
        FFStatus = FONTData.FFStatus;
    catch
        error('FFStatus variable not found. Try reprocessing FONT data')
    end

    FONTDataFFOn = FONTData;
    FONTDataFFOff = FONTData;

    FONTDataFFOn.ADCs = FONTDataFFOn.ADCs(:,(FFStatus==1),:);
    FONTDataFFOn.DACs = FONTDataFFOn.DACs(:,(FFStatus==1),:);
    FONTDataFFOn.BoardStamp = FONTDataFFOn.BoardStamp(:,(FFStatus==1));
    FONTDataFFOn.rbs = FONTDataFFOn.rbs((FFStatus==1),:);
    FONTDataFFOn.Time = FONTDataFFOn.Time((FFStatus==1),:);
    FONTDataFFOn.ChanOffset = FONTDataFFOn.ChanOffset((FFStatus==1),:);
    FONTDataFFOn.FFStatus = FONTDataFFOn.FFStatus((FFStatus==1));

    FONTDataFFOff.ADCs = FONTDataFFOff.ADCs(:,(FFStatus==0),:);
    FONTDataFFOff.DACs = FONTDataFFOff.DACs(:,(FFStatus==0),:);
    FONTDataFFOff.BoardStamp = FONTDataFFOff.BoardStamp(:,(FFStatus==0));
    FONTDataFFOff.rbs = FONTDataFFOff.rbs((FFStatus==0),:);
    FONTDataFFOff.Time = FONTDataFFOff.Time((FFStatus==0),:);
    FONTDataFFOff.ChanOffset = FONTDataFFOff.ChanOffset((FFStatus==0),:);
    FONTDataFFOff.FFStatus = FONTDataFFOff.FFStatus((FFStatus==0));

end

