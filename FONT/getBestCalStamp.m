function [ calStamp,calDir ] = getBestCalStamp( dataSetName )
%getBestCalStamp Looks for the time stamp of the last calibration before dataSetName

        calDir = ['/home/jack/PhaseFeedforward/FONTData/' dataSetName(1:6) '/Calibration/']; % first guess for calDir - calibration directory for the same month as dataSet

        calDirList = dir(calDir);
        calDirList = calDirList(3:end); % first two should be '.' and '..'
        nCals = length(calDirList);
        
        
        if nCals>0 % found some calibrations in the same month
        
            calNames = cell(1,nCals);
            calDate = NaN(1,nCals);
            calTime = NaN(1,nCals);
            for i=1:nCals
                calNames{i} = calDirList(i).name;
                calSplitStr = strsplit(calNames{i},'_');
                calDate(i) = str2double(calSplitStr{1});
                calTime(i) = str2double(calSplitStr{2});
            end

            dataDate = dataSetName(1:8);
            dataTime = dataSetName(10:13);
            dataDate = str2double(dataDate);
            dataTime = str2double(dataTime);

            diffDate = calDate-dataDate;
            diffTime = calTime-dataTime;

            diffDate(diffDate>0) = NaN; % only want to check calibrations before data set

            if (sum(isnan(diffDate))~=length(diffDate)) % if found some calibrations before dataSet in the same month

                [minDiffDate,minDiffDateInd] = min(abs(diffDate));

                if (sum(minDiffDate==diffDate)==1) % only one calibration on the most appropriate day
                    calStamp = calNames{minDiffDateInd};

                else % more than one calibration on the most appropriate day

                    if (minDiffDate==0) % calibration is on the same day as dataSetName. Try to find a calibration before, but if not just use one at a later time on same day 
                        onlyTimesBefore = diffTime;
                        onlyTimesBefore(onlyTimesBefore>0) = NaN;

                        if (sum(isnan(onlyTimesBefore))==length(onlyTimesBefore)) % no time before dataSet, just use the closest calibration after
                            [~,minDiffTimeInd] = min(abs(diffTime));
                            calStamp = calNames{minDiffTimeInd};
                        else % is a calibration before on the same day, so use that one
                            [~,minDiffTimeInd] = min(abs(onlyTimesBefore));
                            calStamp = calNames{minDiffTimeInd};
                        end

                    else % calibration if on a different day. Just use the last calibration on that day
                        diffDate(diffDate ~= minDiffDate) = NaN;
                        bestInd = find(diffDate,1,'last');
                        calStamp = calNames{bestInd};
                    end

                end

            else % no calibrations before this one in the same month, look in the month before and use the last calibration        
                calDir = sprintf('/home/jack/PhaseFeedforward/FONTData/%d/Calibration', str2double(dataSetName(1:6))-1); % first guess for calDir - calibration directory for the same month as dataSet
                calDirList = dir(calDir);
                calDirList = calDirList(3:end); % first two should be '.' and '..'
                nCals = length(calDirList);
                if nCals==0 
                    error('Could not find any calibrations in %s',calDir);
                else
                    calStamp = calDirList(end).name;
                end  
            end
        
        else % no calibrations in the same month as dataSet
            calDir = sprintf('/home/jack/PhaseFeedforward/FONTData/%d/Calibration', str2double(dataSetName(1:6))-1); % first guess for calDir - calibration directory for the same month as dataSet
            calDirList = dir(calDir);
            calDirList = calDirList(3:end); % first two should be '.' and '..'
            nCals = length(calDirList);
            if nCals==0 
                error('Could not find any calibrations in the same month as %s or in the month before (%s)',dataSetName,calDir);
            else
                calStamp = calDirList(end).name;
            end  

        end
end

