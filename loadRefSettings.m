function [deviceString, value, time] = loadRefSettings(dirName)

    fileList = dir([dirName '/*ref_setting']);
    fileList = {fileList.name};
    
    nFiles = length(fileList);
    
    % sort in to time order (not immediate as times before 10:00 are
    % shorter due to no leading 0).
    time = NaN(nFiles,1);
    for i=1:nFiles
        timeStr = fileList{i};
        if (length(fileList{i})==30)
            timeStr = timeStr(1:18);
        elseif (length(fileList{i})==31)
            timeStr = timeStr(1:19);
        else
            error('incorrect time stamp length')
        end
%         timeStr = timeStr(1:19);
        time(i) = datenum(timeStr);
    end
        
    [time,fileOrder] = sort(time);
    
    for i=1:nFiles
        % load file
        tmpFile = [dirName '/' fileList{fileOrder(i)}];
        [tmpDevs,tmpValues] = loadRefSettingFile(tmpFile);
        
        % initialise arrays
        if (i==1)
            deviceString = tmpDevs;
            nDevices = length(deviceString);
            value = NaN(nFiles,nDevices);
        end
        
        if (~isequal(tmpDevs,deviceString))
            % not implemented anything to deal with a file being different to the first one - NaN values for that file
            
        else
            value(i,:) = tmpValues;
        end
        
%         % compare device contents of this file to current arrays, save data
%         % to arrays.
%         if length(tmpDevs)==nDevices % new data has the same number of devices
%             
%             if (isequal(tmpDevs,deviceString)) 
%                 % new data has same devices in same order: save new data
%                 % directly to value array.
%                 value(i,:) = tmpValues;
%             
%             else
%                 % new data has same number of devices, but in different
%                 % order or devices included are different: need to make
%                 % changes to arrays.
%                 newDevices = ~ismember(tmpDevs,deviceString);
%                 
%                 if (sum(newDevices) == 0)
%                 
%                 else
%                     
%                 end
%                 
%             end
%             
%         else
%             
%         end
               
        
        
    end

end
