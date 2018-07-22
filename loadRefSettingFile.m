function [deviceStr, value] = loadRefSettingFile(filename)
% Loads data from a .ref_setting file from the CTF3 monitor.
% Files are found in /user/ctf3op/CTF3Monitor/referenceDirectory
% Majority of this script generated automatically by MatLab.

% duplicate devices:
%   if all same value: use value
%   if one non-zero value only: use non-zero value
%   if multiple non-zero values: NaN

%% Initialize variables.
%tic


delimiter = ' ';
startRow = 2;
endRow = inf;


%% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric strings to numbers.
% Replace non-numeric strings with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

% Converts strings in the input cell array to numbers. Replaced non-numeric
% strings with NaN.
rawData = dataArray{2};
for row=1:size(rawData, 1);
    % Create a regular expression to detect and remove non-numeric prefixes and
    % suffixes.
    regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
    try
        result = regexp(rawData{row}, regexstr, 'names');
        numbers = result.numbers;
        
        % Detected commas in non-thousand locations.
        invalidThousandsSeparator = false;
        if any(numbers==',');
            thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
            if isempty(regexp(thousandsRegExp, ',', 'once'));
                numbers = NaN;
                invalidThousandsSeparator = true;
            end
        end
        % Convert numeric strings to numbers.
        if ~invalidThousandsSeparator;
            numbers = textscan(strrep(numbers, ',', ''), '%f');
            numericData(row, 2) = numbers{1};
            raw{row, 2} = numbers{1};
        end
    catch me
    end
end


%% Split data into numeric and cell columns.
rawNumericColumns = raw(:, 2);
rawCellColumns = raw(:, 1);


%% Allocate imported array to column variable names
deviceStr = rawCellColumns(:, 1);
value = cell2mat(rawNumericColumns(:, 1));
% 
% toc
% %% Split device strings in to device name and acquired parameter
% nDevices = length(deviceStr);
% deviceName = cell(1,nDevices);
% acquiredParameter = cell(1,nDevices);
% 
% for d=1:nDevices
%     tmpDevStr = deviceStr{d};
%     tmpDevStr = strsplit(tmpDevStr,'/');
%     deviceName{d} = tmpDevStr{1};
%     acquiredParameter{d} = tmpDevStr{2};
% end
% 
% toc

%% Deal with duplicate values
[uniqueDevs,uniqueInds] = unique(deviceStr);
uniqueValues = value(uniqueInds); % get array of values with unique devices in same order
duplicateInds = ~ismember(1:length(deviceStr),uniqueInds); % indices that don't appear in list of unique devices - they're the duplicates
duplicateDevs = unique(deviceStr(duplicateInds)); % get the names of duplicate devices (unique needed to deal with cases where device appears more than twice)

for i=1:length(duplicateDevs)
    devInds = strcmp(deviceStr,duplicateDevs{i});
    devValues = value(devInds);
    indForUnique = strcmp(duplicateDevs{i},uniqueDevs);
    
    if (length(unique(devValues))==1) 
        % all duplicates have the same value, don't need to do anything (already in uniqueValues)       
    elseif (sum(devValues~=0)==1)
         % only one non-zero value, assume that is the correct one
         uniqueValues(indForUnique) = devValues(devValues~=0);
    else
        % otherwise just set it to NaN
        uniqueValues(indForUnique) = NaN;
        fprintf('Undefined value: %s\n',duplicateDevs{i});
    end
end

deviceStr = uniqueDevs;
value = uniqueValues;

%toc
end