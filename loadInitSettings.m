function [frascatiTimeStamp,devNames,devValues] = loadInitSettings(filename)
    delimiter = ',';
    formatSpec = '%s%s%s%s%[^\n\r]';
    fileID = fopen(filename,'r');
    
    if (fileID == -1) % file not found
        warning('%s not found. Merged data will not contain calibration or init settings.',filename);
        frascatiTimeStamp = '';
        devNames = {};
        devValues = [];
    else
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
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
        rawData = dataArray{4};
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
                    numericData(row, 4) = numbers{1};
                    raw{row, 4} = numbers{1};
                end
            catch me
            end
        end

        %% Split data into numeric and cell columns.
        rawNumericColumns = raw(:, 4);
        rawCellColumns = raw(:, [1,2,3]);

        %% Replace non-numeric cells with NaN
        R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
        rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

        %% Allocate imported array to column variable names
        devNames = rawCellColumns(1:end-1, 1);
        devValues = cell2mat(rawNumericColumns(1:end-1, 1));
        frascatiTimeStamp = rawCellColumns(end, 2);
        frascatiTimeStamp = frascatiTimeStamp{1};
    end
end
