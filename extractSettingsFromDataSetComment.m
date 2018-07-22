function [ settingHeadings, pulseSettings ] = extractSettingsFromDataSetComment( CTFData )
%extractSettingsFromDataSetComment Extracts settings saved to data set
%comment in comma separated format.
%
%   comment must have even number of rows of form: 
%   settingName1, value1, settingName2, value2, settingName3, value3 etc.
%
%   settingHeadings: cellArray of setting names
%
%   pulseSettings: array with size (nSettings,nPulses) containing setting
%   values (in same order as settingHeadings).

    % use first pulse to get format of data set comment and make arrays,
    % extract names etc.
    pulseComment = CTFData(1).comment;
    splitComment = regexp(pulseComment, ',', 'split');
    splitComment = strtrim(splitComment); % remove spaces

    nColumns = length(splitComment);
    nSettings = nColumns./2; % each setting has heading and value

    settingHeadings = cell(1,nSettings);
    for i=1:2:nColumns
        settingHeadings{ceil(i/2)} = splitComment{i};
    end
    
    nPulses = length(CTFData);
    pulseSettings = NaN*ones(nSettings, nPulses);

    % loop over pulses and extract settings

    for t=1:nPulses

        % take comment of form: Name, value and put value in kickSetting
        pulseComment = CTFData(t).comment;
        splitComment = regexp(pulseComment, ',', 'split');

        for i=2:2:nColumns
            pulseSettings(i/2, t) = str2double(splitComment{i});
        end

    end

end