function signal = extractCTFSignalFromMergedData( signalName, data )
% extractCTFSignalFromMergedData Takes the array of structs data
% containing matlab monitor output data (output from mergeMatMonData) and
% extracts the signal with name signalName in to the array signal.
%  if signal is numeric: signal(triggerIndex, sampleIndex)
%  if signal is struct: signal(triggerIndex)
%  if signal is string: signal{triggerIndex}

    nTriggers = length(data);

    firstSignal = eval(['data(1).' signalName]);

    if (isnumeric(firstSignal))
        nSamples = eval(sprintf('length(data(1).%s)',signalName));        

        % create signal array (numeric array)
        try
            for t=1:nTriggers
                %display(['file no ' num2str(t)]);
                arr = eval(['data(t).' signalName]);
                %t          size(arr)
                signal(t,:) = arr;
            end
        catch % if have some problem, try skipping the first 2 files.
            clear signal; 
            for t=3:nTriggers
                %display(['file no ' num2str(t)]);
                arr = eval(['data(t).' signalName]);
                %t          size(arr)
                signal(t,:) = arr;
            end    
        end

    elseif (isstruct(firstSignal))

        % create signal array (array of structs)
        % NB - look in to pre-allocating empty array of structs for speed
        for t=1:nTriggers
            signal{t} = eval(['data(t).' signalName]);
        end

    elseif (ischar(firstSignal))
        signal = cell(1,nTriggers);

        % create signal array (cell array of strings)
        for t=1:nTriggers
            signal{t} = eval(['data(t).' signalName]);
        end

    else
        error('Input signal data type is not numeric, char or struct. ');
    end

 
end