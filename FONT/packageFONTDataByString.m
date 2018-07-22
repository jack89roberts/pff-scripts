function [ fileADCs, fileValues ] = packageFONTDataByString( dataDir, valueString )
%packageFONTDataByString Take a number of extracted FONT data files and
%create the matrix ADCs (ADCind, fileInd, pulseInd, sampleInd). fileValues
%contains the numeric value x by looking for "valueString_x" in the file
%name. Files will be sorted from min to max numeric value, not in the order
%returned by dir.
    
    fileList = dir([dataDir '/*.mat']);
    nFiles = length(fileList);
    
    fileValues = NaN(1,nFiles);
    
    for f=1:nFiles  % extract shifter values
        datName = fileList(f).name;       
        tmpFileValue = strrep(datName,'.mat','');
        tmpFileValue = strsplit(tmpFileValue,'_');
        valueInd = find(strcmp(tmpFileValue,valueString)) + 1; % Assumes value x saved as valueString_x in file names
        tmpFileValue = tmpFileValue{valueInd};
        fileValues(f) = str2double(tmpFileValue);
    end
    
    % sort phase shifter values in order from lowest to highest, then
    % process files in that order.
    [fileValues,fileOrder] = sort(fileValues);
    
    for f=1:nFiles
        fileInd = fileOrder(f);
        datName = fileList(fileInd).name;
        datFile = load([dataDir '/' datName]);

        fieldName = fields(datFile);        
        ADCs = eval(['datFile.' fieldName{1} '.ADCs']);
        %ADCs = permute(ADCs,[1 3 2]); % Now done by default in my version of data extracting script. my dimension order: ADCindex, pulseIndex, sampleIndex
        
        if f==1
            [nADCs,nPulses,nSamples] = size(ADCs);
            fileADCs = NaN(nADCs,nFiles,nPulses,nSamples);
        end
        
        fileADCs(:,f,:,:) = ADCs;
    
    end
    

end

