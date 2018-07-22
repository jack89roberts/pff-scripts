function [ mixers, diodes, phaseShiftValues ] = packageFONTCalData( dataDir, mixerADCs, diodeIsPresent, diodeADCs )
%packageFONTCalData Take extracted FONT data and use it to create arrays in
%the same format as needed for calibration scripts.
    
    if (nargin<2)
        mixerADCs = [];
    end
    if (nargin<3)
        diodeIsPresent = 1;
    end
    if (nargin<4)
        diodeADCs = [];
    end

    calFiles = dir([dataDir '/*.mat']);
    nCalFiles = length(calFiles);
    
    phaseShiftValues = NaN(1,nCalFiles);
    
    for f=1:nCalFiles  % extract shifter values
        datName = calFiles(f).name;       
        tmpShiftValue = strrep(datName,'.mat','');
        tmpShiftValue = strsplit(tmpShiftValue,'_');
        shftValInd = find(strcmp(tmpShiftValue,'Shft')) + 1; % Assumes shifter value x saved as Shft_x in file names
        tmpShiftValue = tmpShiftValue{shftValInd};
        phaseShiftValues(f) = str2double(tmpShiftValue);
    end
    
    % sort phase shifter values in order from lowest to highest, then
    % process files in that order.
    [phaseShiftValues,fileOrder] = sort(phaseShiftValues);
    
    for f=1:nCalFiles
        fileInd = fileOrder(f);
        datName = calFiles(fileInd).name;
        datFile = load([dataDir '/' datName]);
        
        fieldName = fields(datFile);        
        ADCs = eval(['datFile.' fieldName{1} '.ADCs']);
        [tmpMixers,tmpDiodes] = extractFONTMixerDiode(ADCs,1,mixerADCs,diodeIsPresent,diodeADCs);
        
        if f==1
            [nMons,nPulses,nSamples] = size(tmpMixers);
            mixers = NaN(nMons,nCalFiles,nPulses,nSamples);
            if (isempty(tmpDiodes))
                diodes = [];
            else
                diodes = NaN(nMons,nCalFiles,nPulses,nSamples);
            end
        end
        
        mixers(:,f,:,:) = tmpMixers;
        if (~isempty(tmpDiodes))
            diodes(:,f,:,:) = tmpDiodes; %#ok<AGROW>
        end
    end
    

end

