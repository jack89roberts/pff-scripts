function [FONTData] = saveExtractedFONTData( dataSetName,dataDir )
%saveExtractedFONTData Extracts FONT data from dataSetName in to a struct and
%then saves it to the directory 'Extracted' in dataDir.

    extractedDir = [dataDir '/Extracted'];
    saveName = [extractedDir '/' dataSetName '.mat'];

    fprintf('Saving data to %s...\n',extractedDir);
    
    if (~exist(extractedDir,'dir'))
        mkdir(extractedDir);
    end
    
    
%% Everything below here based on Doug's script.
% Look for "JACK" comments for where I changed things.
    NChannels = 9;       
    %% get data from .dat file
    fid = fopen([dataDir '/' dataSetName '.dat']); % JACK - Changed for Linux paths.

    % extract data from the file header line
    FileHeader = fgetl(fid);
    LocSpace = strfind(FileHeader,' ');
    RunDate = FileHeader(LocSpace(3)+1:LocSpace(4)-1);
    RunTime = FileHeader(LocSpace(4)+1:end);

    ADCs = [];

    % treat the first trigger separately to account for file header 
    TriggerCtr = 1;
    disp([dataSetName ': Processing trigger number ' num2str(TriggerCtr) '...'])              
    BoardStamp(TriggerCtr) = sscanf(fgetl(fid),'%*s %*s %d');
    TimestampLine = fgetl(fid);
    TriggerNo(TriggerCtr) = sscanf(TimestampLine,'%d',1);
    Time(TriggerCtr,:) = sscanf(TimestampLine,'%*d\t%d:%d:%f');
    fgetl(fid);
    ADCs = [ ADCs; fscanf(fid,'%d') ];

    SamplesPerTrigger = length(ADCs);
    SamplesPerChannel = SamplesPerTrigger / 9;

    % process the second trigger and determine how many header lines were in the file
    TriggerCtr = 2;
    disp([dataSetName ': Processing trigger number ' num2str(TriggerCtr) '...'])              

    TestLine = fgetl(fid);
    if strcmp(TestLine(1:3),'Run')
        % this is the file header line and is being appended to every pulse
        NHeaderLines = 4;
        BoardStampLine = fgetl(fid);
    else
        % this line contains the board stamp
        NHeaderLines = 3;
        BoardStampLine = TestLine;
    end

    BoardStamp(TriggerCtr) = sscanf(BoardStampLine,'%*s %*s %d');
    TimestampLine = fgetl(fid);
    TriggerNo(TriggerCtr) = sscanf(TimestampLine,'%d',1);
    Time(TriggerCtr,:) = sscanf(TimestampLine,'%*d\t%d:%d:%f');
    fgetl(fid);
    ADCs = [ ADCs; fscanf(fid,'%d') ];

     % process subsequent triggers
    while ~feof(fid) %&& TriggerCtr < 500 % JACK - not sure why the limit of 500? Maybe DAQ saves to a new file after 500 triggers?
        TriggerCtr = TriggerCtr + 1;
        disp([dataSetName ': Processing trigger number ' num2str(TriggerCtr) '...'])            

        % skip the file header line if it was appended to every trigger
        for lineNo = 1 : (NHeaderLines - 3)
            fgetl(fid);
        end

        BoardStamp(TriggerCtr) = sscanf(fgetl(fid),'%*s %*s %d');
        TimestampLine = fgetl(fid);
        TriggerNo(TriggerCtr) = sscanf(TimestampLine,'%d',1);
        Time(TriggerCtr,:) = sscanf(TimestampLine,'%*d\t%d:%d:%f');

        % skip column headers for ADC data
        fgetl(fid);

        CurrentTriggerData = fscanf(fid,'%d');
        % catch glitched triggers (rare)
        if length(CurrentTriggerData)~=SamplesPerTrigger
            disp('WARNING! Wrong number of ADC samples detected this extraction. Current data replaced by NaN.')
            CurrentTriggerData = nan(SamplesPerTrigger,1);
        end
        ADCs = [ ADCs; CurrentTriggerData ];
    end
    % read in data from the remaining extractions

    % all data extracted from .dat file
    NTotalSamples = length(ADCs);
    ADCs = reshape(ADCs,NChannels,SamplesPerChannel,NTotalSamples/(NChannels*SamplesPerChannel));
    ADCs = permute(ADCs,[1 3 2]); % JACK - I use this order of dimensions: (ADCindex, triggerIndex, sampleIndex)
    fclose(fid);

    %% get data from .rb file

% ORIGINAL CODE FROM DOUG
%     rbPathname = [dataDir '/' dataSetName '.rb']; % JACK - Changed for Linux paths.
%     CheckFileExists = exist(rbPathname,'file');
% 
%     if ~CheckFileExists
%         disp('No readback data found. Readback matrix will be empty.')
%         rbs = [];
%     else
%         fid = fopen(rbPathname);
% 
%         rbs = fscanf(fid,'%d');
% 
%         fclose(fid);
% 
%         NTriggers = TriggerCtr;
% 
%         NReadbacks = length(rbs)/NTriggers;
% 
%         if ~rem(NReadbacks,1)
%             if NReadbacks*NTriggers == length(rbs)
%                 rbs = reshape(rbs,NReadbacks,NTriggers)';
%             else
%                 rbs = [];
%             end
%         else
%             rbs = [];
%             disp('WARNING! Total number of readbacks not as expected.')
%         end
%     end

% JACK CODE - With basic error checking for bad triggers
    NTriggers = TriggerCtr;    
    NRBs = 144;

    rbPathname = [dataDir '/' dataSetName '.rb']; % JACK - Changed for Linux paths.
    CheckFileExists = exist(rbPathname,'file');

    if ~CheckFileExists
        disp('No readback data found. Readback matrix will be empty.')
        rbs = [];
    else        
        fid = fopen(rbPathname);

        TriggerCtr = 0;
        rbs = NaN(NTriggers,NRBs);

        while ~feof(fid)
            TriggerCtr = TriggerCtr + 1;
            disp([dataSetName ': Processing readback for trigger ' num2str(TriggerCtr) '...'])            

            CurrentTriggerRB=sscanf(fgetl(fid),'%d');

            if (length(CurrentTriggerRB) ~= NRBs)
                disp('WARNING! Wrong number of readbacks detected this extaction. Current data replaced by NaN');
            else
                rbs(TriggerCtr,:) = CurrentTriggerRB;
            end
        end

        fclose(fid);
    end

    %% get data from .dac file
    DACPathname = [dataDir '/' dataSetName '.dac']; % JACK - Changed for Linux paths
    CheckFileExists = exist(DACPathname,'file');

    if ~CheckFileExists
        disp('No DAC data found. DAC matrix will be empty.')
        onDACs = [];
    else
        fid = fopen(DACPathname);

        onDACs = fscanf(fid,'%d');

        fclose(fid);

        %NDACs = length(DACs)/TriggerCtr;
        NDACSamp = (SamplesPerChannel/2)+1; % trailing zero
        DACTrigCtr = length(onDACs)/NDACSamp;
        
        if (DACTrigCtr ~= TriggerCtr)
            fprintf('WARNING! DACTrigCtr=%d not equal TriggerCtr=%d. Unsure if this is stable.\n',DACTrigCtr,TriggerCtr);
        end
        
%        if rem(NDACSamp,1) == 0 % JACK - this doesn't do anything?
          if (~isempty(onDACs))
            onDACs = reshape(onDACs,NDACSamp,DACTrigCtr)';
            
            % JACK - split DAC1 and DAC2
            DAC1 = onDACs(1:2:end,1:(end-1));
            DAC2 = onDACs(2:2:end,1:(end-1));
            
            NDACSamp = NDACSamp-1; % to remove trailing zero
            DACTrigCtr = DACTrigCtr/2; %
            
            onDACs = NaN(2,DACTrigCtr,NDACSamp);
            onDACs(1,:,:) = DAC1;
            onDACs(2,:,:) = DAC2;
            %DACs = permute(DACs,[1 3 2]); % JACK - I use this order of dimensions: (DACindex, triggerIndex, sampleIndex)
%               DACs = [];
        else
            onDACs = [];
            DACs = [];
            disp('WARNING! Total number of DAC values not as expected.')
        end
    end
    
    %% Find FF status, match DAC to ADC trigger, zero array when off
    if (~isempty(rbs))
        FFStatus = mod(rbs(:,143),2);

        if (~isempty(onDACs))

            DACs = zeros(2,NTriggers,NDACSamp);
            DACInd = 0;

            for i=1:NTriggers
                
                if (FFStatus(i)==1)
                    DACInd = DACInd+1;
                    DACs(:,i,:) = onDACs(:,DACInd,:);
                end
            end

        end
    else
        FFStatus = [];
    end

    %% Calculate channel offset
    % uses RBs 107 for 6 most sig bits, and 106 for 7 least sig bits (ind
    % 109 and 108 in MatLab)

    if (~isempty(rbs))
        ChanOffset = (rbs(:,109).*(2^7)) + rbs(:,108);
    else
        ChanOffset = [];
    end

    %% assemble data into a structure and save
    % JACK - changes in saved names for structure etc.
    FONTData = struct('RunName',dataSetName,'RunDate',RunDate,'RunTime',RunTime,'ADCs',ADCs,'BoardStamp',BoardStamp,'Time',Time,'rbs',rbs,'DACs',DACs,'FFStatus',FFStatus,'ChanOffset',ChanOffset);
    save(saveName,'FONTData');
    disp(['Extraction complete: files saved to ' saveName ])

end

