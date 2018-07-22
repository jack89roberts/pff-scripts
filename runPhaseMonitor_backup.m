%% INPUTS
close all;
clear all;

cycleName = 'SCT.USER.ALL';
saveData = false;
saveDataPath = 'data/';
nCycles = -5; % no. pulses to save, if nCycles<=0 run continuously
dataSetComment = 'JackComment'; % comment for this data set to be save in all data structs

inputDevicesFile = 'devices.txt';


%%
addpath('/user/ctf3op/davideg/linearFeedback/lib/');
addpath('../common');

fprintf('Loading device list...\n');

if (saveData == true)

    if( ~exist(saveDataPath,'dir') ) 
        mkdir(saveDataPath);
    end;
    
    if(nCycles <= 0)
        fprintf(1,'WARNING: Data saving enabled in continuous mode! All pulses will be saved whilst subscription is on.');
    end;
end

importfile(inputDevicesFile);

signalList={};
monitorIdxs=[];

for i=1:numel(devices)
    if(~isempty(cell2mat(strfind(devices(i),'CE.'))) || ~isempty(cell2mat(strfind(devices(i),'CT.SCOPE'))) )
        % properties for OASIS scopes
        prop1 = [devices{i} '/Acquisition#value'];
        prop2 = [devices{i} '/Sensibility#value'];
        prop3 = [devices{i} '/Offset#value'];
        
        signalList = [signalList prop1 prop2 prop3];
        
        % to create one monitor for OASIS scopes and another for Xeneric
        monitorIdxs = [monitorIdxs 2 2 2];
        
    else
        % properties for Xeneric samplers
        prop1 = [devices{i} '/Samples#samples'];
        
        signalList = [signalList prop1];
        
         % to creatae one monitor for OASIS scopes and another for Xeneric
        monitorIdxs = [monitorIdxs 1];
    end
end

fprintf('Creating monitor object...\n');

global myFactor;
myFactor = 4;
phaseMonitor = matlabMonitor(cycleName,signalList,@(data)phaseDataReceived(data),dataSetComment,monitorIdxs);
phaseMonitor.saveData = saveData;
phaseMonitor.saveDataPath = saveDataPath;

%%

fprintf('Starting subscription...\n');

start(phaseMonitor);

if (nCycles > 0) % user has set a fixed number of pulses to save/view

    if (phaseMonitor.recordedCycles > 0)
        nMissedStart = phaseMonitor.recordedCycles;
        fprintf('Starting subscription took some time. Adding %d cycles to set limit.\n',nMissedStart);
        nCycles = nCycles + nMissedStart;
    end

    recCyc = 0;
    while(phaseMonitor.recordedCycles < nCycles)
        if (phaseMonitor.recordedCycles ~= recCyc)
            recCyc = phaseMonitor.recordedCycles;
            fprintf('Trigger %d out of %d received.\n',recCyc,nCycles);
        end
        pause(0.5);
    end
    fprintf('Trigger %d out of %d received.\n',phaseMonitor.recordedCycles,nCycles);

    fprintf('Stopping subscription...\n');
    stop(phaseMonitor);

else
    fprintf('Subscription started in continuous mode.\nUse stop(phaseMonitor) to end subscription.\n');
    return;
end;