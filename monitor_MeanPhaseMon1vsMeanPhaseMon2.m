close all;
clear all;
addpath('/user/ctf3op/davideg/linearFeedback/lib/');
addpath('../common');

%% INPUTS

cycleName = 'SCT.USER.ALL';
saveData = true;
saveDataPath = 'data/20141121_1830_JustPhaseMonitors/'; % MUST HAVE / AT END
nCycles = 100; % no. pulses to save, if nCycles<=0 run continuously
dataSetComment = ''; % comment for this data set to be save in all data structs

%inputDevicesFile = 'devices.txt';

%% Load device property lists and create monitor object

fprintf('Loading device lists...\n');

% make data directory if needed
if (saveData == true)

    if( ~exist(saveDataPath,'dir') ) 
        mkdir(saveDataPath);
    end;
    
    if(nCycles <= 0)
        fprintf(1,'WARNING: Data saving enabled in continuous mode! All pulses will be saved whilst subscription is on.');
    end;
end

[signalList, monitorIdxs] = combinePropertyLists(... 
                                ...%makePropertyListFromFiles('devices/bpmCC.devs','devices/bpm.props'),...
                                ...%makePropertyListFromFiles('devices/bpmCB.devs','devices/bpm.props'),...
                                ...%makePropertyListFromFiles('devices/bpmCT.devs','devices/bpm.props'),...
                                ...%makePropertyListFromFiles('devices/bpmCM.devs','devices/bpm.props'),...
                                ...%makePropertyListFromFiles('devices/bpr.devs','devices/bpr.props'),...
                                ...%makePropertyListFromFiles('devices/scope.devs','devices/scope.props')...
                                ...%makePropertyListFromFiles('devices/corrs.devs','devices/magnet.props')... 
                            );


% importfile(inputDevicesFile);
% 
% signalList={};
% monitorIdxs=[];
% 
% for i=1:numel(devices)
%     if(~isempty(cell2mat(strfind(devices(i),'CE.'))) || ~isempty(cell2mat(strfind(devices(i),'CT.SCOPE'))) )
%         % properties for OASIS scopes
%         prop1 = [devices{i} '/Acquisition#value'];
%         prop2 = [devices{i} '/Sensibility#value'];
%         prop3 = [devices{i} '/Offset#value'];
%         
%         signalList = [signalList prop1 prop2 prop3];
%         
%         % to create one monitor for OASIS scopes and another for Xeneric
%         monitorIdxs = [monitorIdxs 2 2 2];
%     
%     elseif (~isempty(cell2mat(strfind(devices(i),'BPM'))) || ~isempty(cell2mat(strfind(devices(i),'BPI'))) || ~isempty(cell2mat(strfind(devices(i),'BPS'))))
%         prop1 = [devices{i} 'H/MeanAtCursor#mean'];
%         prop2 = [devices{i} 'S/MeanAtCursor#mean'];
%         
%         signalList = [signalList prop1 prop2];
%         
%         monitorIdxs = [monitorIdxs 3 3];
%         
%     else
%         % properties for Xeneric samplers
%         prop1 = [devices{i} '/Samples#samples'];
%         
%         signalList = [signalList prop1];
%         
%          % to creatae one monitor for OASIS scopes and another for Xeneric
%         monitorIdxs = [monitorIdxs 1];
%     end
% end

fprintf('Creating monitor object...\n');

%global myFactor;
%myFactor = 4;
%phaseMonitor = matlabMonitor(cycleName,signalList,@(data)phaseDataReceived(data),dataSetComment,monitorIdxs);
phaseMonitor = matlabMonitor(cycleName,signalList,@phaseDataReceived,dataSetComment,monitorIdxs);
phaseMonitor.saveData = saveData;
phaseMonitor.saveDataPath = saveDataPath;

%%
if (saveData == true)

    fprintf('Saving initial device settings...\n');

    dateNow = datevec(now);
    initValsFileName = sprintf('%sinitSettings_%02d%02d%02d_%02d%02d_%.0f.dat', saveDataPath, dateNow(1), dateNow(2), dateNow(3) ,dateNow(4), dateNow(5) ,dateNow(6));
    savePropValues('devices/init.devs','devices/magnet.props', initValsFileName);

end

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