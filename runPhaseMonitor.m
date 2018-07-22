close all;
clear all;
%% INPUTS - always need to be changed.

dataSetName = 'Waveform2Change_Gain_0'; % Name of the data set (without date/time - now added automatically)
nCycles = 10; % no. pulses to save, if nCycles<=0 run continuously
dataSetComment = ''; % description of the data set.

saveData = false;
saveInit = true;

%% INPUTS - shouldn't need to be changed.
cycleName = 'SCT.USER.SETUP';
saveDataLocalDir = '/local/home/ctf3op/data/PhaseFeedforward/'; % must have / at end
%saveDataRemoteDir = '/nfs/cs-ccr-nfs6/vol29/share/ctf/data/ctfmod/MatLab/PhaseFeedforward/data/'; % must have / at end
saveDataRemoteDir = '/user/ctf3op/PhaseFeedforward/Data/'; % must have / at end
warning('off','backtrace') % simplifies warning messages (may help to turn back on if debugging needed)
%addpath(genpath('/ctf/data/ctfmod/MatLab/JavaCoInterface2/'));
%addpath('../common');
cd('/nfs/cs-ccr-nfs6/vol29/share/ctf/data/ctfmod/MatLab/PhaseFeedforward');

%% Load device property lists and create monitor object
fprintf('Loading device lists...\n');

[signalList, monitorIdxs] = combinePropertyLists(... 
                                makePropertyListFromFiles('devices/scope.devs','devices/scope.props'),...
                                makePropertyListFromFiles('devices/bpmCC.devs','devices/bpm.props'),...
                                makePropertyListFromFiles('devices/bpmCB.devs','devices/bpm.props'),...
                                makePropertyListFromFiles('devices/bpmCT.devs','devices/bpm.props'),...
                                makePropertyListFromFiles('devices/bpr.devs','devices/bpr.props'),...
                                makePropertyListFromFiles('devices/corrs.devs','devices/magnet.props'),...                                
                                makePropertyListFromFiles('devices/pki.devs','devices/klystron.props'),... 
                                makePropertyListFromFiles('devices/psi.devs','devices/klystron.props'),... 
                                makePropertyListFromFiles('devices/bpmCL.devs','devices/bpm.props'),... 
                                makePropertyListFromFiles('devices/bpmCR.devs','devices/bpm.props'),...
                                makePropertyListFromFiles('devices/oasisScope.devs','devices/oasisScope.props'),...
...%                            makePropertyListFromFiles('devices/frascatiBend.devs','devices/magnet.props'),...
                                makePropertyListFromFiles('devices/gun.devs','devices/gun.props')...
...%                                 ...%makePropertyListFromFiles('devices/septum.devs','devices/septum.props')...
...%                                 makePropertyListFromFiles('devices/bpmAnalogue.devs','devices/bpm.props'),...
...%                                 makePropertyListFromFiles('devices/bpiAnalogue.devs','devices/bpiAnalogue.props')...
...%                                  makePropertyListFromFiles('devices/corrsCT.devs','devices/corrCT.props')...                                
                            );

fprintf('Creating monitor object...\n');

%phaseMonitor = matlabJapcMonitor(cycleName,signalList,@phaseDataReceived,dataSetComment,monitorIdxs);
phaseMonitor = matlabJapcMonitor(cycleName,signalList,@phaseDataReceived,dataSetComment);
phaseMonitor.saveData = saveData;
phaseMonitor.useFastStrategy(1);

%% If saving data, create full paths, make directories, save initial settings

% create full data paths etc. for data saving (save locally first)
if (saveData == true)
    fprintf('Making data directory...\n');
    
    dV = datevec(now);
    dataTimeStamp = sprintf('%02d%02d%02d_%02d%02d',dV(1),dV(2),dV(3),dV(4),dV(5));
    dataSetLocalPath = [saveDataLocalDir dataTimeStamp '_' dataSetName];
    
    % make the local data directory
    % if the directory already exists increment the time stamp in steps
    % of 1 minute until a unique name is found.
    if( ~exist(dataSetLocalPath,'dir') ) 
        mkdir(dataSetLocalPath);
    else
        while(exist(dataSetLocalPath,'dir'))
            dV(5) = dV(5) + 1;
            if (dV(5) == 60)
                dV(4) = dV(4) + 1;
                dV(5) = 0;
                dataTimeStamp = sprintf('%02d%02d%02d_%02d%02d',dV(1),dV(2),dV(3),dV(4),dV(5));
            end
            warning('Data set %s already exists. Changed time stamp to %s',dataSetLocalPath,dataTimeStamp);
            dataSetLocalPath = [saveDataLocalDir dataTimeStamp '_' dataSetName];
        end
        mkdir(dataSetLocalPath);
    end;
    
    phaseMonitor.saveDataPath = dataSetLocalPath;

    fprintf('Created local path: %s\n',dataSetLocalPath);
 
    if(nCycles <= 0)
        warning('Data saving enabled in continuous mode! All pulses will be saved whilst subscription is on.');
    end;

    if (saveInit)
        fprintf('Saving initial device settings...\n');

        initValsFileName = sprintf('%s/initSettings.dat', dataSetLocalPath);
        savePropValues('devices/init.devs','devices/magnet.props', initValsFileName);

        % Time stamp of current reference Frascati calibration
        refFile = fopen('FrascatiCalibrations/frascatiRefCalName');
        tmp = textscan(refFile,'%s');
        fclose(refFile);
        stamp = tmp{1}{1};

        initValsFile = fopen(initValsFileName,'a');
        fprintf(initValsFile,'FrascatiRefCal, %s',stamp);
        fclose(initValsFile);
    end
end

%%

fprintf('Starting subscription...\n');

phaseMonitor.start();

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
    phaseMonitor.stop();
    phaseMonitor.delete();
    
    if (saveData)
        dataSetRemotePath = [saveDataRemoteDir dataTimeStamp '_' dataSetName];        
        fprintf('Copying data to %s...\n',dataSetRemotePath);        
        copyfile(dataSetLocalPath,dataSetRemotePath)
        fprintf('Data Set Name: %s_%s\n',dataTimeStamp,dataSetName);
    end

else
    fprintf('Subscription started in continuous mode.\nUse stop(phaseMonitor) to end subscription.\nData will only be saved locally.\n');
end;

