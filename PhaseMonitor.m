classdef PhaseMonitor
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        saveData = false;
        saveDataPath = 'data/';
        dataSetComment = 'JackComment'; % comment for this data set to be save in all data structs
    end
    
    properties (Access=private)
        cycleName = 'SCT.USER.ALL';
        monitor;
        signalList={};
        monitorIdxs=[];
    end
    
    
    methods
        % constructor
        function obj = PhaseMonitor()
            addpath('/user/ctf3op/davideg/linearFeedback/lib/');
            addpath('/nfs/cs-ccr-nfs6/vol29/share/ctf/data/ctfmod/MatLab/common');
            
            obj.monitor = matlabMonitor(cycleName,signalList,@phaseDataReceived,dataSetComment,monitorIdxs);
            obj.monitor.saveData = saveData;
            obj.monitor.saveDataPath = saveDataPath;

        end
        
        function loadDevices()
            [obj.signalList, obj.monitorIdxs] = combinePropertyLists(... 
                                makePropertyListFromFiles('devices/bpm.devs','devices/bpm.props'),...
                                makePropertyListFromFiles('devices/bpmCT.devs','devices/bpm.props'),...
                                makePropertyListFromFiles('devices/bpr.devs','devices/bpr.props'),...
                                makePropertyListFromFiles('devices/scope.devs','devices/scope.props'),...
                                makePropertyListFromFiles('devices/magnetCC.devs','devices/magnet.props'),...
                                makePropertyListFromFiles('devices/magnetCB.devs','devices/magnet.props'),...
                                makePropertyListFromFiles('devices/corrsCC.devs','devices/magnet.props')...                                 
                            );
        end
        
        function startMonitoring()
            start(obj.monitor);
        end
        
        function stopMonitoring()
            stop(obj.monitor);
        end
        
        function saveData(nCycles)
            if (nCycles > 0) 

                recCyc = 0;
                while(phaseMonitor.recordedCycles < nCycles)
                    if (phaseMonitor.recordedCycles ~= recCyc)
                        recCyc = phaseMonitor.recordedCycles;
                        fprintf('Trigger %d out of %d received.\n',recCyc,nCycles);
                    end
                    pause(0.5);
                end
                fprintf('Trigger %d out of %d received.\n',phaseMonitor.recordedCycles,nCycles);
            
            end

    end
    
end

