classdef subscribeAndPlot < handle
    %subscribeAndPlot 
    %   Small class used to suscribe some signals and plot them online.
    %   Just create an object as
    %     myobj = subscribeAndPlot(<selector>, <signals cell array>);
    %   and start acquiring data
    %     myobj.startAcquisition();
    
    properties
        % cell list with labels for different signals. It is automatically
        % populated, but you can change it later.
        signalsLabels;
        % colors for the different plots. Automatically generated, but you
        % can change it later.
        signalsColors;
        
        % cell with all last data
        lastData;
        lastProcessedPhaseData;
        lastDataTreatedY;
        lastDataTreatedX;
        
        % to be specified to offsets the different signals
        xoffsets = []; 
        % to be specified to offsets the different signals
        yoffsets = []; 
        % to be specified to scale the different signals
        xscalingfactors = [];
        % to be specified to scale the different signals
        yscalingfactors = [];
        
        % ylim plot
        yLims = [];
        % xlim plot
        xLims = [];
        % if to plot in a persistent way
        persistentPlot = false;
        
        % stop acquisition if there are NaNs in the data.
        stopIfNaNs = false;
        
        %matlabJacMonitor object
        datamonitor;
        
        % figure handler (don't touch)
        myFigureHandle;
    end
    properties (Access=private)
        % cell list of signals in use (DON'T CHANGE IT!).
        signals;
    end
    methods
        % constructor
        function obj=subscribeAndPlot(cycleName, signals)
            obj.signals = signals;
            for i=1:numel(signals)
                obj.signalsLabels{i} = matlabDataAndSignalsHelper.decomposeSignal(signals{i});
                obj.signalsColors{i} = rand(1,3);
            end

            obj.datamonitor = matlabJapcMonitor(cycleName, signals, ...
                @(e)updateDispersionFunction(obj,e),'');
        end
        
        function stopAcquisition(obj)
            obj.datamonitor.stop();
            obj.myFigureHandle = [];
        end
        function startAcquisition(obj)
            obj.datamonitor.start(1);
        end
        function pauseAcquisition(obj)
            obj.datamonitor.pauseOn();
        end
        function resumeAcquisition(obj)
            obj.datamonitor.pauseOff();
        end
        function delete(obj)
            obj.stopAcquisition();
            obj.datamonitor.delete();
        end

        function plotSignals(obj)
            if isempty(obj.myFigureHandle)
                obj.myFigureHandle = figure();
            end
            set(obj.myFigureHandle,'DeleteFcn',@(h,e)obj.stopAcquisition);
            
            tmpCurrentFigure = get(0,'CurrentFigure');
            set(0,'CurrentFigure',obj.myFigureHandle);            

            % plot
            if iscell(obj.lastDataTreatedY)
                for i=1:numel(obj.lastDataTreatedY)
                    plot(obj.lastDataTreatedX{i},obj.lastDataTreatedY{i},...
                        '.-','Color',obj.signalsColors{i})
                    hold on
                end
            else
                for i=1:size(obj.lastDataTreatedY,1)
                    plot(obj.lastDataTreatedX(i,:),obj.lastDataTreatedY(i,:),...
                        '.-','Color',obj.signalsColors{i})
                    hold on
                end
            end
            if obj.persistentPlot
                hold on
            else
                hold off
            end
            legend(obj.signalsLabels);
            if ~isempty(obj.xLims)
                xlim(obj.xLims);
            end
            if ~isempty(obj.yLims)
                ylim(obj.yLims);
            end
            grid

            % reset back figure was in use before
            set(0,'CurrentFigure',tmpCurrentFigure);
        end
    end
    methods(Access=private)
        function updateDispersionFunction(obj, dataStruct)
            %acquire data
            auxCell = cell(0);
            
            for i=1:length(obj.signals)
                auxCell{i} = matlabDataAndSignalsHelper.simpleExtractSingleSignal(dataStruct, obj.signals{i});
            end

            if obj.stopIfNaNs
                for i=1:numel(auxCell)
                    if hasInfNaN(auxCell{i})
                        %disp('NaN or Inf in current data. Skip it.')
                        return;
                    end
                end
            end

            obj.lastData = auxCell;
            
            % call data treatment function
            obj.dataTreat();
            
            % call plot function
            obj.plotSignals();
        end
        function dataTreat(obj)
            % if you want to do something on acquired data, here should go
            % your code.
            % The new data is available in obj.lastData

            % apply offsets and scaling factors
            % check
            if numel(obj.xscalingfactors) > 0 && ...
                    numel(obj.signals) ~= numel(obj.xscalingfactors)
                error('Wrong number of x scaling factors.')
            elseif numel(obj.xscalingfactors) == 0
                obj.xscalingfactors = ones(numel(obj.signals),1);
            end
            if numel(obj.yscalingfactors) > 0 && ...
                    numel(obj.signals) ~= numel(obj.yscalingfactors)
                error('Wrong number of y scaling factors.')
            elseif numel(obj.yscalingfactors) == 0
                obj.yscalingfactors = ones(numel(obj.signals),1);
            end
            if numel(obj.xoffsets) > 0 && ...
                    numel(obj.signals) ~= numel(obj.xoffsets)
                error('Wrong number of x offsets.')
            elseif numel(obj.xoffsets) == 0
                obj.xoffsets = zeros(numel(obj.signals),1);
            end
            if numel(obj.yoffsets) > 0 && ...
                    numel(obj.signals) ~= numel(obj.yoffsets)
                error('Wrong number of x offsets.')
            elseif numel(obj.yoffsets) == 0
                obj.yoffsets = zeros(numel(obj.signals),1);
            end
            
            % apply to data
            obj.lastDataTreatedY = obj.lastData;
            obj.lastDataTreatedX = [];
            if iscell(obj.lastDataTreatedY)
                for i=1:numel(obj.lastDataTreatedY)
                    obj.lastDataTreatedY{i} = obj.yoffsets(i)+...
                        obj.yscalingfactors(i)*obj.lastDataTreatedY{i};
                    obj.lastDataTreatedX{i} = obj.xoffsets(i)+...
                        obj.xscalingfactors(i)*(1:numel(obj.lastDataTreatedY{i}));
                end
            else
                for i=1:size(obj.lastDataTreatedY,1)
                    obj.lastDataTreatedY(i,:) = obj.yoffsets(i)+...
                        obj.yscalingfactors(i)*obj.lastDataTreatedY(i,:);
                    obj.lastDataTreatedX(i,:) = obj.xoffsets(i)+...
                        obj.xscalingfactors(i)*(1:numel(obj.lastDataTreatedY(i,:)));
                end
            end
        end
    end
    
end

