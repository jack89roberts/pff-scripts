dataSetName= '20150930_1846_Gun135';
dataDir = '/home/jack/PhaseFeedforward/Analysis/201509';

dataSetName= '20150715_1414_FF_Gain63_Gate185_350_Int_Odd';
dataDir = '/home/jack/PhaseFeedforward/Analysis/201507';

%plotXRange = [535 702];

%% load data
addpath('../');

addpath('/home/jack/Documents/MATLAB/madx2Matlab');

load([dataDir '/' dataSetName '.mat']);

%%
figure;
for i=1:nBPMs
    if (useBPM(i))
        if (~isempty(bpmPulseRange{i}))
            ax1 = subplot(1,3,1);
            plot(meanBPMHAlongPulse{i});
            title(bpmNames{i})
            xlabel('Sample No.');
            ylabel('Horizontal Position [mm]');
            grid
            hold all;
            yLim = get(gca,'YLim');
            plot([bpmPulseRange{i}(1) bpmPulseRange{i}(1)],yLim,'k');
            plot([bpmPulseRange{i}(end) bpmPulseRange{i}(end)],yLim,'k');
            hold off;

            ax2 = subplot(1,3,2);
            plot(meanBPMVAlongPulse{i});
            title(bpmNames{i})
            xlabel('Sample No.');
            ylabel('Vertical Position [mm]');
            grid
            hold all;
            yLim = get(gca,'YLim');
            plot([bpmPulseRange{i}(1) bpmPulseRange{i}(1)],yLim,'k');
            plot([bpmPulseRange{i}(end) bpmPulseRange{i}(end)],yLim,'k');
            hold off;

            ax3 = subplot(1,3,3);
            plot(meanBPMSAlongPulse{i});
            title(bpmNames{i})
            xlabel('Sample No.');
            ylabel('Current [A]');
            grid
            yLim = get(gca,'YLim');
            hold all;
            plot([bpmPulseRange{i}(1) bpmPulseRange{i}(1)],yLim,'k');
            plot([bpmPulseRange{i}(end) bpmPulseRange{i}(end)],yLim,'k');
            hold off;

            linkaxes([ax3,ax2,ax1],'x');
%             xlim([bpmPulseRange{i}(1)+10 bpmPulseRange{i}(end)-5]);
            
            input(sprintf('BPM %d of %d',i,nBPMs));
        end
    end
end