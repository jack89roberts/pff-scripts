function bpmSampleRange = getBPMSampleRange( bpmS,bpmH )
%bpmSampleRange Prompt for bpm sample range to use
%   bpmS and bpmH should contain only one pulse, arrays in format
%   bpmIndex,sampleNo.
    
    if (nargin<2)
        bpmNames = {...
            'CC.SVBPM0435',...
            'CC.SVBPI0535',...
            'CC.SVBPI0645',...
            'CC.SVBPI0685',...
            'CC.SVBPI0735',...
            'CC.SVBPM0845',...
            'CC.SVBPM0930',...
            'CB.SVBPM0150',...
            'CB.SVBPS0210',...
            'CB.SVBPS0250',...
            'CB.SVBPS0310'...
        };
        nBPMs = length(bpmNames);
        nSamples = [];
        
        for i=1:nBPMs
            tmpBPMS = matlabJapc.staticGetSignal('SCT.USER.SETUP',[bpmNames{i} 'S/Samples#samples']);
            tmpBPMH = matlabJapc.staticGetSignal('SCT.USER.SETUP',[bpmNames{i} 'H/Samples#samples']);
            if (isempty(nSamples))
                nSamples = length(tmpBPMS);
            end
            
            bpmS(i,:) = tmpBPMS;
            bpmH(i,:) = tmpBPMH;
        end

    else
        [nBPMs,~] = size(bpmS);
    end
    
    figure;
    
    for i=1:nBPMs
        [tmpBPMS,tmpBPMH] = getAlignedXCorr(bpmS(i,:),'end',{bpmH(i,:)});
        tmpBPMH = tmpBPMH{1};
        subplot(1,2,1)
        plot(tmpBPMS);
        hold all;
        title('BPMS');
        xlabel('Sample No,');
        ylabel('Current [A]');
        subplot(1,2,2);
        plot(tmpBPMH);
        hold all;
        title('BPMH');
        xlabel('Sample No,');
        ylabel('Position [mm]');
    end

    startSamp = input('Start sample: ');
    endSamp = input('End sample: ');
    bpmSampleRange = startSamp:endSamp;
    
    subplot(1,2,1);
    yLim = get(gca,'YLim');
    plot([startSamp startSamp],yLim,'k');
    plot([endSamp endSamp],yLim,'k');
    
    subplot(1,2,2);
    yLim = get(gca,'YLim');
    plot([startSamp startSamp],yLim,'k');
    plot([endSamp endSamp],yLim,'k');

    
end

