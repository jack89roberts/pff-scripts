function [ pulseSampleRange, staticSampleDelay ] = pickPulseRange( mixers,diodes,phases )
%pickPulseRange Plot and prompt user to pick sample range of signal pulse

    if (isempty(mixers))
        useMixer = 0;
    else
        useMixer = 1;
    end
    if (nargin<2 || isempty(diodes))
        useDiode = 0;
    else
        useDiode = 1;
    end
    if (nargin<3 || isempty(phases))
        usePhase = 0;
    else
        usePhase = 1;
    end
    
    [nMons,~,~] = size(mixers);
    
    % plot signals
    nSubPlots = useMixer+useDiode+usePhase;
    plotSize = 10*(10+nSubPlots);
    for mon=1:nMons
        figs(mon) = figure;
        subPlotInd = 1;
        plotIndex = plotSize+subPlotInd;
        if (useMixer)
            axes(plotIndex) = subplot(plotIndex);
            plot(squeeze(mixers(mon,:,:))');
            xlabel('Sample No.')
            ylabel('Output [V]')
            title(sprintf('MON %d MIXER',mon))
            subPlotInd = subPlotInd+1;
            plotIndex = plotSize+subPlotInd;
        end
        if (useDiode)
            axes(plotIndex) = subplot(plotIndex);
            plot(squeeze(diodes(mon,:,:))');
            xlabel('Sample No.')
            ylabel('Output [V]')
            title(sprintf('MON %d Diode',mon))
            subPlotInd = subPlotInd+1;
            plotIndex = plotSize+subPlotInd;
        end
        if (usePhase)
            axes(plotIndex) = subplot(plotIndex);
            plot(squeeze(phases(mon,:,:))');
            xlabel('Sample No.')
            ylabel('Phase [degrees]')
            title(sprintf('MON %d Phase',mon))
        end
        linkaxes(axes,'x');
    end
    
    % prompt for sample range
    pulseSampleRange = cell(1,nMons);
    staticSampleDelay = zeros(1,nMons);
    for mon=1:nMons
        fprintf('MONITOR %d\n',mon)
        pulseSampleRange{mon} = input('Sample Range: ');
        fprintf('---------------------------------\n');
        
        staticSampleDelay(mon) = pulseSampleRange{1}(end)-pulseSampleRange{mon}(end);
    end
    
    % add sample range to plots
    for mon=1:nMons
        figure(figs(mon));
        subPlotInd = 1;
        plotIndex = plotSize+subPlotInd;
        if (useMixer)
            subplot(plotIndex);
            hold all;
            plot([pulseSampleRange{mon}(1) pulseSampleRange{mon}(1)],get(gca,'YLim'),'k');
            plot([pulseSampleRange{mon}(end) pulseSampleRange{mon}(end)],get(gca,'YLim'),'k');
            subPlotInd = subPlotInd+1;
            plotIndex = plotSize+subPlotInd;
        end
        if (useDiode)
            subplot(plotIndex);
            hold all;
            plot([pulseSampleRange{mon}(1) pulseSampleRange{mon}(1)],get(gca,'YLim'),'k');
            plot([pulseSampleRange{mon}(end) pulseSampleRange{mon}(end)],get(gca,'YLim'),'k');
            subPlotInd = subPlotInd+1;
            plotIndex = plotSize+subPlotInd;
        end
        if (usePhase)
            subplot(plotIndex);
            hold all;
            plot([pulseSampleRange{mon}(1) pulseSampleRange{mon}(1)],get(gca,'YLim'),'k');
            plot([pulseSampleRange{mon}(end) pulseSampleRange{mon}(end)],get(gca,'YLim'),'k');
            subPlotInd = subPlotInd+1;
            plotIndex = plotSize+subPlotInd;
        end
    end

end

