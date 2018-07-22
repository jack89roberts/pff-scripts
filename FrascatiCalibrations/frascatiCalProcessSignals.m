function [alignedMixers, alignedDiodes, sampleRange] = frascatiCalProcessSignals( rawMixers, rawDiodes, rawMixerSensitivities, rawDiodeSensitivities, sampleRange, monIndices, useAlignment, useOffset )
%frascatiCalProcessSignals
%   Aligns calibration signals, converts them to volts.
%   Asks user to select what sample range to use.
    
    if (strcmp(getenv('USER'),'jack') == 1)
        addpath('/home/jack/Documents/MATLAB/PhaseFeedforward/ctfmod/PhaseFeedforward');
    else
        addpath('/nfs/cs-ccr-nfs6/vol29/share/ctf/data/ctfmod/MatLab/PhaseFeedforward');
    end
       
    [nMons, nScanPoints, nPulses, nSamples] = size(rawMixers);
    
    if (nargin<6 || isempty(monIndices))
        monIndices = 1:nMons;
    end
    if (nargin<7 || isempty(useAlignment))
        useAlignment = 1;
    end
    if (nargin<8 || isempty(useOffset))
        useOffset = 2;
    end
    
    % convert signals to volts, subtract offsets
    for mon=1:nMons
        if (~isempty(rawMixerSensitivities))
            rawMixers(mon,:,:,:) = rawMixers(mon,:,:,:).*rawMixerSensitivities(mon);
            if (~isempty(rawDiodes))
                rawDiodes(mon,:,:,:) = rawDiodes(mon,:,:,:).*rawDiodeSensitivities(mon);
            end
        end
        
        for point=1:nScanPoints
            for pulse=1:nPulses
                if (useOffset==1)
                    rawMixers(mon,point,pulse,:) = rawMixers(mon,point,pulse,:) + mixerOffsets(mon,p);   
                    if (~isempty(rawDiodes))
                        rawDiodes(mon,point,pulse,:) = rawDiodes(mon,point,pulse,:) + diodeOffsets(mon,p);
                    end                   
                end
            end
        end
    
    end
       
    % align signals
    %nPointsAligned = 300;
    %alignedMixers = NaN(nMons,nScanPoints,nPulses,nPointsAligned);
    %alignedDiodes = NaN(nMons,nScanPoints,nPulses,nPointsAligned);   
    if (useAlignment)
        if (isempty(rawDiodes))
            error('Diode must be present to use alignment function');
        end
        alignedDiodes = NaN(nMons,nScanPoints,nPulses,nSamples);
        alignedMixers = NaN(nMons,nScanPoints,nPulses,nSamples);

        for mon=monIndices
            for point=1:nScanPoints
                % align pulses
                tmpMixRaw = squeeze(rawMixers(mon,point,:,:));
                tmpDioRaw = squeeze(rawDiodes(mon,point,:,:));
                %[tmpDioA,tmpMixA] = getAligned2(tmpDioRaw,20,300,1,{tmpMixRaw}); 
                [tmpDioA,tmpMixA] = getAlignedXCorr(tmpDioRaw,'end',{tmpMixRaw});
                alignedDiodes(mon,point,:,:) = tmpDioA;
                alignedMixers(mon,point,:,:) = tmpMixA{1};
             end
        end
    else
        alignedMixers = rawMixers;
        alignedDiodes = rawDiodes;
    end
    
    figRefs = cell(1,nMons);
    if (~useAlignment && (isempty(sampleRange) || nargin<5))
        sampleRange = cell(1,nMons);
    end
    
    % plot signals, then ask for sample range to use
    for mon=monIndices
        figRefs{mon} = figure();
        for point=1:nScanPoints
            if (~isempty(alignedDiodes))
                subplot(1,2,1)
                plot(squeeze(alignedDiodes(mon,point,1,:)));
                hold all;
                subplot(1,2,2)
            end
            plot(squeeze(alignedMixers(mon,point,1,:)));
            hold all;
        end

        if (~isempty(alignedDiodes))
            subplot(1,2,1)
            title(sprintf('MIX %d DIODE',mon))
            xlabel('Sample No.')
            ylabel('Output [V]')
            subplot(1,2,2)
        end
        title(sprintf('MIX %d MIXER',mon))
        xlabel('Sample No.')
        ylabel('Output [V]')
        hold off;

        % different sample range for each monitor
        if (~useAlignment)
            if (isempty(sampleRange{mon}))
                fprintf('MONITOR %d\n',mon);
                startSample = input('Start calibration at sample: ');
                endSample = input('End calibration at sample: ');
                sampleRange{mon} = startSample:endSample;
                fprintf('------------------------------\n');
            end
        end

    end

    % same sample range for all monitors
    if (useAlignment)
        if (isempty(sampleRange))
            startSample = input('Start calibration at sample: ');
            endSample = input('End calibration at sample: ');
            sampleRange = startSample:endSample;
        end
    end
    

    % remove bad pulses
    for mon=monIndices
        if (useAlignment)
            tmpSampleRange = sampleRange;
        else
            tmpSampleRange = sampleRange{mon};
        end
        
        for point=1:nScanPoints
            tmpMix = squeeze(alignedMixers(mon,point,:,:));
            
            % cut on diode
            if (~isempty(alignedDiodes))
                tmpDio = squeeze(alignedDiodes(mon,point,:,:));
                [tmpDio,tmpMix] = removeBadPulses(tmpDio,tmpSampleRange,{tmpMix});
                alignedDiodes(mon,point,:,:) = tmpDio;
                tmpMix= tmpMix{1};
            end
            
            % cut on mixer
            tmpMix = removeBadPulses(tmpMix,tmpSampleRange);
            alignedMixers(mon,point,:,:) = tmpMix;

         end
    end

    % plot final results
    for mon=monIndices
        if (~isempty(figRefs{mon}))
            figure(figRefs{mon});
        else
            figRefs{mon} = figure();
        end
        
        if (useAlignment)
            tmpSampleRange = sampleRange;
        else
            tmpSampleRange = sampleRange{mon};
        end
        
        for point=1:nScanPoints
            if (~isempty(alignedDiodes))
                subplot(1,2,1)
                plot(squeeze(alignedDiodes(mon,point,1,:)));
                hold all;
                subplot(1,2,2)
            end
            plot(squeeze(alignedMixers(mon,point,1,:)));
            hold all;
        end
        
        if (~isempty(alignedDiodes))
            subplot(1,2,1)
            plot([tmpSampleRange(1) tmpSampleRange(1)],get(gca,'YLim'),'k')
            plot([tmpSampleRange(end) tmpSampleRange(end)],get(gca,'YLim'),'k')
            title(sprintf('MIX %d DIODE',mon))
            xlabel('Sample No.')
            ylabel('Output [V]')
            subplot(1,2,2)
        end
        plot([tmpSampleRange(1) tmpSampleRange(1)],get(gca,'YLim'),'k')
        plot([tmpSampleRange(end) tmpSampleRange(end)],get(gca,'YLim'),'k')
        title(sprintf('MIX %d MIXER',mon))
        xlabel('Sample No.')
        ylabel('Output [V]')
        hold off;
    end
    
end

