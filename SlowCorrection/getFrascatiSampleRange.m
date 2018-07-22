function sampleRange = getFrascatiSampleRange( mixers,diodes )
%getFrascatiSampleRange Prompt for frascati sample range to use
%   mixers and diodes should contain only one pulse, arrays in format
%   monIndex,sampleNo.

    if (nargin<2)
        nMons = 3;
        nSamples = [];
        signals = makePropertyListFromFiles('../devices/scope.devs','../devices/scope.props');
        for i=1:nMons
            tmpMixer = matlabJapc.staticGetSignal('SCT.USER.SETUP',signals{2*i - 1});
            tmpDiode = matlabJapc.staticGetSignal('SCT.USER.SETUP',signals{2*i});
            
            if (isempty(nSamples))
                nSamples = length(tmpMixer);
                mixers = NaN(nMons,nSamples);
                diodes = NaN(nMons,nSamples);
            end
            
            mixers(i,:) = tmpMixer;
            diodes(i,:) = tmpDiode;
        end
    else
        [nMons,~] = size(mixers);   
    end

    figure;
    
    
    
    for i=1:nMons
%         [tmpDiode,tmpMixer] = getAlignedXCorr(diodes(i,:),'end',{mixers(i,:)});
%         tmpMixer = tmpMixer{1};
        tmpDiode = diodes(i,:);
        tmpMixer = mixers(i,:);
        subplot(1,2,1)
        plot(tmpDiode);
        hold all;
        title('DIODES');
        xlabel('Sample No,');
        ylabel('Output [mV]');
        subplot(1,2,2);
        plot(tmpMixer);
        hold all;
        title('MIXERS');
        xlabel('Sample No,');
        ylabel('Output [mV]');
    end

    startSamp = input('Start sample: ');
    endSamp = input('End sample: ');
    sampleRange = startSamp:endSamp;
    
    subplot(1,2,1);
    yLim = get(gca,'YLim');
    plot([startSamp startSamp],yLim,'k');
    plot([endSamp endSamp],yLim,'k');
    
    subplot(1,2,2);
    yLim = get(gca,'YLim');
    plot([startSamp startSamp],yLim,'k');
    plot([endSamp endSamp],yLim,'k');

    
end

