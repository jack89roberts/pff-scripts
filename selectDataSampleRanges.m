function [ sampleRanges ] = selectDataSampleRanges( processedPhaseData )
%selectDataSampleRanges Makes some plots for user to select sample range to
%use for each monitor.
%   Detailed explanation goes here

    frascatiNMonitors = processedPhaseData.frascatiNMonitors;
    frascatiStartSamples = NaN*ones(1,frascatiNMonitors);
    frascatiEndSamples = NaN*ones(1,frascatiNMonitors);
  
    sampleSelectFigureHandle = figure();
    set(0,'CurrentFigure', sampleSelectFigureHandle);      

    frascatiGuessStartSamples = NaN*ones(frascatiNMonitors);
    frascatiGuessEndSamples = NaN*ones(frascatiNMonitors);
    
    for mon=1:frascatiNMonitors
        
        [frascatiGuessStartSamples(mon),frascatiGuessEndSamples(mon)] = getSignalRange(processedPhaseData.frascatiDiodes(mon,1,:), 10);    
        
        subplot(1,3,1);
        plot(squeeze(processedPhaseData.frascatiDiodes(mon,1,:)));
        title(sprintf('Frascati Mon %d Diode',mon));
        xlabel('Sample No.');
        ylabel('Output [V]');

        subplot(1,3,2);
        plot(squeeze(processedPhaseData.frascatiMixers(mon,1,:)));
        title(sprintf('Frascati Mon %d Mixer',mon));
        xlabel('Sample No.');
        ylabel('Output [V]');
        
        subplot(1,3,3);
        plot(squeeze(processedPhaseData.frascatiPhases(mon,1,:)));
        title(sprintf('Frascati Mon %d Phase',mon));
        xlabel('Sample No.');
        ylabel('Output [degrees]');

        fprintf(1, 'Monitor %d\n',mon);
        frascatiStartSamples(mon) = input('Start calibration at sample: ');
        frascatiEndSamples(mon) = input('End calibration at sample: ');
        fprintf(1,'------------------------------------\n');
    end

    subplot(1,3,1)
    plot(squeeze(processedPhaseData.petsI(1,:)));
    title('PETS I');
    xlabel('Sample No.');
    ylabel('Output [V]');

    subplot(1,3,2)
    plot(squeeze(processedPhaseData.petsQ(1,:)));
    title('PETS Q');
    xlabel('Sample No.');
    ylabel('Output [V]');

    subplot(1,3,2)
    plot(squeeze(processedPhaseData.petsPhase(1,:)));
    title('PETS Phase');
    xlabel('Sample No.');
    ylabel('Output [degrees]');

    
    fprintf(1, 'PETS Phase\n');
    petsStartSample = input('Start calibration at sample: ');
    petsEndSample = input('End calibration at sample: ');

    close(sampleSelectFigureHandle);
    
    sampleRanges = struct();
    sampleRanges.frascatiStartSamples = frascatiStartSamples;
    sampleRanges.frascatiEndSamples = frascatiEndSamples;
    sampleRanges.petsStartSample = petsStartSample;
    sampleRanges.petsEndSample = petsEndSample;

end

