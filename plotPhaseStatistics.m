function [ ] = plotPhaseStatistics( processedData, saveDir )
%plotPhaseStatistics Makes a few plots from phase statistics struct in
%processedData
%   vargin(1) : saveDir -> will save plots if given

    if nargin<2
        savePlots = false;
        
    else
        if(~exist(saveDir,'dir'))
            mkdir(saveDir);
        end
        savePlots = true;
    end

    phaseStatistics = processedData.phaseStatistics;
    frascatiNMonitors = processedData.phaseData.frascatiNMonitors;
    frascatiTimeAxes = processedData.phaseData.frascatiTimeAxes;
    petsTimeAxis = processedData.phaseData.petsTimeAxis;
    frascatiStartSamples = processedData.phaseData.frascatiStartSamples;
    frascatiEndSamples = processedData.phaseData.frascatiEndSamples;
    petsStartSample = processedData.phaseData.petsStartSample;
    petsEndSample = processedData.phaseData.petsEndSample;
    dataName = strrep(processedData.dataName,'_',' ');
    
    frascatiNames = {'Mon1','Mon2','Mon3'};
    petsName = 'PETS';
    
    % mean phase plots
    figure;            
    for mon=1:frascatiNMonitors
        tmpPhase = phaseStatistics.frascatiMeanPulsePhases(mon,:);
        plot(squeeze(tmpPhase));
        hold on
    end
    plot(phaseStatistics.petsMeanPulsePhase);
    tmpLegStr = cell(1,frascatiNMonitors+1);
    for mon=1:frascatiNMonitors
        tmpLegStr{mon} = [frascatiNames{mon} sprintf(' (std %.1f^{o})', phaseStatistics.frascatiStdPulsePhases(mon))];
    end
    tmpLegStr{frascatiNMonitors+1} = [petsName sprintf(' (std %.1f^{o})', phaseStatistics.petsStdPulsePhase)];
    legend(tmpLegStr,'Location','best');
    title({dataName, 'Mean Pulse Phase vs. Time'});
    xlabel('Time [pulse no.]');
    ylabel('Phase [12GHz degrees]');
    grid
    hold off
    if (savePlots)
        saveStr = [saveDir '/meanPulsePhaseVsTime'];
        print([saveStr '.png'],'-dpng');
    end

    figure;
    for mon=1:frascatiNMonitors
        tmpPhase = phaseStatistics.frascatiMeanSamplePhases(mon,frascatiStartSamples(mon):frascatiEndSamples(mon));
        tmpTimeAxis = frascatiTimeAxes{mon};
        tmpTimeAxis = tmpTimeAxis(frascatiStartSamples(mon):frascatiEndSamples(mon));
        plot(tmpTimeAxis,squeeze(tmpPhase));
        hold on
    end
    tmpTimeAxis = petsTimeAxis(petsStartSample:petsEndSample);
    plot(tmpTimeAxis,phaseStatistics.petsMeanSamplePhase(petsStartSample:petsEndSample));

    title({dataName,'Mean Phase Along the Pulse'});
    xlabel('Time [ns]');
    ylabel('Phase [12GHz degrees]');
    legend([frascatiNames petsName],'Location','best');
    grid
    hold off
    if (savePlots)
        saveStr = [saveDir '/meanSamplePhase'];
        print([saveStr '.png'],'-dpng');
    end

    figure;
    for mon=1:frascatiNMonitors
        tmpPhase = phaseStatistics.frascatiStdSamplePhases(mon,frascatiStartSamples(mon):frascatiEndSamples(mon));
        tmpTimeAxis = frascatiTimeAxes{mon};
        tmpTimeAxis = tmpTimeAxis(frascatiStartSamples(mon):frascatiEndSamples(mon));
        plot(tmpTimeAxis,squeeze(tmpPhase));
        hold on
    end
    tmpTimeAxis = petsTimeAxis(petsStartSample:petsEndSample);
    plot(tmpTimeAxis,phaseStatistics.petsStdSamplePhase(petsStartSample:petsEndSample));

    tmpLegStr = cell(1,frascatiNMonitors+1);
    for mon=1:frascatiNMonitors
        tmpLegStr{mon} = [frascatiNames{mon} sprintf(' (mean %.1f^{o})', phaseStatistics.frascatiMeanStdSamplePhases(mon))];
    end
    tmpLegStr{frascatiNMonitors+1} = [petsName sprintf(' (mean %.1f^{o})', phaseStatistics.petsMeanStdSamplePhase)];
    legend(tmpLegStr,'Location','best');%legend([frascatiNames petsName]);
    title({dataName,'Std Phase Along the Pulse'});
    xlabel('Time [ns]');
    ylabel('Phase [12GHz degrees]');
    grid
    hold off
    if (savePlots)
        saveStr = [saveDir '/stdSamplePhase'];
        print([saveStr '.png'],'-dpng');
    end

    % correlation plots
    figure;
    tmpPhaseA = phaseStatistics.frascatiMeanPulsePhases(1,:);
    tmpPhaseB = phaseStatistics.frascatiMeanPulsePhases(2,:);
    plot(squeeze(tmpPhaseA),squeeze(tmpPhaseB),'o');  
    tmpCorr = phaseStatistics.corr12_meanPulsePhase;
    tmpGrad = phaseStatistics.fit12_meanPulsePhase(1);
    titString = sprintf('Mon1 vs. Mon2: corr=%.2f, grad=%.2f',tmpCorr,tmpGrad);
    title({dataName,titString});
    xlabel('Mon1 Phase [12GHz degrees]');
    ylabel('Mon2 Phase [12GHz degrees]');
    grid
    if (savePlots)
        saveStr = [saveDir '/corr12_meanPulsePhase'];
        print([saveStr '.png'],'-dpng');
    end
    
    figure;
    tmpPhaseA = phaseStatistics.frascatiMeanPulsePhases(1,:);
    tmpPhaseB = phaseStatistics.frascatiMeanPulsePhases(3,:);
    plot(squeeze(tmpPhaseA),squeeze(tmpPhaseB),'o');
    tmpCorr = phaseStatistics.corr13_meanPulsePhase;
    tmpGrad = phaseStatistics.fit13_meanPulsePhase(1);
    titString = sprintf('Mon1 vs. Mon3: corr=%.2f, grad=%.2f',tmpCorr,tmpGrad);
    title({dataName,titString});
    xlabel('Mon1 Phase [12GHz degrees]');
    ylabel('Mon3 Phase [12GHz degrees]');
    grid
    if (savePlots)
        saveStr = [saveDir '/corr13_meanPulsePhase'];
        print([saveStr '.png'],'-dpng');
    end
    
    figure;
    tmpPhaseA = phaseStatistics.frascatiMeanPulsePhases(3,:);
    tmpPhaseB = phaseStatistics.petsMeanPulsePhase;
    plot(squeeze(tmpPhaseA),squeeze(tmpPhaseB),'o');
    tmpCorr = phaseStatistics.corrPetsFrascati3_meanPulsePhase;
    tmpGrad = phaseStatistics.fitPetsFrascati3_meanPulsePhase(1);
    titString = sprintf('Mon3 vs. PETS: corr=%.2f, grad=%.2f',tmpCorr,tmpGrad);
    title({dataName,titString});
    xlabel('Mon3 Phase [12GHz degrees]');
    ylabel('PETS Phase [12GHz degrees]');
    grid
    if (savePlots)
        saveStr = [saveDir '/corrPetsFrascati3_meanPulsePhase'];
        print([saveStr '.png'],'-dpng');
    end
    

end

