function [ calibrationFactors, zeroCrossings, calRSquare, calFitConf] = frascatiCalCalculateCalibrations( phases, mixers, diodes, sampleRange, useDiode, phaseShiftFreq, savePlotName, monIndices )
%frascatiCalCalculateCalibrations Calculated calibration factors from
%scan data.
%   This function should be called using the output of
%   frascatiCalProcessSignals.
%   useDiode=0: mixer only, useDiode=1: mixer/sqrt(diode), useDiode=2:
%   mixer/diode.
    if (nargin<6)
        phaseShiftFreq = 3;
    end
    if (nargin<7)
        plotResults = 1;
        savePlotName = '';
    elseif(isempty(savePlotName))
        plotResults=0;
    else
        plotResults = 1;
    end
    
    if (plotResults && ~isempty(savePlotName))
        saveFitPlot = 1;
    end
    

    [nMons, ~, nPulses, ~] = size(mixers);
    
    if (nargin<8)
        monIndices = 1:nMons;
    end
    
    if (plotResults)
        allMonFig = figure;
        mixCols = [0 0 1; 1 0 0 ;0 0.8 0];
        plotHand=[];
    end

    
    zeroCrossings = NaN(nMons,1);
    calibrationFactors = NaN(nMons,4);
    calRSquare = NaN(nMons,1);
    calFitConf = NaN(nMons,4);
    
    % Set frascatiOutput depending on whether diode should be used or not
    if (useDiode==1)
        frascatiOutput = mixers./sqrt(abs(diodes));
        
        % deal with any Infs coming from 0 diode output.
        frascatiOutput(frascatiOutput == Inf) = NaN;
        frascatiOutput(frascatiOutput == -Inf) = NaN;
    
    elseif (useDiode==2)
        frascatiOutput = mixers./abs(diodes);
        
        % deal with any Infs coming from 0 diode output.
        frascatiOutput(frascatiOutput == Inf) = NaN;
        frascatiOutput(frascatiOutput == -Inf) = NaN;
       
    else
        frascatiOutput = mixers;
    end
    
    % calculate means, fit weights
    if (~iscell(sampleRange))
        % one sample range for all monitors
        frascatiOutput = frascatiOutput(:,:,:,sampleRange);
        meanFrascatiOutput = nanmean(nanmean(frascatiOutput,4),3);
        errFrascatiOutput = nanstd(nanmean(frascatiOutput,4),0,3)./sqrt(nPulses);
        weightFrascatiOutput = 1./(nanstd(nanmean(frascatiOutput,4),0,3).^2);
    end
    
    for mon=monIndices
        try
            % calculate means, fit weights
            if (iscell(sampleRange))
                % different sample range for each monitor
                monFrascatiOutput = frascatiOutput(:,:,:,sampleRange{mon});
                meanFrascatiOutput = nanmean(nanmean(monFrascatiOutput,4),3);
                errFrascatiOutput = nanstd(nanmean(monFrascatiOutput,4),0,3)./sqrt(nPulses);
                weightFrascatiOutput = 1./(nanstd(nanmean(monFrascatiOutput,4),0,3).^2);   
            end
            
            monMean = squeeze(meanFrascatiOutput(mon,:,:,:));
            monWeight = squeeze(weightFrascatiOutput(mon,:,:,:));
            monErr = squeeze(errFrascatiOutput(mon,:,:,:));

            % remove any NaN entries (if there were points with only NaNs)
            goodPoints = ~isnan(monMean);
            goodPhases = phases(goodPoints);
            monWeight = monWeight(goodPoints);
            monMean = monMean(goodPoints);
            monErr = monErr(goodPoints);

            % find max and min points in scan (max and min on sine curve). 
            [maxFrascatiOutput,maxFrascatiOutputIndex] = max(monMean);
            [minFrascatiOutput,minFrascatiOutputIndex] = min(monMean);

            % find zero crossing of the fitted curve (guess half way between
            % min/max)
            guessZeroCrossIndex = abs(maxFrascatiOutputIndex-minFrascatiOutputIndex)./2; % index difference between min and max points
            guessZeroCrossIndex = guessZeroCrossIndex + min([minFrascatiOutputIndex maxFrascatiOutputIndex]); % plus the starting index
            guessZeroCrossIndex = round(guessZeroCrossIndex); % rounded to the nearest integer

            % sensible guesses for fit parameters
            initA = (maxFrascatiOutput-minFrascatiOutput)/2; % amplitude
            initB = (12/phaseShiftFreq)*(pi/180); % phase units, 4x for 3GHz -> 12GHz
            initD = (maxFrascatiOutput+minFrascatiOutput)/2; % vertical offset
            % calculate phase offset using other guesses for initial
            % parameters and phase/output at first point.
            sinArg = (monMean-initD)./initA;
            sinArg(sinArg>1)= NaN; % avoid imaginary component if data point above amplitude guess
            sinArg(sinArg<-1)= NaN;
            sinArg = asin(sinArg);
            % need to select 1st quadrant where asin valid - positive
            % gradient in difference
            sinArg(diff(monMean)<0) = NaN;
            initC = sinArg-(initB.*goodPhases);
            initC = mod(initC,2*pi); % make everything between 0 and 2pi
            initC(initC>pi) = initC(initC>pi)-2*pi;
            initC(initC<-pi) = initC(initC<-pi)+2*pi;
            initC = nanmean(initC);
            
%             gradFirst = sign(diff(monMean(1:2)));
%             signFirst = sign(monMean(1));
%             if (gradFirst<0)
%                 if (signFirst<0)
%                     initC = initC-(pi/2);
%                 else
%                     initC = initC + pi;
%                 end
%             end
%             initC = initC-(initB.*goodPhases(1));
%             initGuessFunc = @(x) initA.*sin(initB.*x + initC) + initD;
%             
%             gradientFit = sign(diff(initGuessFunc(goodPhases(1:2))));
%             if (gradientFit ~= gradientData)
%                 if (initC > 0)
%                     initC = initC-pi;
%                 else
%                     initC = initC+pi;
%                 end
%             end
            initParams = [initA, initB, initC, initD];

            % calculate fit
            [calibrationFactors(mon,:), calRSquare(mon), calFitConf(mon,:)] = offsetSinFit(goodPhases', monMean', initParams, monWeight);
            
            % save fit parameters
%             calibrationFactors(mon,1) = scanSinFit.a;
%             calibrationFactors(mon,2) = scanSinFit.b;
%             calibrationFactors(mon,3) = scanSinFit.c;
%             calibrationFactors(mon,4) = scanSinFit.d;
            fitFunc = @(x) (calibrationFactors(mon,1)).*sin( (calibrationFactors(mon,2)).*x + (calibrationFactors(mon,3)) ) + calibrationFactors(mon,4);
            initGuessFunc = @(x) initA.*sin(initB.*x + initC) + initD;

            % find zero crossing of the fitted curve (guess half way between
            % min/max)
            guessZeroCrossValue = goodPhases(guessZeroCrossIndex); % phase shifter value at this index
            % to get zero crossing use fit function without vertical offset,
            % otherwise zero crossing fit func ~= zero crossing sine curve
            % (actually looking for sin(x) = d, not sin(x) = 0)
            zeroFunc = @(x) (calibrationFactors(mon,1)).*sin( (calibrationFactors(mon,2)).*x + (calibrationFactors(mon,3)) ); 
            zeroCrossings(mon) = fzero(zeroFunc, guessZeroCrossValue);

            % always want to set the zero crossing to the falling slope so check
            % this here
            dydxAtZero = calibrationFactors(mon,1)*calibrationFactors(mon,2)*cos(calibrationFactors(mon,2)*zeroCrossings(mon) + calibrationFactors(mon,3)); % Ab*cos(bx+c)
            if (dydxAtZero > 0)
                zeroCrossings(mon) = zeroCrossings(mon) + (180/(12/phaseShiftFreq)); % divide by 4 because of 3 GHz phase shifter vs. 12 GHz mixer
            end
        
            
            % plot results
            if (plotResults)
                figure;
                %plot(goodPhases,monMean,'.b');
                errorbar(goodPhases,monMean,monErr,'ob');
                hold all;
                fitPlotX = linspace(min(goodPhases),max(goodPhases),720);
                fitPlotY = initGuessFunc(fitPlotX);
                plot(fitPlotX, fitPlotY,'r--');
                fitPlotY = fitFunc(fitPlotX);
                plot(fitPlotX, fitPlotY,'r','LineWidth',2);                
                plot([zeroCrossings(mon) zeroCrossings(mon)], [minFrascatiOutput maxFrascatiOutput],'k','LineWidth',2);
                title(sprintf('MIXER %d',mon));
                xlabel('Phase shifter setting');
                ylabel('Output');
                grid on;
                legend('Data','Initial Guess','Fit','Zero Crossing');
                
                if (saveFitPlot)
                    savePlot(sprintf('%s_Mon%d',savePlotName,mon),'',[1 0 0]);
                end
                
                figure(allMonFig);
                tmpHand=errorbar(goodPhases,monMean,monErr,'o','Color',mixCols(mon,:));
                plotHand = [plotHand tmpHand];
                hold all;
                plot(fitPlotX, fitPlotY,'Color',mixCols(mon,:),'LineWidth',1.5);                

            end
    
        catch e
            warning('Calibration of Monitor %d failed. %s',mon,e.message);
        end
    end
    
    if (plotResults)
        figure(allMonFig)
        title('Calibration Results')
        ylabel('Mixer Output [V]')
        xlabel('Phase Shifter Setting')
        xlim([min(phases) max(phases)])
        plotLab = {'Mon 1','Mon 2','Mon 3'};
        tmpInd = 1:nMons;
        plotLab = plotLab(tmpInd(~isnan(squeeze(calibrationFactors(:,1))))');
        legend(plotHand,plotLab);
%         format_plots;
    end
end

