function [] = plotFitResult(mon,phases,mixers,diodes,useMixerOverSqrtDiode,sampleRange,calibrationFactors,phaseShiftFreq)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    [nMons,~,nPulses,~] = size(mixers);

    % Set frascatiOutput depending on whether diode should be used or not
    if (useMixerOverSqrtDiode)
        frascatiOutput = mixers./sqrt(abs(diodes));        
        % deal with any Infs coming from 0 diode output.
        frascatiOutput(frascatiOutput == Inf) = NaN;
        frascatiOutput(frascatiOutput == -Inf) = NaN;       
    else
        frascatiOutput = mixers;
    end
    
    frascatiOutput = frascatiOutput(:,:,:,sampleRange);
    meanFrascatiOutput = nanmean(nanmean(frascatiOutput,4),3);
    errFrascatiOutput = nanstd(nanmean(frascatiOutput,4),0,3)./sqrt(nPulses); 
    
%     for mon=1:nMons
        try
            monMean = squeeze(meanFrascatiOutput(mon,:,:,:));
            monErr = squeeze(errFrascatiOutput(mon,:,:,:));

            % remove any NaN entries (if there were points with only NaNs)
            goodPoints = ~isnan(monMean);
            goodPhases = phases(goodPoints);
            monMean = monMean(goodPoints);
            monErr = monErr(goodPoints);
           
            % sensible guesses for fit parameters
            maxFrascatiOutput = max(monMean);
            minFrascatiOutput = min(monMean);
            % sensible guesses for fit parameters
            initA = mean(abs([maxFrascatiOutput minFrascatiOutput])); % amplitude
            initB = (12/phaseShiftFreq)*(pi/180); % phase units, 4x for 3GHz -> 12GHz
            initD = (abs(maxFrascatiOutput) - abs(minFrascatiOutput))./2; % vertical offset

            % calculate suitable guess for sine phase shift based on offset
            % of first point from zero and gradient of calibration curve
            % between first and second points.
            sinArg = (monMean(1)-initD)./initA;
            if (sinArg>1); sinArg = 1; end; % avoid imaginary component if data point above amplitude guess
            if (sinArg<-1); sinArg = -1; end;
            initC = asin(sinArg)-(initB.*goodPhases(1));
            initC = mod(initC,2*pi);

            initParams = [initA, initB, initC, initD];
            initGuessFunc = @(x) initA.*sin(initB.*x + initC) + initD;
            
            fitFunc = @(x) (calibrationFactors(mon,1)).*sin( (calibrationFactors(mon,2)).*x + (calibrationFactors(mon,3)) ) + calibrationFactors(mon,4);

            figure;
            errorbar(goodPhases,monMean,monErr,'ob');
            hold all;
            fitPlotX = linspace(min(goodPhases),max(goodPhases),720);
            fitPlotY = fitFunc(fitPlotX);
            plot(fitPlotX, fitPlotY,'r','LineWidth',2);
            fitPlotY = initGuessFunc(fitPlotX);
            plot(fitPlotX, fitPlotY,'r--');
            title(sprintf('MIXER %d',mon));
            xlabel('Phase Shifter Setting');
            ylabel('Output [V]');
            grid on;
            legend('Data','Fit','Initial Guess','Location','best');
        catch
        end
%     end
end

