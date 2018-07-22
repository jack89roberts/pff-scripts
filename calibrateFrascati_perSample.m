% Calibration per sample point using previously saved data file from
% calibrateFrascati.m script.

dataFile = 'data/frascatiCalibrations/frascatiCalibration_20141002_1453.mat';

%%

load(dataFile);

% some plots to help pick sample range
for i=1:nMons
    figure;
    subplot(1,2,1);
    for n=1:scanNPoints
        plot(squeeze(scanMixers(i,n,1,:)));
        title(sprintf('Monitor %d Mixer',i));
        xlabel('Sample No.');
        ylabel('Output [V]');
        hold all;
    end
    subplot(1,2,2);
    for n=1:scanNPoints
        plot(squeeze(scanDiodes(i,n,1,:)));
        title(sprintf('Monitor %d Diode',i));
        xlabel('Sample No.');
        ylabel('Output [V]');
        hold all;
    end
end

%% Calibrate all samples

fprintf(1, 'Calibration with all Samples in range:\n');
calStartSample = input('Start calibration at sample: ');
calEndSample = input('End calibration at sample: ');

nCalSamples = length(calStartSample:calEndSample);

zeroCrossings = NaN*ones(nMons,nSamples);
calibrationFactors = NaN*ones(nMons,nSamples,4);

for s=calStartSample:calEndSample

    if (useMixerOverSqrtDiode == true)
        meanFrascatiOutput = nanmean(scanMixerSqrtDiode(:,:,:,s),3);
        stdFrascatiOutput = nanstd(scanMixerSqrtDiode(:,:,:,s),0,3); % NOTE - this not currently used for anything. Could/should use it for fitting weights.
    else
        meanFrascatiOutput = nanmean(scanMixers(:,:,:,s),3);
        stdFrascatiOutput = nanstd(scanMixers(:,:,:,s),0,3); % NOTE - this not currently used for anything. Could/should use it for fitting weights.
    end

  
    
    % find max and min points in scan (max and min on sine curve). 
    [maxFrascatiOutput,maxFrascatiOutputIndex] = max(meanFrascatiOutput,[],2);
    [minFrascatiOutput,minFrascatiOutputIndex] = min(meanFrascatiOutput,[],2);

    for mon=1:nMons
        % sensible guesses for fit parameters
        initA = mean(abs([maxFrascatiOutput(mon) minFrascatiOutput(mon)])); % amplitude
        initB = 4*(pi./180); % phase units, 4x for 3GHz -> 12GHz
        initC = (meanFrascatiOutput(mon,1)./initA)*(pi/2); % phase shift
        initD = (abs(maxFrascatiOutput(mon)) - abs(minFrascatiOutput(mon)))./2; % vertical offset
        initParams = [initA, initB, initC, initD];

        % calculate fit
        [scanSinFit,scanSinGoF] = offsetSinFit(scanPhShiftValues', meanFrascatiOutput', initParams);

        % save fit parameters
        calibrationFactors(mon,s,1) = scanSinFit.a;
        calibrationFactors(mon,s,2) = scanSinFit.b;
        calibrationFactors(mon,s,3) = scanSinFit.c;
        calibrationFactors(mon,s,4) = scanSinFit.d;
        
        % find zero crossing of the fitted curve (guess half way between
        % min/max)
        guessZeroCrossIndex = abs(maxFrascatiOutputIndex-minFrascatiOutputIndex)./2; % index difference between min and max points
        guessZeroCrossIndex = guessZeroCrossIndex + min([minFrascatiOutputIndex maxFrascatiOutputIndex]); % plus the starting index
        guessZeroCrossIndex = round(guessZeroCrossIndex); % rounded to the nearest integer
        guessZeroCrossValue = scanPhShiftValues(guessZeroCrossIndex); % phase shifter value at this index
        % to get zero crossing use fit function without vertical offset,
        % otherwise zero crossing fit func ~= zero crossing sine curve
        % (actually looking for sin(x) = d, not sin(x) = 0)
        zeroFunc = @(x) (calibrationFactors(mon,s,1)).*sin( (calibrationFactors(mon,s,2)).*x + (calibrationFactors(mon,s,3)) ); 
        zeroCrossings(mon,s) = fzero(zeroFunc, guessZeroCrossValue);     

    end

end

zeroCrossings = zeroCrossings(:,calStartSample:calEndSample);
calibrationFactors = calibrationFactors(:,calStartSample:calEndSample,:);

% plot results
monitorNames = {'Monitor 1','Monitor 2', 'Monitor 3'};
calibFactorNames = {'a','b','c','d'};

for c=1:4 % calibration factor index
    figure;
    for mon=1:nMons
        plot(calStartSample:calEndSample, squeeze(calibrationFactors(mon,:,c)));
        hold all;
    end
    title(calibFactorNames{c});
    legend(monitorNames);
    xlabel('Sample No.');
    ylabel('Parameter Value');
end
figure;
for mon=1:nMons
    plot(calStartSample:calEndSample, squeeze(zeroCrossings(mon,:)));
    hold all;
end
title('Zero Crossing');
legend(monitorNames);
xlabel('Sample No.');
ylabel('Phase [degrees]');