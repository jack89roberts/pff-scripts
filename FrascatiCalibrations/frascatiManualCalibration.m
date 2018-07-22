close all; clear all;
%%
dataBaseDir = '/home/jack/PhaseFeedforward/CTFData/201504/FrascatiCalibrations';
dataSetName = '20150429_1153_MechanicalShifterOnMixer1_Calibration';
sampleRange = 49:51; %49:51; % if this is left empty, will prompt to ask for sample range
phaseChangePulses = [];%[16, 24, 31, 36, 45, 49, 61, 66, 73, 80, 88, 95, 101, 109, 115, 124, 129]; % sample numbers when phase was changed. if this is left empty, will prompt to ask for samples
    
quickCal = true; % just using min/max value.
saveData = true; % save calibration results

useMixerOverSqrtDiode = 0; % NB: This script only works correctly with useMixerOverSqrtDiode=0.

%% load data
addpath('/home/jack/Documents/MATLAB/PhaseFeedforward/ctfmod/PhaseFeedforward');
[CTFData, ~, dataDir] = loadMergedData(dataSetName,'',dataBaseDir,'');

%% signal processing

[ mixers, diodes ] = extractMixerDiode(CTFData);
[nMons,nPulses,nSamples] = size(mixers);

% align mixers, diodes, phases
for mon=1:nMons
    Diode = squeeze(diodes(mon,:,:));
    Mixer = squeeze(mixers(mon,:,:));
    
    [Diode,Mixer] = getAligned2(Diode,20,nSamples,1,{Mixer});
    diodes(mon,:,:) = Diode;
    mixers(mon,:,:) = Mixer{1};
end

% prompt for sample range to use if not given
if (isempty(sampleRange))
    a = figure();
    
    for mon=1:nMons
        subplot(1,2,1)
        plot(squeeze(mixers(mon,1,:)));
        hold all;
        subplot(1,2,2)
        plot(squeeze(diodes(mon,1,:)));
        hold all;
    end
    
    subplot(1,2,1);
    title('MIXERS');
    subplot(1,2,2);
    title('DIODES')
    
    startSample = input('Start sample: ');
    endSample = input('End sample: ');
    sampleRange = startSample:endSample;

end

% remove duplicate pulses
for mon=1:nMons
    Diode = squeeze(diodes(mon,:,:));
    Mixer = squeeze(mixers(mon,:,:));

    diodes(mon,:,:) = removeDuplicatePulses(Diode);
    mixers(mon,:,:) = removeDuplicatePulses(Mixer);
end

% remove bad pulses: outside 3 sigma on diode
for mon=1:nMons
    Diode = squeeze(diodes(mon,:,:));
    Mixer = squeeze(mixers(mon,:,:));
    
    [Diode,Mixer] = removeBadPulses(Diode,sampleRange,{Mixer});
    diodes(mon,:,:) = Diode;
    mixers(mon,:,:) = Mixer{1};
end

% quick plot of scan
figure;
plot(squeeze(nanmean(mixers(:,:,sampleRange),3))');
xlabel('Pulse No.')
ylabel('Mean Output [mV]')
legend('Mon1','Mon2','Mon3');

%%
if (~quickCal)
    %% split the mixer output in to the different phase shifter settings

    % prompt for pulse numbers of phase changes
    % pulse numbers entered in input should be the first pulse after the
    % phase change.
    if (isempty(phaseChangePulses))
        phaseChangePulses = input('Phase change samples ([a,b,c] format): ');
    end

    nPoints = length(phaseChangePulses);
    frascatiOutput = NaN(nMons,nPoints,nPulses,nSamples); % very memory inefficient, but allow for each phase shifter setting being up to nPulses long

    for p=1:nPoints

        if (p==1)
            startPulse = 1;
        else
            startPulse = phaseChangePulses(p-1); 
        end
        endPulse = phaseChangePulses(p)-1;

        tmpMixers = mixers(:,startPulse:endPulse,:);
        [~,nPulsesPerPoint,~] = size(tmpMixers);

        frascatiOutput(:,p,1:nPulsesPerPoint,:) = tmpMixers(:,:,:);

    end

    %% calculate fit
    meanFrascatiOutput = nanmean(nanmean(frascatiOutput(:,:,:,sampleRange),4),3);
    errFrascatiOutput = nanstd(nanmean(frascatiOutput,4),0,3)./sqrt(nPulses);
    weightFrascatiOutput = 1./(nanstd(nanmean(frascatiOutput,4),0,3).^2);

    calibrationFactors = NaN(nMons,4);
    zeroCrossings = NaN(1,nMons);

    for mon=1:nMons
        monMean = squeeze(meanFrascatiOutput(mon,:,:,:));
        monWeight = squeeze(weightFrascatiOutput(mon,:,:,:));
        monErr = squeeze(errFrascatiOutput(mon,:,:,:));

        % find max and min points in scan (max and min on sine curve). 
        [maxFrascatiOutput,maxFrascatiOutputIndex] = max(monMean);
        [minFrascatiOutput,minFrascatiOutputIndex] = min(monMean);

        nPoints180Deg = abs(maxFrascatiOutputIndex-minFrascatiOutputIndex);
        nDegsPerPoint = 180/nPoints180Deg;
        shifterPhases = nDegsPerPoint:nDegsPerPoint:(nDegsPerPoint*nPoints);

        % find zero crossing of the fitted curve (guess half way between
        % min/max)
        guessZeroCrossIndex = abs(maxFrascatiOutputIndex-minFrascatiOutputIndex)./2; % index difference between min and max points
        guessZeroCrossIndex = guessZeroCrossIndex + min([minFrascatiOutputIndex maxFrascatiOutputIndex]); % plus the starting index
        guessZeroCrossIndex = round(guessZeroCrossIndex); % rounded to the nearest integer

        % sensible guesses for fit parameters
        initA = mean(abs([maxFrascatiOutput minFrascatiOutput])); % amplitude
        initB = pi./180; % phase units, 4x for 3GHz -> 12GHz
        initC = (monMean(1)./initA)*(pi/2); % phase shift
        initD = (abs(maxFrascatiOutput) - abs(minFrascatiOutput))./2; % vertical offset
        initParams = [initA, initB, initC, initD];

        try
            % calculate fit
            [scanSinFit,calibrationGoF] = offsetSinFit(shifterPhases', monMean', initParams, monWeight);

            % save fit parameters
            calibrationFactors(mon,1) = scanSinFit.a;
            calibrationFactors(mon,2) = scanSinFit.b;
            calibrationFactors(mon,3) = scanSinFit.c;
            calibrationFactors(mon,4) = scanSinFit.d;
            fitFunc = @(x) (calibrationFactors(mon,1)).*sin( (calibrationFactors(mon,2)).*x + (calibrationFactors(mon,3)) ) + calibrationFactors(mon,4);

            % find zero crossing of the fitted curve (guess half way between
            % min/max)
            guessZeroCrossValue = shifterPhases(guessZeroCrossIndex); % phase shifter value at this index
            % to get zero crossing use fit function without vertical offset,
            % otherwise zero crossing fit func ~= zero crossing sine curve
            % (actually looking for sin(x) = d, not sin(x) = 0)
            zeroFunc = @(x) (calibrationFactors(mon,1)).*sin( (calibrationFactors(mon,2)).*x + (calibrationFactors(mon,3)) ); 
            zeroCrossings(mon) = fzero(zeroFunc, guessZeroCrossValue);

            % always want to set the zero crossing to the falling slope so check
            % this here
            dydxAtZero = calibrationFactors(mon,1)*calibrationFactors(mon,2)*cos(calibrationFactors(mon,2)*zeroCrossings(mon) + calibrationFactors(mon,3)); % Ab*cos(bx+c)
            if (dydxAtZero > 0)
                zeroCrossings(mon) = zeroCrossings(mon) + (180/4); % divide by 4 because of 3 GHz phase shifter vs. 12 GHz mixer
            end


            % plot results
            figure;
            %plot(phases,monMean,'.b');
            errorbar(shifterPhases,monMean,monErr,'ob');
            hold all;
            plot([zeroCrossings(mon) zeroCrossings(mon)], [minFrascatiOutput maxFrascatiOutput],'k','LineWidth',2);
            fitPlotX = linspace(min(shifterPhases),max(shifterPhases),720);
            fitPlotY = fitFunc(fitPlotX);
            plot(fitPlotX, fitPlotY,'r','LineWidth',2);
            title(sprintf('MIXER %d',mon));
            xlabel('Phase shifter setting (degrees)');
            ylabel('Frascati Output');
            grid on;
            legend('Data','Zero Crossing','Fit');

        catch e
            warning('Calibration of Monitor %d failed. %s',mon,e.message);
        end
    end
    
else     
%% quick calibration: find max and min, use to calculate amplitude and offset
    maxPulses = cell(1,nMons);
    minPulses = cell(1,nMons);
    for mon=1:nMons
        fprintf('--------------------------------------------------\n');
        fprintf('MONITOR %d\n',mon);
        fprintf('--------------------------------------------------\n');
        maxPulses{mon} = input('Pulse window max output (a:b): ');
        minPulses{mon} = input('Pulse window min output (a:b): ');
    end
    fprintf('--------------------------------------------------\n');

    calAmp = NaN(1,nMons);
    calOffset = NaN(1,nMons);
    
    for mon=1:nMons
        meanMaxOutput = nanmean(nanmean(mixers(:,maxPulses{mon},sampleRange),3),2);
        meanMinOutput = nanmean(nanmean(mixers(:,minPulses{mon},sampleRange),3),2);

        calAmp(mon) = (meanMaxOutput(mon)-meanMinOutput(mon))/2;
        calOffset(mon) = (meanMaxOutput(mon)+meanMinOutput(mon))/2;
    end
    
    calibrationFactors = NaN(nMons,4);
    calibrationFactors(:,1) = calAmp;
    calibrationFactors(:,4) = calOffset
    zeroCrossings = NaN(1,nMons);
end

%% save data
if (saveData)
    fprintf('Saving data...\n');

    calibConstName = sprintf('%s%s%s',dataBaseDir,'/frascatiCalibrationConstants_',dataSetName);
    calibConstFileBackup = fopen(calibConstName,'w');
    for mon=1:nMons 
        fprintf(calibConstFileBackup,'Mixer%d, %f, %f, %f, %f\n',...
            mon,...
            calibrationFactors(mon,1),...
            calibrationFactors(mon,2),...
            calibrationFactors(mon,3),...
            calibrationFactors(mon,4)...
        );
    end
    fprintf(calibConstFileBackup, 'useMixerOverSqrtDiode, %d', useMixerOverSqrtDiode);
    fclose(calibConstFileBackup);
end