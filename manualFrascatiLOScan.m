dataDir = '/nfs/cs-ccr-nfs6/vol29/share/ctf/data/ctfmod/MatLab/PhaseFeedforward/data/frascatiCalibrations/20141002_1635_manualShifterCal';

useMixerOverSqrtDiode = true;

phaseShifterFreq = 2; % frequency of phase shifter degrees in GHz, for some reason needed 2 GHz not 1 GHz here for manual phase shifter to get fit to work??

mixerName = 'CT.SCOPE01.CH04.Acquisition.value.value';
diodeName = 'CT.SCOPE01.CH03.Acquisition.value.value';
mixerSensitivityName = 'CT.SCOPE01.CH04.Acquisition.sensitivity.value';
diodeSensitivityName = 'CT.SCOPE01.CH03.Acquisition.sensitivity.value';

nSamples = 500;

%%
phaseMultiplier = 12./phaseShifterFreq;

CTFData = mergeMatMonData(dataDir);
nPulses = length(CTFData);

% assumes sensitivity does not change within data (just use first pulse)
eval(['mixerSensitivity = CTFData(1).' mixerSensitivityName ';']);
eval(['diodeSensitivity = CTFData(1).' diodeSensitivityName ';']);

scanMixers = NaN*ones(nPulses,nSamples);
scanDiodes = NaN*ones(nPulses,nSamples);
scanPhShiftValues = NaN*ones(1,nPulses);

for t=1:nPulses
    
    % take comment of form: Name, value and put value in kickSetting
    pulseComment = CTFData(t).comment;
    splitComment = regexp(pulseComment, ',', 'split');
    scanPhShiftValues(t) = str2double(splitComment{2});
    
    eval(['scanMixers(t,:) = CTFData(t).' mixerName ';']);
    eval(['scanDiodes(t,:) = CTFData(t).' diodeName ';']);
end

% convert in to volts and calculate mixer/sqrt(diode)
% convert to volts
scanMixers = scanMixers.*mixerSensitivity;
scanDiodes = scanDiodes.*diodeSensitivity;

% subtract offset mean first 20 samples
for t=1:nPulses
    scanMixers(t,:) = scanMixers(t,:) - mean(scanMixers(t,1:20));    
    scanDiodes(t,:) = scanDiodes(t,:) - mean(scanDiodes(t,1:20));
end

scanMixerSqrtDiode = scanMixers./abs(sqrt(scanDiodes));

% subplot(1,2,1)
% scatter(scanPhShiftValues,squeeze(scanMixers(:,150)))
% subplot(1,2,2)
% scatter(scanPhShiftValues,squeeze(scanDiodes(:,150)))



%%

figure;
subplot(1,2,1)
plot(scanMixers(1,:))
subplot(1,2,2)
plot(scanDiodes(1,:))

calStartSample = input('Start calibration at sample: ');
calEndSample = input('End calibration at sample: ');

calibrationFactors = NaN*ones(1,4);

if (useMixerOverSqrtDiode == true)
    meanFrascatiOutput = nanmean(scanMixerSqrtDiode(:, calStartSample:calEndSample),2);
else
    meanFrascatiOutput = nanmean(scanMixers(:, calStartSample:calEndSample),2);
end

% find max and min points in scan (max and min on sine curve). 
[maxFrascatiOutput,maxFrascatiOutputIndex] = max(meanFrascatiOutput);
[minFrascatiOutput,minFrascatiOutputIndex] = min(meanFrascatiOutput);
    

    
% sensible guesses for fit parameters
initA = mean(abs([maxFrascatiOutput minFrascatiOutput])); % amplitude
initB = phaseMultiplier*(pi./180); % phase units, 4x for 3GHz -> 12GHz
initC = (meanFrascatiOutput(1)./initA)*(pi/2); % phase shift
initD = (abs(maxFrascatiOutput) - abs(minFrascatiOutput))./2; % vertical offset
initParams = [initA, initB, initC, initD];

% calculate fit
[scanSinFit,scanSinGoF] = offsetSinFit(scanPhShiftValues', meanFrascatiOutput, initParams);

% save fit parameters
calibrationFactors(1) = scanSinFit.a;
calibrationFactors(2) = scanSinFit.b;
calibrationFactors(3) = scanSinFit.c;
calibrationFactors(4) = scanSinFit.d;
calibrationFactors
fitFunc = @(x) (calibrationFactors(1)).*sin( (calibrationFactors(2)).*x + (calibrationFactors(3)) ) + calibrationFactors(4);

% find zero crossing of the fitted curve (guess half way between
% min/max)
guessZeroCrossIndex = abs(maxFrascatiOutputIndex-minFrascatiOutputIndex)./2; % index difference between min and max points
guessZeroCrossIndex = guessZeroCrossIndex + min([minFrascatiOutputIndex maxFrascatiOutputIndex]); % plus the starting index
guessZeroCrossIndex = round(guessZeroCrossIndex); % rounded to the nearest integer
guessZeroCrossValue = scanPhShiftValues(guessZeroCrossIndex); % phase shifter value at this index
% to get zero crossing use fit function without vertical offset,
% otherwise zero crossing fit func ~= zero crossing sine curve
% (actually looking for sin(x) = d, not sin(x) = 0)
zeroFunc = @(x) (calibrationFactors(1)).*sin( (calibrationFactors(2)).*x + (calibrationFactors(3)) ); 
zeroCrossing = fzero(zeroFunc, guessZeroCrossValue)

% plot results
figure;
plot(scanPhShiftValues,meanFrascatiOutput,'.b');
hold all;
plot([zeroCrossing zeroCrossing], [minFrascatiOutput maxFrascatiOutput],'k','LineWidth',2);
fitPlotX = linspace(min(scanPhShiftValues),max(scanPhShiftValues),720);
fitPlotY = fitFunc(fitPlotX);
plot(fitPlotX, fitPlotY,'r','LineWidth',2);
xlabel('Phase shifter setting (degrees)');
ylabel('Frascati Output');
grid on;
legend('Data','Zero Crossing','Fit');
    
