% NEED TO CHECK CURSORS CORRECT FIRST!!!

overwriteCalibrations = true; % if true, new zero crossings will be set and reference calibration files will be overwritten 
saveData = true; % if true all matlab data and files will be saved in the datadir directory

scanNPoints = 361; % no. of shifter settings in the scan
scanNTriggersPerPoint = 1; % no. of triggers collected for each shifter setting

datadir = 'data/frascatiCalibrations/';

phShiftNames = { ...
   'CC.PHAS-PHM0251'...
};


frascatiNames = {...
    'CT.SCOPE01.CH02'...
};

%%

% CHANGING THIS RANGE WILL PROBABLY CAUSE THE CALIBRATION TO FAIL
% PHASE SHIFTERS AT 3GHz, MULTIPLIED TO 12GHz FOR PHASE MONITORS
% 90 degree range shifters = 360 degree LO scan
scanStartPhase = 90.0; % initial value for phase shifters in the scan
scanEndPhase = 180.0; % final value for phase shifters in the scan

%%

addpath('../JavaCoInterface/');

phShiftGetPropStr = 'Acquisition';
phShiftGetValStr = 'phase_Val';

phShiftSetPropStr = 'Setting';
phShiftSetValStr = 'phase';

frascatiGetPropStr = 'Acquisition';
frascatiGetValStr = 'value';



nFrascatiMons = length(phShiftNames);

% get the initial values for all the phase shifters so we can set them back
% afterwards if necessary
fprintf(1,'Acquiring initial phase shifter settings...\n');
initPhShiftValues = NaN*ones(1,nFrascatiMons);
for i=1:nFrascatiMons
    initPhShiftValues(i) = JGetCoValueFesa(phShiftNames{i}, phShiftGetPropStr, phShiftGetValStr);
end

scanPhShiftValues = linspace(scanStartPhase,scanEndPhase,scanNPoints);
scanLOPhases = 4*scanPhShiftValues - min(4*scanPhShiftValues);

tmpPulse = JGetCoArrayFesa(frascatiNames{1},frascatiGetPropStr, frascatiGetValStr);
nFrascatiSamples = length(tmpPulse);
scanFrascatiPulses = NaN*ones(nFrascatiMons,scanNPoints,scanNTriggersPerPoint, nFrascatiSamples);

% get units: Volts = acquisition*sensitivity + offset ?
frascatiSensitivity = JGetCoValueFesa(frascatiNames{1}, 'Acquisition', 'sensitivity');
frascatiOffset = JGetCoValueFesa(frascatiNames{1},'Acquisition','offset');

%% start the scan

fprintf(1,'Starting scan...\n');
for ph=1:scanNPoints % loop over phase shifter settings
    newPhVal = scanPhShiftValues(ph);
    fprintf(1,'Phase Shifter setting: %.2f degrees\n',newPhVal);
    
    for mon=1:nFrascatiMons % change phase shifter values
        JSetCoValueFesa(phShiftNames{mon}, phShiftSetPropStr, phShiftSetValStr, newPhVal);
        JSetCoValueFesa(phShiftNames{mon}, phShiftSetPropStr, phShiftSetValStr, newPhVal);
        JSetCoValueFesa(phShiftNames{mon}, phShiftSetPropStr, phShiftSetValStr, newPhVal);
    end
    
    pause(5); % give time for phase shifters to change
    
    % get scanNTriggersPerPoint pulses and save them in
    % scanBPRMeansAtCursor
    for t = 1:scanNTriggersPerPoint
        fprintf(1,'Acquiring trigger %d of %d\n',t, scanNTriggersPerPoint);
        
        nLoops = 0; % just used to avoid printing messages every iteration
        while( ~isBeamPresent() ) % if the beam is not on, wait for it
            if (mod(nLoops,10) == 0)
                fprintf(1,'Beam not enabled, waiting...\n');
            end
            
            nLoops = nLoops + 1;
            pause(10);
        end
        
        for mon=1:nFrascatiMons
            scanFrascatiPulses(mon,ph,t,:) = JGetCoArrayFesa(frascatiNames{mon},frascatiGetPropStr, frascatiGetValStr);
        end
        pause(1.5); % wait some time for next pulse
    end
    
end
fprintf(1,'Scan finished.\n');

% reinstate initial phase shifter values
fprintf(1,'Reinstating initial phase shifter values...\n');
for mon=1:nFrascatiMons
    JSetCoValueFesa(phShiftNames{mon}, phShiftSetPropStr, phShiftSetValStr, initPhShiftValues(mon));
end

% quickly plot a trace for each phase shifter setting in the scan
figure;
for n=1:15:scanNPoints
    plot(squeeze(scanFrascatiPulses(1,n,1,:)));
    hold all;
end

%% calculate calibration factors and find zero crossings

fprintf(1,'Calculating calibration factors and zero crossings...\n');

zeroCrossings = NaN*ones(nFrascatiMons,1);
calibrationFactors = NaN*ones(nFrascatiMons,3);

% just use means in some range for now
calStartSample = input('Start calibration at sample: ');
calEndSample = input('End calibration at sample: ');

meanFrascatiOutput = mean(mean(scanFrascatiPulses(:,:,:,calStartSample:calEndSample),4),3);
stdFrascatiOutput = std(mean(scanFrascatiPulses(:,:,:,calStartSample:calEndSample),4),0,3); % NOTE - this not currently used for anything. Could/should use it for fitting weights.

% find max and min points in scan (max and min on sine curve). 
[maxFrascatiOutput,maxFrascatiOutputIndex] = max(meanFrascatiOutput,[],2);
[minFrascatiOutput,minFrascatiOutputIndex] = min(meanFrascatiOutput,[],2);
    
for mon=1:nFrascatiMons
    
    % sensible guesses for fit parameters
    initA = mean(abs([maxFrascatiOutput(mon) minFrascatiOutput(mon)]));
    initB = 4*(pi./180); % 4x for 3GHz -> 12GHz
    initC = (meanFrascatiOutput(mon,1)./initA)*(pi/2);
    initParams = [initA, initB, initC];
    
    % calculate fit
    [scanSinFit,scanSinGoF] = sinFit(scanPhShiftValues', meanFrascatiOutput', initParams);
    
    % save fit parameters
    calibrationFactors(mon,1) = scanSinFit.a1;
    calibrationFactors(mon,2) = scanSinFit.b1;
    calibrationFactors(mon,3) = scanSinFit.c1;
    fitFunc = @(x) (calibrationFactors(mon,1)).*sin( (calibrationFactors(mon,2)).*x + (calibrationFactors(mon,3)) );
    
    % find zero crossing of the fitted curve (guess half way between
    % min/max)
    guessZeroCrossIndex = abs(maxFrascatiOutputIndex-minFrascatiOutputIndex)./2; % index difference between min and max points
    guessZeroCrossIndex = guessZeroCrossIndex + min([minFrascatiOutputIndex maxFrascatiOutputIndex]); % plus the starting index
    guessZeroCrossIndex = round(guessZeroCrossIndex); % rounded to the nearest integer
    guessZeroCrossValue = scanPhShiftValues(guessZeroCrossIndex); % phase shifter value at this index
    zeroCrossings(mon) = fzero(fitFunc, guessZeroCrossValue);
    
    % plot results
    figure;
    plot(scanPhShiftValues,meanFrascatiOutput(mon,:),'.b');
    hold all;
    plot([zeroCrossings(mon) zeroCrossings(mon)], [minFrascatiOutput(mon) maxFrascatiOutput(mon)],'k','LineWidth',2);
    fitPlotX = linspace(min(scanPhShiftValues),max(scanPhShiftValues),720);
    fitPlotY = fitFunc(fitPlotX);
    plot(fitPlotX, fitPlotY,'r','LineWidth',2);
    title(frascatiNames{mon});
    xlabel('Phase shifter setting (degrees)');
    ylabel('Frascati Output');
    grid on;
    legend('Data','Zero Crossing','Fit');
    
end

%% set zero crossings, save calibration factors

if (overwriteCalibrations) % set phase shifters to new zero crossings and overwrite old calibration files
    zeroCrossFile = fopen('frascatiZeroCrossings','w');
    calibConstFile = fopen('frascatiCalibrationConstants','w');

    for mon=1:nFrascatiMons 
        % set phase shifters to zero crossings, and save them to file
        fprintf(1,'Setting zero crossing of %s to %.2f degrees...\n', phShiftNames{mon}, zeroCrossings(mon));
        JSetCoValueFesa(phShiftNames{mon}, phShiftSetPropStr, phShiftSetValStr, zeroCrossings(mon));

        fprintf(zeroCrossFile,'%s, %.2f\n', phShiftNames{mon}, zeroCrossings(mon));

        % save calibration constants to file
        fprintf(calibConstFile,'%s, %f, %f, %f\n', frascatiNames{mon}, calibrationFactors(mon,1), calibrationFactors(mon,2), calibrationFactors(mon,3));

    end

    fclose(zeroCrossFile);
    fclose(calibConstFile);
end

if (saveData) % save data (backup)
    
    dV = datevec(now);
    save(sprintf('%sfrascatiCalibration_%02d%02d%02d_%02d%02d.mat',datadir,dV(1),dV(2),dV(3),dV(4),dV(5)));
    
    zeroCrossFileBackup = fopen(sprintf('%sfrascatiZeroCrossings_%02d%02d%02d_%02d%02d',datadir,dV(1),dV(2),dV(3),dV(4),dV(5)),'w');
    calibConstFileBackup = fopen(sprintf('%sfrascatiCalibrationConstants_%02d%02d%02d_%02d%02d',datadir,dV(1),dV(2),dV(3),dV(4),dV(5)),'w');
    
    for mon=1:nFrascatiMons 
        fprintf(zeroCrossFileBackup,'%s, %.2f\n', phShiftNames{mon}, zeroCrossings(mon));
        fprintf(calibConstFileBackup,'%s, %f, %f, %f\n', frascatiNames{mon}, calibrationFactors(mon,1), calibrationFactors(mon,2), calibrationFactors(mon,3));
    end
    
    fclose(zeroCrossFileBackup);
    fclose(calibConstFileBackup);
end