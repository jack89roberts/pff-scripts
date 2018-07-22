error('Deprecated. Please use FrascatiCalibrations/runFrascatiCalibration.m')

%%
clear all;
close all;
%%
useMixerOverSqrtDiode = false; % if true calibrate using mixer/sqrt(diode), if false calibrate using mixer only
overwriteCalibrations = true; % if true, new zero crossings will be set and reference calibration files will be overwritten 
saveData = true; % if true all matlab data and files will be saved in the datadir directory

scanNPoints = 37; % no. of shifter settings in the scan
scanNTriggersPerPoint = 5; % no. of triggers collected for each shifter setting

datadir = 'data/frascatiCalibrations/';
dataSetComment = '';

phShiftNames = { ...
   'CT.PHAS-PHM0370',...
   'CC.PHAS-PHM0251',...
   'CB.PHAS-PHM0160'...
};


monMixerNames = {...
    'CT.SCOPE01.CH02',...
    'CT.SCOPE01.CH04',...
    'CT.SCOPE01.CH06'...
};

monDiodeNames = {...
    'CT.SCOPE01.CH01',...
    'CT.SCOPE01.CH03',...
    'CT.SCOPE01.CH05'...
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

phAcqPropStr = 'Acquisition';
phAcqValStr = 'value';

nMons = length(phShiftNames);

% get the initial values for all the phase shifters so we can set them back
% afterwards if necessary
fprintf(1,'Acquiring initial phase shifter settings...\n');
initPhShiftValues = NaN*ones(1,nMons);
for i=1:nMons
    initPhShiftValues(i) = JGetCoValueFesa(phShiftNames{i}, phShiftGetPropStr, phShiftGetValStr);
end

% initialise arrays
scanPhShiftValues = linspace(scanStartPhase,scanEndPhase,scanNPoints);
scanLOPhases = 4*scanPhShiftValues - min(4*scanPhShiftValues);

tmpPulse = JGetCoArrayFesa(monMixerNames{1},phAcqPropStr, phAcqValStr);
nSamples = length(tmpPulse);
scanMixers = NaN*ones(nMons,scanNPoints,scanNTriggersPerPoint, nSamples);
scanDiodes = NaN*ones(nMons,scanNPoints,scanNTriggersPerPoint, nSamples);

% units: Volts = acquisition*sensitivity
mixerSensitivities = NaN*ones(1,nMons);
diodeSensitivities = NaN*ones(1,nMons);
for mon=1:nMons
    mixerSensitivities(mon) = JGetCoValueFesa(monMixerNames{mon}, 'Acquisition', 'sensitivity');
    diodeSensitivities(mon) = JGetCoValueFesa(monDiodeNames{mon}, 'Acquisition', 'sensitivity');
end


%% start the scan

fprintf(1,'Starting scan...\n');

% if saving data, get the date now to make the file names
if (saveData) % save data (backup)   
    dV = datevec(now);
    saveFileName_data = sprintf('%sfrascatiCalibration_%02d%02d%02d_%02d%02d.mat',datadir,dV(1),dV(2),dV(3),dV(4),dV(5));
    saveFileName_zeroCross = sprintf('%sfrascatiZeroCrossings_%02d%02d%02d_%02d%02d',datadir,dV(1),dV(2),dV(3),dV(4),dV(5));
    saveFileName_constants = sprintf('%sfrascatiCalibrationConstants_%02d%02d%02d_%02d%02d',datadir,dV(1),dV(2),dV(3),dV(4),dV(5));      
end

nBeamDisabled = 0; % increment if had to stop acquiring data at some point.

for ph=1:scanNPoints % loop over phase shifter settings
    newPhVal = scanPhShiftValues(ph);
    fprintf(1,'Phase Shifter setting: %.2f degrees\n',newPhVal);
    
    for mon=1:nMons % change phase shifter values
        JSetCoValueFesa(phShiftNames{mon}, phShiftSetPropStr, phShiftSetValStr, newPhVal);
        JSetCoValueFesa(phShiftNames{mon}, phShiftSetPropStr, phShiftSetValStr, newPhVal);
        JSetCoValueFesa(phShiftNames{mon}, phShiftSetPropStr, phShiftSetValStr, newPhVal);
    end
    
    pause(5); % give time for phase shifters to change
    
    % get scanNTriggersPerPoint pulses and save them in arrays
    for t = 1:scanNTriggersPerPoint
        fprintf(1,'Acquiring trigger %d of %d\n',t, scanNTriggersPerPoint);
        
        nLoops = 0; % just used to avoid printing messages every iteration
        beamWasDisabled = false;
        
        while( ~isBeamPresent() ) % if the beam is not on, wait for it
            beamWasDisabled = true;
            if (mod(nLoops,10) == 0)
                fprintf(1,'Beam not enabled, waiting...\n');
            end
            
            nLoops = nLoops + 1;
            pause(10);
        end
        
        if (beamWasDisabled)
            nBeamDisabled = nBeamDisabled + 1;
        end
        
        for mon=1:nMons
            scanMixers(mon,ph,t,:) = JGetCoArrayFesa(monMixerNames{mon},phAcqPropStr, phAcqValStr);
            scanDiodes(mon,ph,t,:) = JGetCoArrayFesa(monDiodeNames{mon},phAcqPropStr, phAcqValStr);
        end
        
        pause(1.5); % wait some time for next pulse
    end
    
    if (saveData) % save data after each shifter setting so something is on disk if matlab crashes etc.
        save(saveFileName_data);

%         zeroCrossFileBackup = fopen(saveFileName_zeroCross,'w');
%         calibConstFileBackup = fopen(saveFileName_constants,'w');
%         for mon=1:nMons 
%             fprintf(zeroCrossFileBackup,'%s, %.2f\n', phShiftNames{mon}, zeroCrossings(mon));
%             fprintf(calibConstFileBackup,'%s, %f, %f, %f\n', monMixerNames{mon}, calibrationFactors(mon,1), calibrationFactors(mon,2), calibrationFactors(mon,3));
%         end
%         fclose(zeroCrossFileBackup);
%         fclose(calibConstFileBackup);
    end
    
end

% convert in to volts and calculate mixer/sqrt(diode)
scanMixerSqrtDiode = NaN*ones(nMons,scanNPoints,scanNTriggersPerPoint, nSamples);
for mon=1:nMons
    % convert to volts
    scanMixers(mon,:,:,:) = scanMixers(mon,:,:,:).*mixerSensitivities(mon);
    scanDiodes(mon,:,:,:) = scanDiodes(mon,:,:,:).*diodeSensitivities(mon);
    
    % subtract offset mean first 20 samples
    for ph=1:scanNPoints
        for t=1:scanNTriggersPerPoint
            scanMixers(mon,ph,t,:) = scanMixers(mon,ph,t,:) - mean(scanMixers(mon,ph,t,1:20));    
            scanDiodes(mon,ph,t,:) = scanDiodes(mon,ph,t,:) - mean(scanDiodes(mon,ph,t,1:20));
        end
    end
    
    % calculate mixer/sqrt(diode)
    scanMixerSqrtDiode(mon,:,:,:) = scanMixers(mon,:,:,:)./sqrt(abs(scanDiodes(mon,:,:,:)));
end

% deal with infs coming from diode=0 in mixer/sqrt(diode)
for mon=1:nMons
    for ph=1:scanNPoints
        for t=1:scanNTriggersPerPoint
            for s=1:nSamples
                if (scanMixerSqrtDiode(mon,ph,t,s) == Inf || scanMixerSqrtDiode(mon,ph,t,s) == -Inf);
                    scanMixerSqrtDiode(mon,ph,t,s) = NaN;
                end
            end
        end
    end
end

fprintf(1,'Scan finished.\n');

% reinstate initial phase shifter values
fprintf(1,'Reinstating initial phase shifter values...\n');
for mon=1:nMons
    JSetCoValueFesa(phShiftNames{mon}, phShiftSetPropStr, phShiftSetValStr, initPhShiftValues(mon));
end
%% plot a trace for each phase shifter setting in the scan and use it to prompt user to pick sample range to use.

calStartSamples = NaN*ones(1,nMons);
calEndSamples = NaN*ones(1,nMons);

figure;
for i=1:nMons    
    subplot(1,2,1);
    for n=1:scanNPoints
        plot(squeeze(scanMixers(i,n,1,:)));
        title(sprintf('Monitor %d Mixer',i));
        xlabel('Sample No.');
        ylabel('Output [V]');
        hold all;
    end
    hold off;
    subplot(1,2,2);
    for n=1:scanNPoints
        plot(squeeze(scanDiodes(i,n,1,:)));
        title(sprintf('Monitor %d Diode',i));
        xlabel('Sample No.');
        ylabel('Output [V]');
        hold all;
    end
    hold off;
    
    fprintf(1, 'Monitor %d, Shifter %s, Diode %s, Mixer%s:\n',mon,phShiftNames{i},monDiodeNames{i},monMixerNames{i});
    calStartSamples(i) = input('Start calibration at sample: ');
    calEndSamples(i) = input('End calibration at sample: ');
end

if (saveData) % save data before fitting
    save(saveFileName_data);

%     zeroCrossFileBackup = fopen(saveFileName_zeroCross,'w');
%     calibConstFileBackup = fopen(saveFileName_constants,'w');
%     for mon=1:nMons 
%         fprintf(zeroCrossFileBackup,'%s, %.2f\n', phShiftNames{mon}, zeroCrossings(mon));
%         fprintf(calibConstFileBackup,'%s, %f, %f, %f\n', monMixerNames{mon}, calibrationFactors(mon,1), calibrationFactors(mon,2), calibrationFactors(mon,3));
%     end
%     fclose(zeroCrossFileBackup);
%     fclose(calibConstFileBackup);
end

%% calculate calibration factors and find zero crossings

fprintf(1,'Calculating calibration factors and zero crossings...\n');

%% CALIBRATION AT ONE SAMPLE POINT
% fprintf(1, 'Calibration with one sample:\n');
% calSample = input('Calibrate using sample: ');
% 
% zeroCrossings = NaN*ones(nMons,1);
% calibrationFactors = NaN*ones(nMons,3);
% 
% meanFrascatiOutput = mean(scanMixers(:,:,:,calSample),3);
% stdFrascatiOutput = std(scanMixers(:,:,:,calSample),0,3); % NOTE - this not currently used for anything. Could/should use it for fitting weights.
% 
% % find max and min points in scan (max and min on sine curve). 
% [maxFrascatiOutput,maxFrascatiOutputIndex] = max(meanFrascatiOutput,[],2);
% [minFrascatiOutput,minFrascatiOutputIndex] = min(meanFrascatiOutput,[],2);
%     
% for mon=1:nMons
%     
%     % sensible guesses for fit parameters
%     initA = mean(abs([maxFrascatiOutput(mon) minFrascatiOutput(mon)]));
%     initB = 4*(pi./180); % 4x for 3GHz -> 12GHz
%     initC = (meanFrascatiOutput(mon,1)./initA)*(pi/2);
%     initParams = [initA, initB, initC];
%     
%     % calculate fit
%     [scanSinFit,scanSinGoF] = sinFit(scanPhShiftValues', meanFrascatiOutput(mon,:)', initParams);
%     
%     % save fit parameters
%     calibrationFactors(mon,1) = scanSinFit.a1;
%     calibrationFactors(mon,2) = scanSinFit.b1;
%     calibrationFactors(mon,3) = scanSinFit.c1;
%     fitFunc = @(x) (calibrationFactors(mon,1)).*sin( (calibrationFactors(mon,2)).*x + (calibrationFactors(mon,3)) );
%     
%     % find zero crossing of the fitted curve (guess half way between
%     % min/max)
%     guessZeroCrossIndex = abs(maxFrascatiOutputIndex-minFrascatiOutputIndex)./2; % index difference between min and max points
%     guessZeroCrossIndex = guessZeroCrossIndex + min([minFrascatiOutputIndex maxFrascatiOutputIndex]); % plus the starting index
%     guessZeroCrossIndex = round(guessZeroCrossIndex); % rounded to the nearest integer
%     guessZeroCrossValue = scanPhShiftValues(guessZeroCrossIndex); % phase shifter value at this index
%     zeroCrossings(mon) = fzero(fitFunc, guessZeroCrossValue);
%     
%     % plot results
%     figure;
%     plot(scanPhShiftValues,meanFrascatiOutput(mon,:),'.b');
%     hold all;
%     plot([zeroCrossings(mon) zeroCrossings(mon)], [minFrascatiOutput(mon) maxFrascatiOutput(mon)],'k','LineWidth',2);
%     fitPlotX = linspace(min(scanPhShiftValues),max(scanPhShiftValues),720);
%     fitPlotY = fitFunc(fitPlotX);
%     plot(fitPlotX, fitPlotY,'r','LineWidth',2);
%     title(monMixerNames{mon});
%     xlabel('Phase shifter setting (degrees)');
%     ylabel('Frascati Output');
%     grid on;
%     legend('Data','Zero Crossing','Fit');
%     
% end

%% CALIBRATION WITH MEAN PULSE PHASE
fprintf(1, 'Calibration with mean pulse phase:\n');

zeroCrossings = NaN*ones(nMons,1);
calibrationFactors = NaN*ones(nMons,4);

for mon=1:nMons
    
    % get mean output
    if (useMixerOverSqrtDiode == true)
        meanFrascatiOutput = nanmean(nanmean(scanMixerSqrtDiode(mon,:,:,calStartSamples(mon):calEndSamples(mon)),4),3);
        stdFrascatiOutput = nanstd(nanmean(scanMixerSqrtDiode(mon,:,:,calStartSamples(mon):calEndSamples(mon)),4),0,3); % NOTE - this not currently used for anything. Could/should use it for fitting weights.
    else
        meanFrascatiOutput = nanmean(nanmean(scanMixers(mon,:,:,calStartSamples(mon):calEndSamples(mon)),4),3);
        stdFrascatiOutput = nanstd(nanmean(scanMixers(mon,:,:,calStartSamples(mon):calEndSamples(mon)),4),0,3); % NOTE - this not currently used for anything. Could/should use it for fitting weights.
    end

    % find max and min points in scan (max and min on sine curve). 
    [maxFrascatiOutput,maxFrascatiOutputIndex] = max(meanFrascatiOutput);
    [minFrascatiOutput,minFrascatiOutputIndex] = min(meanFrascatiOutput);


    % find zero crossing of the fitted curve (guess half way between
    % min/max)
    guessZeroCrossIndex = abs(maxFrascatiOutputIndex-minFrascatiOutputIndex)./2; % index difference between min and max points
    guessZeroCrossIndex = guessZeroCrossIndex + min([minFrascatiOutputIndex maxFrascatiOutputIndex]); % plus the starting index
    guessZeroCrossIndex = round(guessZeroCrossIndex); % rounded to the nearest integer
    
    
    
    
    % sensible guesses for fit parameters
    initA = mean(abs([maxFrascatiOutput minFrascatiOutput])); % amplitude
    initB = 4*(pi./180); % phase units, 4x for 3GHz -> 12GHz
    initC = (meanFrascatiOutput(1)./initA)*(pi/2); % phase shift
    initD = (abs(maxFrascatiOutput) - abs(minFrascatiOutput))./2; % vertical offset
    initParams = [initA, initB, initC, initD];
    
    % calculate fit
    [scanSinFit,scanSinGoF] = offsetSinFit(scanPhShiftValues', meanFrascatiOutput', initParams);
    
    % save fit parameters
    calibrationFactors(mon,1) = scanSinFit.a;
    calibrationFactors(mon,2) = scanSinFit.b;
    calibrationFactors(mon,3) = scanSinFit.c;
    calibrationFactors(mon,4) = scanSinFit.d;
    fitFunc = @(x) (calibrationFactors(mon,1)).*sin( (calibrationFactors(mon,2)).*x + (calibrationFactors(mon,3)) ) + calibrationFactors(mon,4);
    
    % find zero crossing of the fitted curve (guess half way between
    % min/max)
    guessZeroCrossValue = scanPhShiftValues(guessZeroCrossIndex); % phase shifter value at this index
    % to get zero crossing use fit function without vertical offset,
    % otherwise zero crossing fit func ~= zero crossing sine curve
    % (actually looking for sin(x) = d, not sin(x) = 0)
    zeroFunc = @(x) (calibrationFactors(mon,1)).*sin( (calibrationFactors(mon,2)).*x + (calibrationFactors(mon,3)) ); 
    zeroCrossings(mon) = fzero(zeroFunc, guessZeroCrossValue);
    
    % plot results
    figure;
    plot(scanPhShiftValues,meanFrascatiOutput,'.b');
    hold all;
    plot([zeroCrossings(mon) zeroCrossings(mon)], [minFrascatiOutput maxFrascatiOutput],'k','LineWidth',2);
    fitPlotX = linspace(min(scanPhShiftValues),max(scanPhShiftValues),720);
    fitPlotY = fitFunc(fitPlotX);
    plot(fitPlotX, fitPlotY,'r','LineWidth',2);
    title(monMixerNames{mon});
    xlabel('Phase shifter setting (degrees)');
    ylabel('Frascati Output');
    grid on;
    legend('Data','Zero Crossing','Fit');
    
end

%% CALIBRATION FOR ALL SAMPLES
% fprintf(1, 'Calibration with mean pulse phase:\n');
% calStartSample = input('Start calibration at sample: ');
% calEndSample = input('End calibration at sample: ');
% 
% nCalSamples = length(calStartSample:calEndSample);
% 
% zeroCrossings = NaN*ones(nMons,nCalSamples);
% calibrationFactors = NaN*ones(nMons,nCalSamples,4);
% 
% for s=calStartSample:calEndSample
% 
%     if (useMixerOverSqrtDiode == true)
%         meanFrascatiOutput = mean(scanMixerSqrtDiode(:,:,:,s),3);
%         stdFrascatiOutput = std(scanMixerSqrtDiode(:,:,:,s),0,3); % NOTE - this not currently used for anything. Could/should use it for fitting weights.
%     else
%         meanFrascatiOutput = mean(scanMixers(:,:,:,s),3);
%         stdFrascatiOutput = std(scanMixers(:,:,:,s),0,3); % NOTE - this not currently used for anything. Could/should use it for fitting weights.
%     end
% 
%     % find max and min points in scan (max and min on sine curve). 
%     [maxFrascatiOutput,maxFrascatiOutputIndex] = max(meanFrascatiOutput,[],2);
%     [minFrascatiOutput,minFrascatiOutputIndex] = min(meanFrascatiOutput,[],2);
% 
%     for mon=1:nMons
% 
%         % sensible guesses for fit parameters
%         initA = mean(abs([maxFrascatiOutput(mon) minFrascatiOutput(mon)])); % amplitude
%         initB = 4*(pi./180); % phase units, 4x for 3GHz -> 12GHz
%         initC = (meanFrascatiOutput(mon,1)./initA)*(pi/2); % phase shift
%         initD = (abs(maxFrascatiOutput(mon)) - abs(minFrascatiOutput(mon)))./2; % vertical offset
%         initParams = [initA, initB, initC, initD];
% 
%         % calculate fit
%         [scanSinFit,scanSinGoF] = offsetSinFit(scanPhShiftValues', meanFrascatiOutput(mon,:)', initParams);
% 
%         % save fit parameters
%         calibrationFactors(mon,s,1) = scanSinFit.a;
%         calibrationFactors(mon,s,2) = scanSinFit.b;
%         calibrationFactors(mon,s,3) = scanSinFit.c;
%         calibrationFactors(mon,s,4) = scanSinFit.d;
%         
%         % find zero crossing of the fitted curve (guess half way between
%         % min/max)
%         guessZeroCrossIndex = abs(maxFrascatiOutputIndex-minFrascatiOutputIndex)./2; % index difference between min and max points
%         guessZeroCrossIndex = guessZeroCrossIndex + min([minFrascatiOutputIndex maxFrascatiOutputIndex]); % plus the starting index
%         guessZeroCrossIndex = round(guessZeroCrossIndex); % rounded to the nearest integer
%         guessZeroCrossValue = scanPhShiftValues(guessZeroCrossIndex); % phase shifter value at this index
%         % to get zero crossing use fit function without vertical offset,
%         % otherwise zero crossing fit func ~= zero crossing sine curve
%         % (actually looking for sin(x) = d, not sin(x) = 0)
%         zeroFunc = @(x) (calibrationFactors(mon,s,1)).*sin( (calibrationFactors(mon,s,2)).*x + (calibrationFactors(mon,s,3)) ); 
%         zeroCrossings(mon,s) = fzero(zeroFunc, guessZeroCrossValue);     
% 
%     end
% 
% end
% 
% % plot results
% monitorNames = {'Monitor 1','Monitor 2', 'Monitor 3'};
% calibFactorNames = {'a','b','c','d'};
% 
% for c=1:4 % calibration factor index
%     figure;
%     for mon=1:nMons
%         plot(squeeze(calibrationFactors(mon,:,c)));
%         hold all;
%     end
%     title(calibFactorNames{c});
%     legend(monitorNames);
%     xlabel('Sample No.');
%     ylabel('Parameter Value');
% end
% figure;
% for mon=1:nMons
%     plot(squeeze(zeroCrossings(mon,:)));
%     hold all;
% end
% title('Zero Crossing');
% legend(monitorNames);
% xlabel('Sample No.');
% ylabel('Phase [degrees]');

%% set zero crossings, save calibration factors

if (overwriteCalibrations) % set phase shifters to new zero crossings and overwrite old calibration files
    overwriteDefaults = input('Save this calibration as the new reference (1=yes/0=no)? ');
    zeroCrossFile = fopen('frascatiZeroCrossings','w');
    calibConstFile = fopen('frascatiCalibrationConstants','w');

    for mon=1:nMons 
        % set phase shifters to zero crossings, and save them to file
        
        % always want to set the zero crossing to the falling slope so check
        % this here
        dydxAtZero = calibrationFactors(mon,1)*calibrationFactors(mon,2)*cos(calibrationFactors(mon,2)*zeroCrossings(mon) + calibrationFactors(mon,3)); % Ab*cos(bx+c)
        if (dydxAtZero > 0)
            zeroCrossingToSet = zeroCrossings(mon) + (180/4); % divide by 4 because of 3 GHz phase shifter vs. 12 GHz mixer
        else
            zeroCrossingToSet = zeroCrossings(mon);
        end
        
        fprintf(1,'Setting zero crossing of %s to %.2f degrees...\n', phShiftNames{mon}, zeroCrossingToSet);
        JSetCoValueFesa(phShiftNames{mon}, phShiftSetPropStr, phShiftSetValStr, zeroCrossingToSet);
         
        if (overwriteDefaults)
            fprintf(1,'Overwriting calibration reference files...\n');
            fprintf(zeroCrossFile,'%s, %.2f\n', phShiftNames{mon}, zeroCrossings(mon));

            % save calibration constants to file
            fprintf(calibConstFile,'%s, %f, %f, %f, %f\n', monMixerNames{mon}, calibrationFactors(mon,1), calibrationFactors(mon,2), calibrationFactors(mon,3), calibrationFactors(mon,4));
        end
    end

    fclose(zeroCrossFile);
    fclose(calibConstFile);
end

if (saveData) % save data (backup)
    
    save(saveFileName_data);
    
    zeroCrossFileBackup = fopen(saveFileName_zeroCross,'w');
    calibConstFileBackup = fopen(saveFileName_constants,'w');
    for mon=1:nMons 
        fprintf(zeroCrossFileBackup,'%s, %.2f\n', phShiftNames{mon}, zeroCrossings(mon));
        fprintf(calibConstFileBackup,'%s, %f, %f, %f, %f\n', monMixerNames{mon}, calibrationFactors(mon,1), calibrationFactors(mon,2), calibrationFactors(mon,3), calibrationFactors(mon,4));
    end
    fclose(zeroCrossFileBackup);
    fclose(calibConstFileBackup);
    
end

calibrationFactors