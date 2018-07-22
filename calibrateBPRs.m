% NEED TO CHECK CURSORS CORRECT FIRST!!!

overwriteCalibrations = true; % if true, new zero crossings will be set and reference calibration files will be overwritten 
saveData = true; % if true all matlab data and files will be saved in the datadir directory

scanNPoints = 17; % no. of shifter settings in the scan
scanNTriggersPerPoint = 5; % no. of triggers collected for each shifter setting

datadir = '/user/ctf3op/PhaseFeedforward/BPRCalibrations/';

phShiftNames = { ...
   % 'CL.PHAS-BPR0290',...
   %'CL.PHAS-BPR0475',...
   'CT.PHAS-BPR0532'...
%    'CR.PHAS-BPR0505',...
%    'CR.PHAS-BPR0915',...
};


bprNames = {...
    %'CL.STBPR0290S',...
    %'CL.STBPR0475S',...
    'CT.STBPR0532S'...
%     'CR.STBPR0505S',...
%     'CC.STBPR0915S',...
};

%%

% CHANGING THIS RANGE WILL PROBABLY CAUSE THE CALIBRATION TO FAIL
scanStartPhase = 0.0; % initial value for phase shifters in the scan
scanEndPhase = 360.0; % final value for phase shifters in the scan

%%

addpath('../JavaCoInterface/');

phShiftGetPropStr = 'Acquisition';
phShiftGetValStr = 'phase_Val';

phShiftSetPropStr = 'Setting';
phShiftSetValStr = 'phase';

bprGetPropStr = 'MeanAtCursor';
bprGetValStr = 'mean';

cycleName = 'SCT.USER.SETUP';

nBPRs = length(phShiftNames);

% get the initial values for all the phase shifters so we can set them back
% afterwards if necessary
fprintf(1,'Acquiring initial phase shifter settings...\n');
initPhShiftValues = NaN*ones(1,nBPRs);
for i=1:nBPRs
    initPhShiftValues(i) = ...
        matlabJapc.staticGetSignal(cycleName,...
            sprintf('%s/%s#%s',phShiftNames{i},phShiftGetPropStr,phShiftGetValStr)...
        );
                
        %JGetCoValueFesa(phShiftNames{i}, phShiftGetPropStr, phShiftGetValStr);
end

scanPhShiftValues = linspace(scanStartPhase,scanEndPhase,scanNPoints);
scanBPRMeansAtCursor = NaN*ones(nBPRs,scanNPoints,scanNTriggersPerPoint);

%% start the scan

fprintf(1,'Starting scan...\n');
for ph=1:scanNPoints % loop over phase shifter settings
    newPhVal = scanPhShiftValues(ph);
    fprintf(1,'Phase Shifter setting: %d degrees\n',newPhVal);
    
    for bpr=1:nBPRs % change phase shifter values
        %JSetCoValueFesa(phShiftNames{bpr}, phShiftSetPropStr, phShiftSetValStr, newPhVal);
        
        matlabJapc.staticSetSignal(cycleName,...
                        sprintf('%s/%s#%s',phShiftNames{bpr},phShiftSetPropStr,phShiftSetValStr),...
                        newPhVal...
                    );
    end
    
    pause(6); % give time for phase shifters to change
    
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
            pause(5);
        end
        
        for bpr=1:nBPRs
            %scanBPRMeansAtCursor(bpr,ph,t) = JGetCoValueFesa(bprNames{bpr},bprGetPropStr, bprGetValStr);
            
            scanBPRMeansAtCursor(bpr,ph,t) = ...
                matlabJapc.staticGetSignal(cycleName,...
                    sprintf('%s/%s#%s',bprNames{bpr},bprGetPropStr,bprGetValStr)...
                );
            
        end
        pause(1.5); % wait some time for next pulse
    end
    
end
fprintf(1,'Scan finished.\n');

% reinstate initial phase shifter values
fprintf(1,'Reinstating initial phase shifter values...\n');
for bpr=1:nBPRs
    %JSetCoValueFesa(phShiftNames{bpr}, phShiftSetPropStr, phShiftSetValStr, initPhShiftValues(bpr));
    
    matlabJapc.staticSetSignal(cycleName,...
                        sprintf('%s/%s#%s',phShiftNames{bpr},phShiftSetPropStr,phShiftSetValStr),...
                        initPhShiftValues(bpr)...
                    );
end

%% calculate calibration factors and find zero crossings

% need to find calibration factors and which slope of the sine curve to use
% to avoid non-linear part of phase shifter, as well as zero crossings to
% give best resolution.

fprintf(1,'Calculating calibration factors and zero crossings...\n');

zeroCrossings = NaN*ones(nBPRs,1);
calibrationFactors = NaN*ones(nBPRs,3);


meanBPROutput = mean(scanBPRMeansAtCursor,3);
stdBPROutput = std(scanBPRMeansAtCursor,0,3); % NOTE - this not currently used for anything. Could/should use it for fitting weights.

% find max and min points in scan (max and min on sine curve). Then use
% this to define location of two slopes of sine curve.
[maxBPROutput,maxBPROutputIndex] = max(meanBPROutput,[],2);
[minBPROutput,minBPROutputIndex] = min(meanBPROutput,[],2);

% first slope between min/max points
slope1Start = min([maxBPROutputIndex minBPROutputIndex],[],2);
slope1End = max([maxBPROutputIndex minBPROutputIndex],[],2);

slopes = NaN*ones(1,nBPRs); % which slope fit has decided to use for each bpm
    
for bpr=1:nBPRs
    
    slope1Range = slope1Start(bpr):slope1End(bpr); % index list for slope 1
    slope1Phases = scanPhShiftValues(slope1Range); % phase shifter values for slope 1
    meanBPROutputSlope1 = meanBPROutput(bpr,slope1Range);
    stdBPROutputSlope1 = stdBPROutput(bpr,slope1Range);
    
    
    % second slope: from 1st point in scan to first min/max, then from second max/min to
    % last point in scan.
    meanBPROutputSlope2 = [ meanBPROutput(bpr,slope1End(bpr):scanNPoints) meanBPROutput(bpr,1:slope1Start(bpr)) ];
    stdBPROutputSlope2 = [ stdBPROutput(bpr,slope1End(bpr):scanNPoints) stdBPROutput(bpr,1:slope1Start(bpr)) ];
    
    slope2Range = [slope1End(bpr):scanNPoints 1:slope1Start(bpr)]; % index list for slope 2
    
    slope2Phases = NaN*zeros(1,length(slope2Range)); % phase shifter values for slope 2
    tmpAddPhase = 0; % used to wrap around 360 degrees
    for i=1:length(slope2Range)
        if (slope2Range(i) == 1) % phShifter = 0 degrees
            tmpAddPhase = 360; % add 360 degrees on to values in first half of scan to get a continuous list of phase values 
        end
        slope2Phases(i) = scanPhShiftValues(slope2Range(i)) + tmpAddPhase;
    end

    % sine fit on slope 1 and slope 2. One with worst rsquared value should
    % be the one with the non-linear effects from the phase shifter.
    
    % sensible guesses for fit parameters
    initA = mean(abs([maxBPROutput(bpr) minBPROutput(bpr)]));
    initB = pi./180;
    initC = (meanBPROutput(bpr,1)./initA)*(pi/2);
    initParams = [initA, initB, initC];
    
    [slope1Fit,slope1GoF] = sinFit(slope1Phases', meanBPROutputSlope1', initParams);
    [slope2Fit,slope2GoF] = sinFit(slope2Phases', meanBPROutputSlope2', initParams);
    
    if (slope1GoF.rsquare > slope2GoF.rsquare)
        slopes(bpr) = 1;
        
        calibrationFactors(bpr,1) = slope1Fit.a1;
        calibrationFactors(bpr,2) = slope1Fit.b1;
        calibrationFactors(bpr,3) = slope1Fit.c1;
        
        fitFunc = @(x) (slope1Fit.a1).*sin( (slope1Fit.b1).*x + (slope1Fit.c1) );
        
        % find zero crossing of the fitted curve (guess half way between
        % min/max)
        guessZeroCross = scanPhShiftValues(round(((slope1End(bpr)-slope1Start(bpr))./2)+slope1Start(bpr)));
        zeroCrossings(bpr) = fzero(fitFunc, guessZeroCross);
        
    else
        slopes(bpr) = 2;
        
        calibrationFactors(bpr,1) = slope2Fit.a1;
        calibrationFactors(bpr,2) = slope2Fit.b1;
        calibrationFactors(bpr,3) = slope2Fit.c1;
        
        fitFunc = @(x) (slope2Fit.a1).*sin( (slope2Fit.b1).*x + (slope2Fit.c1) );
        
        % find zero crossing of the fitted curve (guess half way between
        % min/max + 180 degrees)
        guessZeroCross = scanPhShiftValues(round(((slope1End(bpr)-slope1Start(bpr))./2)+slope1Start(bpr)))+180;
        zeroCrossings(bpr) = fzero(fitFunc, guessZeroCross);
        if (zeroCrossings(bpr) > 360) % deal with zero crossing outside of range
            zeroCrossings(bpr) = zeroCrossings(bpr) - 360;
        end
    end
    
     % plot results
    figure;
    plot(scanPhShiftValues,scanBPRMeansAtCursor(bpr,:,1),'.b');
    fitFunc = @(x) (calibrationFactors(bpr,1)).*sin( (calibrationFactors(bpr,2)).*x + (calibrationFactors(bpr,3)) );
    hold all;
    plot([zeroCrossings(bpr) zeroCrossings(bpr)], [minBPROutput(bpr) maxBPROutput(bpr)],'k','LineWidth',2);
    if (slopes(bpr) == 1)
        plot(slope1Phases,fitFunc(slope1Phases),'r','LineWidth',3);
    else
        fitValues = fitFunc(slope2Phases);
        firstRange = 1:length(slope1End(bpr):scanNPoints);
        secondRange = (length(slope1End(bpr):scanNPoints)+1):length(slope2Phases);
        plot(slope2Phases(firstRange),fitValues(firstRange),'r','LineWidth',3);
        plot(slope2Phases(secondRange)-360,fitValues(secondRange),'r','LineWidth',3);
    end
    title(bprNames{bpr});
    xlabel('Phase shifter setting (degrees)');
    ylabel('BPR Output (A)');
    grid on;
    legend('Data','Zero Crossing','Fit');
    
end

%% set zero crossings, save calibration factors

if (overwriteCalibrations) % set phase shifters to new zero crossings and overwrite old calibration files
    zeroCrossFile = fopen('bprZeroCrossings','w');
    calibConstFile = fopen('bprCalibrationConstants','w');

    for bpr=1:nBPRs 
        % set phase shifters to zero crossings, and save them to file
        fprintf(1,'Setting zero crossing of %s to %.2f degrees...\n', phShiftNames{bpr}, zeroCrossings(bpr));
        %JSetCoValueFesa(phShiftNames{bpr}, phShiftSetPropStr, phShiftSetValStr, zeroCrossings(bpr));

        matlabJapc.staticSetSignal(cycleName,...
                        sprintf('%s/%s#%s',phShiftNames{bpr},phShiftSetPropStr,phShiftSetValStr),...
                        zeroCrossings(bpr)...
                    );
        
        fprintf(zeroCrossFile,'%s, %.2f\n', phShiftNames{bpr}, zeroCrossings(bpr));

        % save calibration constants to file
        fprintf(calibConstFile,'%s, %f, %f, %f\n', bprNames{bpr}, calibrationFactors(bpr,1), calibrationFactors(bpr,2), calibrationFactors(bpr,3));

    end

    fclose(zeroCrossFile);
    fclose(calibConstFile);
end

if (saveData) % save data (backup)
    
    dV = datevec(now);
    save(sprintf('%sbprCalibration_%02d%02d%02d_%02d%02d.mat',datadir,dV(1),dV(2),dV(3),dV(4),dV(5)));
    
    zeroCrossFileBackup = fopen(sprintf('%sbprZeroCrossings_%02d%02d%02d_%02d%02d',datadir,dV(1),dV(2),dV(3),dV(4),dV(5)),'w');
    calibConstFileBackup = fopen(sprintf('%sbprCalibrationConstants_%02d%02d%02d_%02d%02d',datadir,dV(1),dV(2),dV(3),dV(4),dV(5)),'w');
    
    for bpr=1:nBPRs 
        fprintf(zeroCrossFileBackup,'%s, %.2f\n', phShiftNames{bpr}, zeroCrossings(bpr));
        fprintf(calibConstFileBackup,'%s, %f, %f, %f\n', bprNames{bpr}, calibrationFactors(bpr,1), calibrationFactors(bpr,2), calibrationFactors(bpr,3));
    end
    
    fclose(zeroCrossFileBackup);
    fclose(calibConstFileBackup);
end