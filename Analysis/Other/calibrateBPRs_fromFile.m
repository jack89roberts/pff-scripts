close all; clearvars;
%%
calStamp = '20151123_1644';
calDir = '/home/jack/PhaseFeedforward/CTFData/201511/BPRCalibrations';
saveNewData=1;
%%
calFile = [calDir '/bprCalibration_' calStamp '.mat'];
load(calFile);

%% calculate calibration factors and find zero crossings

% need to find calibration factors and which slope of the sine curve to use
% to avoid non-linear part of phase shifter, as well as zero crossings to
% give best resolution.

fprintf(1,'Calculating calibration factors and zero crossings...\n');

zeroCrossings = NaN*ones(nBPRs,1);
calibrationFactors = NaN*ones(nBPRs,4);


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
%     initA = mean(abs([maxBPROutput(bpr) minBPROutput(bpr)]));
%     initB = pi./180;
%     initC = (meanBPROutput(bpr,1)./initA)*(pi/2);
%     initParams = [initA, initB, initC];   
    initA = (maxBPROutput(bpr)-minBPROutput(bpr))/2; % amplitude
    initB = pi./180; % phase units, 4x for 3GHz -> 12GHz
    initD = (maxBPROutput(bpr)+minBPROutput(bpr))/2; % vertical offset
    % calculate phase offset using other guesses for initial
    % parameters and phase/output at first point.
    sinArg = (meanBPROutput(bpr,:)-initD)./initA;
    sinArg(sinArg>1)= NaN; % avoid imaginary component if data point above amplitude guess
    sinArg(sinArg<-1)= NaN;
    sinArg = asin(sinArg);
    % need to select 1st quadrant where asin valid - positive
    % gradient in difference
    sinArg(diff(meanBPROutput(bpr,:))<0) = NaN;
    initC = sinArg-(initB.*scanPhShiftValues);
    initC = mod(initC,2*pi); % make everything between 0 and 2pi
    initC(initC>pi) = initC(initC>pi)-2*pi;
    initC(initC<-pi) = initC(initC<-pi)+2*pi;
    initC = nanmean(initC);
    initParams = [initA, initB, initC, initD];
  
    
%     [slope1Fit,slope1GoF] = sinFit(slope1Phases', meanBPROutputSlope1', initParams);
%     [slope2Fit,slope2GoF] = sinFit(slope2Phases', meanBPROutputSlope2', initParams);

    [slope1Fit,slope1RSquare,slope1Conf] = offsetSinFit(slope1Phases', meanBPROutputSlope1', initParams);
    [slope2Fit,slope2RSquare,slope2Conf] = offsetSinFit(slope2Phases', meanBPROutputSlope2', initParams);
    
    if (slope1RSquare > slope2RSquare)
        slopes(bpr) = 1;
        
        calibrationFactors(bpr,1) = slope1Fit(1);
        calibrationFactors(bpr,2) = slope1Fit(2);
        calibrationFactors(bpr,3) = slope1Fit(3);
        calibrationFactors(bpr,4) = slope1Fit(4);

        fitFunc = @(x) slope1Fit(1).*sin( slope1Fit(2).*x + slope1Fit(3) ) + slope1Fit(4);
        
        % find zero crossing of the fitted curve (guess half way between
        % min/max)
        guessZeroCross = scanPhShiftValues(round(((slope1End(bpr)-slope1Start(bpr))./2)+slope1Start(bpr)));
        zeroCrossings(bpr) = fzero(fitFunc, guessZeroCross);
        
    else
        slopes(bpr) = 2;
        
        calibrationFactors(bpr,1) = slope2Fit(1);
        calibrationFactors(bpr,2) = slope2Fit(2);
        calibrationFactors(bpr,3) = slope2Fit(3);
        calibrationFactors(bpr,4) = slope2Fit(4);

        fitFunc = @(x) slope2Fit(1).*sin( slope2Fit(2).*x + slope2Fit(3) ) + slope2Fit(4);
        
        % find zero crossing of the fitted curve (guess half way between
        % min/max)        
        guessZeroCross = scanPhShiftValues(round(((slope1End(bpr)-slope1Start(bpr))./2)+slope1Start(bpr)))+180;
        zeroCrossings(bpr) = fzero(fitFunc, guessZeroCross);        
    end
    
    if (zeroCrossings(bpr) > 360) % deal with zero crossing outside of range
        zeroCrossings(bpr) = zeroCrossings(bpr) - 360;
    end
    
     % plot results
    figure;
    plot(scanPhShiftValues,scanBPRMeansAtCursor(bpr,:,1),'.b');
    fitFunc = @(x) (calibrationFactors(bpr,1)).*sin( (calibrationFactors(bpr,2)).*x + (calibrationFactors(bpr,3)) ) + calibrationFactors(bpr,4);
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

%% save results
if (saveNewData) % save data (backup)
    
    zeroCrossFileBackup = fopen(sprintf('%s/bprZeroCrossings_%s',calDir,calStamp),'w');
    calibConstFileBackup = fopen(sprintf('%s/bprCalibrationConstants_%s',calDir,calStamp),'w');
    
    for bpr=1:nBPRs 
        fprintf(zeroCrossFileBackup,'%s, %.2f\n', phShiftNames{bpr}, zeroCrossings(bpr));
        fprintf(calibConstFileBackup,'%s, %f, %f, %f, %f\n', bprNames{bpr}, calibrationFactors(bpr,1), calibrationFactors(bpr,2), calibrationFactors(bpr,3),calibrationFactors(bpr,4));
    end
    
    fclose(zeroCrossFileBackup);
    fclose(calibConstFileBackup);
    
end