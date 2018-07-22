clear all
close all
%%
% commented out all savePngEpsFig for plot saving
% NOTE - PARTS WILL NEED TO BE CHANGED TO ALLOW FOR DIFFERENT SAMPLE RANGE
% FOR 3RD MONITOR (only mean sine fit changed)

dataDir = '/home/jack/PhaseFeedforward/CTFData/201412/frascatiCalibrations';
dataFile = 'frascatiCalibration_20141210_1844'; % WITHOUT .MAT EXTENSION
% load file
load([dataDir '/' dataFile '.mat'], '-regexp', '^(?!dataFile$)\w'); % load file, but don't load dataDir variable (would overwrite dataDir value above)

useMixerOverSqrtDiode = true; % if true calibrate using mixer/sqrt(diode), if false calibrate using mixer only
overwriteCalibrations = false; % if true, new zero crossings will be set and reference calibration files will be overwritten 
selectCalSamples = false; % if true, will ask for new sample ranges to use (and make plots)
doLinearCalibrations = false; % if true will perform linear calibrations
doPerSampleCalibrations = false; % if true will perform calibrations point by point for each sample

removeBadPulses = true;
removeThreshold = 0.01; % remove pulses more than removeThreshold V below mean transmission (based on diode channel).

linearRange = 60; % linear fit in range +/- range degrees (of 12 GHz);

monitorNames = {'Mon1','Mon2', 'Mon3'};

%%

saveDir = [dataDir '/plots_' dataFile];

[nMons,scanNPoints,nAvg,nSamples] = size(scanDiodes);


if (~exist('calStartSamples','var')) % deal with old data where only one sample range used for all monitors
    calStartSamples = NaN*ones(1,nMons);
    calEndSamples = NaN*ones(1,nMons);
    for mon=1:nMons
        calStartSamples(mon) = calStartSample; 
        calEndSamples(mon) = calEndSample;
    end
end



if (selectCalSamples)
% plot a trace for each phase shifter setting in the scan and use it to prompt user to pick sample range to use.
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
end


if (removeBadPulses) % remove pulses outside transmission window (based on diode signal amplitude)
    transmission = NaN*ones(nMons,scanNPoints,nAvg);
    isGoodPulse = NaN*ones(nMons,scanNPoints,nAvg);
    nBadPulses = NaN*ones(1,nMons);
    for mon=1:nMons
        monTransmission = transmission(mon,:,:);
        meanMonTransmission = nanmean(monTransmission(:));
        for p=1:scanNPoints
            for n=1:nAvg
                transmission(mon,p,n) = nanmean(squeeze(scanDiodes(mon,p,n,calStartSamples(mon):calEndSamples(mon))));
                isGoodPulse(mon,p,n) = ~((transmission(mon,p,n)-meanMonTransmission) > removeThreshold);
                if (~isGoodPulse(mon,p,n))
                    scanDiodes(mon,p,n,:) = NaN;
                    scanMixers(mon,p,n,:) = NaN;
                    scanMixerSqrtDiode(mon,p,n,:) = NaN;
                end
            end
        end
        
        tmpBadPulses = ~isGoodPulse(mon,:,:);
        nBadPulses(mon) = sum(tmpBadPulses(:));
    end
    
    nBadPulses
end


% figure;
% plot(max(scanMixers));
% hold all;
% plot(mean(scanMixers));
% plot(min(scanMixers));
% xlabel('Sample No.');
% ylabel('Output [V]');
% title('Mixer: Max, Min, Mean Output Per Sample');
% legend('Max','Mean','Min');
% %savePngEpsFig(saveDir,'Mixer_MaxMeanMin',1);
% 
% figure;
% plot(max(scanDiodes));
% hold all;
% plot(mean(scanDiodes));
% plot(min(scanDiodes));
% xlabel('Sample No.');
% ylabel('Output [V]');
% title('Diode: Max, Min, Mean Output Per Sample');
% legend('Max','Mean','Min');
% %savePngEpsFig(saveDir,'Diode_MaxMeanMin',1);


fprintf(1,'Calculating calibration factors and zero crossings...\n');


for i=1:nMons
    figure;
    for n=1:scanNPoints
        plot(squeeze(scanMixers(i,n,1,:)));
        title(sprintf('%s Mixer',monitorNames{i}));
        xlabel('Sample No.');
        ylabel('Output [V]');
        hold all;
    end
    %savePngEpsFig(saveDir,sprintf('Mixers_%s',monitorNames{mon}),1);
    figure;
    for n=1:scanNPoints
        plot(squeeze(scanDiodes(i,n,1,:)));
        title(sprintf('%s Diode',monitorNames{i}));
        xlabel('Sample No.');
        ylabel('Output [V]');
        hold all;
    end
    %savePngEpsFig(saveDir,sprintf('Diodes_%s',monitorNames{mon}),1);
end

%nCalSamples = length(calStartSample:calEndSample);



%% CALIBRATION WITH MEAN PULSE PHASE
% fprintf(1, 'Calibration with mean pulse phase:\n');
% 
% meanZeroCrossings = NaN*ones(nMons,1);
% meanCalibrationFactors = NaN*ones(nMons,4);
% meanFitRSquared = NaN*ones(nMons,1);
% 
% if (useMixerOverSqrtDiode == true)
%     meanFrascatiOutput = nanmean(nanmean(scanMixerSqrtDiode(:,:,:,calStartSample:calEndSample),4),3);
%     stdFrascatiOutput = nanstd(nanmean(scanMixerSqrtDiode(:,:,:,calStartSample:calEndSample),4),0,3); % NOTE - this not currently used for anything. Could/should use it for fitting weights.
% else
%     meanFrascatiOutput = nanmean(nanmean(scanMixers(:,:,:,calStartSample:calEndSample),4),3);
%     stdFrascatiOutput = nanstd(nanmean(scanMixers(:,:,:,calStartSample:calEndSample),4),0,3); % NOTE - this not currently used for anything. Could/should use it for fitting weights.
% end
% 
% % find max and min points in scan (max and min on sine curve). 
% [maxFrascatiOutput,maxFrascatiOutputIndex] = max(meanFrascatiOutput,[],2);
% [minFrascatiOutput,minFrascatiOutputIndex] = min(meanFrascatiOutput,[],2);
%     
% for mon=1:nMons
%     
%     % sensible guesses for fit parameters
%     initA = nanmean(abs([maxFrascatiOutput(mon) minFrascatiOutput(mon)])); % amplitude
%     initB = 4*(pi./180); % phase units, 4x for 3GHz -> 12GHz
%     initC = (meanFrascatiOutput(mon,1)./initA)*(pi/2); % phase shift
%     initD = (abs(maxFrascatiOutput(mon)) - abs(minFrascatiOutput(mon)))./2; % vertical offset
%     initParams = [initA, initB, initC, initD];
%     
%     % calculate fit
%     [scanSinFit,scanSinGoF] = offsetSinFit(scanPhShiftValues', meanFrascatiOutput(mon,:)', initParams);
%     meanFitRSquared(mon) = scanSinGoF.rsquare;
%     
%     % save fit parameters
%     meanCalibrationFactors(mon,1) = scanSinFit.a;
%     meanCalibrationFactors(mon,2) = scanSinFit.b;
%     meanCalibrationFactors(mon,3) = scanSinFit.c;
%     meanCalibrationFactors(mon,4) = scanSinFit.d;
%     fitFunc = @(x) (meanCalibrationFactors(mon,1)).*sin( (meanCalibrationFactors(mon,2)).*x + (meanCalibrationFactors(mon,3)) ) + meanCalibrationFactors(mon,4);
%     
%     % find zero crossing of the fitted curve (guess half way between
%     % min/max)
%     guessZeroCrossIndex = abs(maxFrascatiOutputIndex(mon)-minFrascatiOutputIndex(mon))./2; % index difference between min and max points
%     guessZeroCrossIndex = guessZeroCrossIndex + min([minFrascatiOutputIndex(mon) maxFrascatiOutputIndex(mon)]); % plus the starting index
%     guessZeroCrossIndex = round(guessZeroCrossIndex); % rounded to the nearest integer
%     guessZeroCrossValue = scanPhShiftValues(guessZeroCrossIndex); % phase shifter value at this index
%     % to get zero crossing use fit function without vertical offset,
%     % otherwise zero crossing fit func ~= zero crossing sine curve
%     % (actually looking for sin(x) = d, not sin(x) = 0)
%     zeroFunc = @(x) (meanCalibrationFactors(mon,1)).*sin( (meanCalibrationFactors(mon,2)).*x + (meanCalibrationFactors(mon,3)) ); 
%     meanZeroCrossings(mon) = fzero(zeroFunc, guessZeroCrossValue);
%     
%     % plot results
%     figure;
%     plot(scanPhShiftValues,meanFrascatiOutput(mon,:),'.b');
%     hold all;
%     plot([meanZeroCrossings(mon) meanZeroCrossings(mon)], [minFrascatiOutput(mon) maxFrascatiOutput(mon)],'k','LineWidth',2);
%     fitPlotX = linspace(min(scanPhShiftValues),max(scanPhShiftValues),720);
%     fitPlotY = fitFunc(fitPlotX);
%     plot(fitPlotX, fitPlotY,'r','LineWidth',2);
%     title(monMixerNames{mon});
%     xlabel('Phase shifter setting [3GHz degrees]');
%     if (useMixerOverSqrtDiode == 1)
%         ylabel('Mixer/sqrt(Diode)');
%     else
%         ylabel('Mixer [V]');
%     end
%     grid on;
%     legend('Data','Zero Crossing','Fit','Location','best');
%     title(sprintf('Amplitude: %.3f, Zero Crossing: %.3f degrees',meanCalibrationFactors(mon,1),meanZeroCrossings(mon)));
%     %savePngEpsFig(saveDir,sprintf('MeanPhaseFit_%s',monitorNames{mon}),1);
% end
% 
% meanCalibrationFactors

fprintf(1, 'Calibration with mean pulse phase:\n');

meanZeroCrossings = NaN*ones(nMons,1);
meanCalibrationFactors = NaN*ones(nMons,4);
meanFitRSquared = NaN*ones(nMons,1);

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
    meanFitRSquared(mon) = scanSinGoF.rsquare;
    
    % save fit parameters
    meanCalibrationFactors(mon,1) = scanSinFit.a;
    meanCalibrationFactors(mon,2) = scanSinFit.b;
    meanCalibrationFactors(mon,3) = scanSinFit.c;
    meanCalibrationFactors(mon,4) = scanSinFit.d;
    fitFunc = @(x) (meanCalibrationFactors(mon,1)).*sin( (meanCalibrationFactors(mon,2)).*x + (meanCalibrationFactors(mon,3)) ) + meanCalibrationFactors(mon,4);
    
    % find zero crossing of the fitted curve (guess half way between
    % min/max)
    guessZeroCrossValue = scanPhShiftValues(guessZeroCrossIndex); % phase shifter value at this index
    % to get zero crossing use fit function without vertical offset,
    % otherwise zero crossing fit func ~= zero crossing sine curve
    % (actually looking for sin(x) = d, not sin(x) = 0)
    zeroFunc = @(x) (meanCalibrationFactors(mon,1)).*sin( (meanCalibrationFactors(mon,2)).*x + (meanCalibrationFactors(mon,3)) ); 
    meanZeroCrossings(mon) = fzero(zeroFunc, guessZeroCrossValue);
    
    % plot results
    figure;
    plot(scanPhShiftValues,meanFrascatiOutput,'.b');
    hold all;
    plot([meanZeroCrossings(mon) meanZeroCrossings(mon)], [minFrascatiOutput maxFrascatiOutput],'k','LineWidth',2);
    fitPlotX = linspace(min(scanPhShiftValues),max(scanPhShiftValues),720);
    fitPlotY = fitFunc(fitPlotX);
    plot(fitPlotX, fitPlotY,'r','LineWidth',2);
    title(monMixerNames{mon});
    xlabel('Phase shifter setting (degrees)');
    if (useMixerOverSqrtDiode == 1)
        ylabel('Mixer/sqrt(Diode)');
    else
        ylabel('Mixer [V]');
    end
    grid on;
    legend('Data','Zero Crossing','Fit');
    
end

% save results
save([dataDir '/' dataFile '_calResults_mean.mat'],'calStartSamples','calEndSamples','meanCalibrationFactors','meanZeroCrossings','meanFitRSquared','dataSetComment');   
zeroCrossFileBackup = fopen([dataDir '/' dataFile '_NewZeroCrossings'],'w');
calibConstFileBackup = fopen([dataDir '/' dataFile '_NewCalibrationConstants'],'w');
for mon=1:nMons 
    fprintf(zeroCrossFileBackup,'%s, %.2f\n', phShiftNames{mon}, meanZeroCrossings(mon));
    fprintf(calibConstFileBackup,'%s, %f, %f, %f, %f\n', monMixerNames{mon}, meanCalibrationFactors(mon,1), meanCalibrationFactors(mon,2), meanCalibrationFactors(mon,3), meanCalibrationFactors(mon,4));
end
fclose(zeroCrossFileBackup);
fclose(calibConstFileBackup);
    
meanCalibrationFactors


%% Linear calibration on mean phase
% NB: NEED TO BE ABLE TO DEAL WITH ZERO CROSSING BEING ON DIFFERENT
% SLOPES/DATA FOR FIT BEING SPLIT IN HALF AT EACH END OF SCAN.
if (doLinearCalibrations)
    fprintf(1, 'Linear calibration with mean pulse phase:\n');

    meanLinearCalibrationFactors = NaN*ones(nMons,2);
    meanLinearFitRSquared = NaN*ones(nMons,1);

    if (useMixerOverSqrtDiode == true)
        meanFrascatiOutput = nanmean(nanmean(scanMixerSqrtDiode(:,:,:,calStartSample:calEndSample),4),3);
        stdFrascatiOutput = nanstd(nanmean(scanMixerSqrtDiode(:,:,:,calStartSample:calEndSample),4),0,3); % NOTE - this not currently used for anything. Could/should use it for fitting weights.
    else
        meanFrascatiOutput = nanmean(nanmean(scanMixers(:,:,:,calStartSample:calEndSample),4),3);
        stdFrascatiOutput = nanstd(nanmean(scanMixers(:,:,:,calStartSample:calEndSample),4),0,3); % NOTE - this not currently used for anything. Could/should use it for fitting weights.
    end

    rangeShifterOffset = (linearRange.*3)./12; % convert linear range to use in to frequency of phase shifter (3GHz)


    for mon=1:nMons
        isInRange = (scanPhShiftValues >= meanZeroCrossings(mon)-rangeShifterOffset) & (scanPhShiftValues <= meanZeroCrossings(mon)+rangeShifterOffset);
        rangePhaseValues = scanPhShiftValues(isInRange)-meanZeroCrossings(mon);
        rangePhaseValues = (rangePhaseValues.*12)./3; % phase offset from zero crossing in 12GHz degrees
        rangeOutput = meanFrascatiOutput(mon, isInRange);

        [fitresult,gof] = linearFit(rangePhaseValues, rangeOutput');
        meanLinearCalibrationFactors(mon,1) = fitresult.p1;
        meanLinearCalibrationFactors(mon,2) = fitresult.p2;
        meanLinearFitRSquared(mon) = gof.rsquare;

        % plot results
        figure;
        plot(rangePhaseValues,rangeOutput,'.b');
        hold all;
        fitPlotX = linspace(-linearRange,linearRange,720);
        fitPlotY = meanLinearCalibrationFactors(mon,1).*fitPlotX + meanLinearCalibrationFactors(mon,2);
        plot(fitPlotX, fitPlotY,'r','LineWidth',2);
        xlabel('Offset From Zero Crossing [12GHz degrees]');
        if (useMixerOverSqrtDiode == 1)
            ylabel('Mixer/sqrt(Diode)');
        else
            ylabel('Mixer [V]');
        end
        grid on;
        legend('Data','Fit','Location','best');
        title(sprintf('%s Mean, Grad: %.4f, RSq: %.3f',monitorNames{mon},meanLinearCalibrationFactors(mon,1),meanLinearFitRSquared(mon)));
        xlim([-linearRange linearRange]);
        ylim([-1 1]);
        hold off;
        %savePngEpsFig(saveDir, sprintf('%s_MeanLinearFit',monitorNames{mon}), 1);
    end


    % save results
    save([dataDir '/' dataFile '_linearCalResults_mean.mat'],'calStartSample','calEndSample','meanLinearCalibrationFactors','meanZeroCrossings','meanLinearFitRSquared','dataSetComment');
end


%% CALIBRATION FOR ALL SAMPLES

if (doPerSampleCalibrations)
    fprintf(1, 'Calibration with for all samples:\n');

    meanZeroCrossings = NaN*ones(nMons,nCalSamples);
    meanCalibrationFactors = NaN*ones(nMons,nCalSamples,4);
    fitRSquared = NaN*ones(nMons,nSamples);

    figure;
    for s=calStartSample:calEndSample
        fprintf('%d\n',s);
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
            initA = nanmean(abs([maxFrascatiOutput(mon) minFrascatiOutput(mon)])); % amplitude
            initB = 4*(pi./180); % phase units, 4x for 3GHz -> 12GHz
            initC = (meanFrascatiOutput(mon,1)./initA)*(pi/2); % phase shift
            initD = (abs(maxFrascatiOutput(mon)) - abs(minFrascatiOutput(mon)))./2; % vertical offset
            initParams = [initA, initB, initC, initD];

            % calculate fit
            [scanSinFit,scanSinGoF] = offsetSinFit(scanPhShiftValues', meanFrascatiOutput(mon,:)', initParams);

            % save fit parameters
            fitRSquared(mon,s) = scanSinGoF.rsquare;

            meanCalibrationFactors(mon,s,1) = scanSinFit.a;
            meanCalibrationFactors(mon,s,2) = scanSinFit.b;
            meanCalibrationFactors(mon,s,3) = scanSinFit.c;
            meanCalibrationFactors(mon,s,4) = scanSinFit.d;

            fitFunc = @(x) (meanCalibrationFactors(mon,s,1)).*sin( (meanCalibrationFactors(mon,s,2)).*x + (meanCalibrationFactors(mon,s,3)) ) + meanCalibrationFactors(mon,s,4);


            % find zero crossing of the fitted curve (guess half way between
            % min/max)
            guessZeroCrossIndex = abs(maxFrascatiOutputIndex(mon)-minFrascatiOutputIndex(mon))./2; % index difference between min and max points
            guessZeroCrossIndex = guessZeroCrossIndex + min([minFrascatiOutputIndex(mon) maxFrascatiOutputIndex(mon)]); % plus the starting index
            guessZeroCrossIndex = round(guessZeroCrossIndex); % rounded to the nearest integer
            guessZeroCrossValue = scanPhShiftValues(guessZeroCrossIndex); % phase shifter value at this index
            % to get zero crossing use fit function without vertical offset,
            % otherwise zero crossing fit func ~= zero crossing sine curve
            % (actually looking for sin(x) = d, not sin(x) = 0)
            zeroFunc = @(x) (meanCalibrationFactors(mon,s,1)).*sin( (meanCalibrationFactors(mon,s,2)).*x + (meanCalibrationFactors(mon,s,3)) ); 
            meanZeroCrossings(mon,s) = fzero(zeroFunc, guessZeroCrossValue);

            % plot results
    %         plot(scanPhShiftValues,meanFrascatiOutput(mon,:),'.b');
    %         hold all;
    %         plot([zeroCrossings(mon,s) zeroCrossings(mon,s)], [minFrascatiOutput(mon) maxFrascatiOutput(mon)],'k','LineWidth',2);
    %         fitPlotX = linspace(min(scanPhShiftValues),max(scanPhShiftValues),720);
    %         fitPlotY = fitFunc(fitPlotX);
    %         plot(fitPlotX, fitPlotY,'r','LineWidth',2);
    %         title(sprintf('%s, Samp%03d, Ampl: %.3f',monitorNames{mon},s,calibrationFactors(mon,s,1)));
    %         xlabel('Phase shifter setting [3GHz degrees]');
    %         if (useMixerOverSqrtDiode == 1)
    %             ylabel('Mixer/sqrt(Diode)');
    %         else
    %             ylabel('Mixer [V]');
    %         end
    %         grid on;
    %         legend('Data','Zero Crossing','Fit','Location','best');
    %         title(sprintf('%s, Samp %03d, Ampl: %.3f',monitorNames{mon},s,calibrationFactors(mon,s,1)));
    %         hold off;
    %         %savePngEpsFig([saveDir '/allSamples'],sprintf('%s_Sample_%03d',monitorNames{mon},s),1);

        end

    end

    % save results
    save([dataDir '/' dataFile '_calResults_allSamples.mat'],'calStartSample','calEndSample','calibrationFactors','zeroCrossings','fitRSquared','dataSetComment');

    % plot results
    calibFactorNames = {'a','b','c','d'};

    for c=1:4 % calibration factor index
        figure;
        for mon=1:nMons
            plot(squeeze(meanCalibrationFactors(mon,:,c)));
            hold all;
        end
        title(calibFactorNames{c});
        legend(monitorNames,'Location','best');
        xlabel('Sample No.');
        ylabel('Parameter Value');
        %savePngEpsFig(saveDir,sprintf('fitParam_%s',calibFactorNames{c}),1);
    end

    figure;
    for mon=1:nMons
        plot(squeeze(meanZeroCrossings(mon,:)));
        hold all;
    end
    title('Zero Crossing');
    legend(monitorNames,'Location','best');
    xlabel('Sample No.');
    ylabel('Phase [degrees]');
    %savePngEpsFig(saveDir,'ZeroCrossing',1);

    figure;
    for mon=1:nMons
        plot(fitRSquared(mon,:));
        hold all;
    end
    title('Fit R Squared');
    xlabel('Sample No.');
    ylabel('R squared');
    legend(monitorNames,'Location','best');
    %%savePngEpsFig(saveDir,'RSquared',1);
end

%% Linear calibration on all samples
    if (doLinearCalibrations && doPerSampleCalibrations)

    fprintf(1, 'Linear calibration with mean pulse phase:\n');

    rangeShifterOffset = (linearRange.*3)./12; % convert 12GHz offset in to 3GHz of phase shifter

    linearCalibrationFactors = NaN*ones(nMons,nCalSamples,2);
    linearFitRSquared = NaN*ones(nMons,nSamples);

    figure;
    for s=calStartSample:calEndSample
        fprintf('%d\n',s);
        if (useMixerOverSqrtDiode == true)
            meanFrascatiOutput = nanmean(scanMixerSqrtDiode(:,:,:,s),3);
            stdFrascatiOutput = nanstd(scanMixerSqrtDiode(:,:,:,s),0,3); % NOTE - this not currently used for anything. Could/should use it for fitting weights.
        else
            meanFrascatiOutput = nanmean(scanMixers(:,:,:,s),3);
            stdFrascatiOutput = nanstd(scanMixers(:,:,:,s),0,3); % NOTE - this not currently used for anything. Could/should use it for fitting weights.
        end

        for mon=1:nMons
            % linear fit

            isInRange = (scanPhShiftValues >= meanZeroCrossings(mon,s)-rangeShifterOffset) & (scanPhShiftValues <= meanZeroCrossings(mon,s)+rangeShifterOffset);
            rangePhaseValues = scanPhShiftValues(isInRange)-meanZeroCrossings(mon,s);
            rangePhaseValues = (rangePhaseValues.*12)./3; % phase offset from zero crossing in 12GHz degrees
            rangeOutput = meanFrascatiOutput(mon,isInRange);

            [fitresult,gof] = linearFit(rangePhaseValues, rangeOutput);
            linearCalibrationFactors(mon,s,1) = fitresult.p1;
            linearCalibrationFactors(mon,s,2) = fitresult.p2;
            linearFitRSquared(mon,s) = gof.rsquare;

            % plot results
            plot(rangePhaseValues,rangeOutput,'.b');
            hold all;
            fitPlotX = linspace(-linearRange,linearRange,720);
            fitPlotY = linearCalibrationFactors(mon,s,1).*fitPlotX + linearCalibrationFactors(mon,s,2);
            plot(fitPlotX, fitPlotY,'r','LineWidth',2);
            xlabel('Offset From Zero Crossing [12GHz degrees]');
            if (useMixerOverSqrtDiode == 1)
                ylabel('Mixer/sqrt(Diode)');
            else
                ylabel('Mixer [V]');
            end
            grid on;
            legend('Data','Fit','Location','best');
            title(sprintf('%s, Samp: %03d, Grad: %.4f, RSq: %.3f',monitorNames{mon},s,linearCalibrationFactors(mon,s,1),linearFitRSquared(mon,s)));
            xlim([-linearRange linearRange]);
            ylim([-2 2]);
            hold off;
            %%savePngEpsFig([saveDir '/linearFitAllSamples'], sprintf('%s_Sample_%03d',monitorNames{mon},s), 1);

        end

    end

    % save results
    save([dataDir '/' dataFile '_linearCalResults_allSamples.mat'],'calStartSample','calEndSample','linearCalibrationFactors','zeroCrossings','linearFitRSquared','dataSetComment');

    % plot results
    calibFactorNames = {'Gradient','Offset'};

    for c=1:2 % calibration factor index
        figure;
        for mon=1:nMons
            plot(squeeze(linearCalibrationFactors(mon,:,c)));
            hold all;
        end
        title(calibFactorNames{c});
        legend(monitorNames,'Location','best');
        xlabel('Sample No.');
        ylabel('Parameter Value');
        %savePngEpsFig(saveDir,sprintf('linearFitParam_%s',calibFactorNames{c}),1);
    end

    figure;
    for mon=1:nMons
        plot(linearFitRSquared(mon,:));
        hold all;
    end
    title('Linear Fit R Squared');
    xlabel('Sample No.');
    ylabel('R squared');
    legend(monitorNames,'Location','best');
    %savePngEpsFig(saveDir,'linearFitRSquared',1);
end
%% Other plots
if (doPerSampleCalibrations)
    for mon=1:nMons
        % max, mean and min fits, sine
        figure;
        meanFrascatiOutput = nanmean(nanmean(scanMixerSqrtDiode(mon,:,:,calStartSample:calEndSample),4),3);
        legMean = sprintf('Mean (ampl=%.3f)',meanCalibrationFactors(mon,1));
        fitFunc = @(x) (meanCalibrationFactors(mon,1)).*sin( (meanCalibrationFactors(mon,2)).*x + (meanCalibrationFactors(mon,3)) ) + meanCalibrationFactors(mon,4);
        plot(scanPhShiftValues,meanFrascatiOutput,'.b');
        hold all;
        fitPlotX = linspace(min(scanPhShiftValues),max(scanPhShiftValues),720);
        fitPlotY = fitFunc(fitPlotX);
        meanPlot = plot(fitPlotX, fitPlotY,'b','LineWidth',2);

        [~,maxSample] = max(squeeze(meanCalibrationFactors(mon,:,1)));
        legMax = sprintf('Max (samp%d, ampl=%.3f)',maxSample,meanCalibrationFactors(mon,maxSample,1));
        meanFrascatiOutput = nanmean(scanMixerSqrtDiode(mon,:,:,maxSample),3);
        fitFunc = @(x) (meanCalibrationFactors(mon,maxSample,1)).*sin( (meanCalibrationFactors(mon,maxSample,2)).*x + (meanCalibrationFactors(mon,maxSample,3)) ) + meanCalibrationFactors(mon,maxSample,4);
        plot(scanPhShiftValues,meanFrascatiOutput,'.r');
        fitPlotX = linspace(min(scanPhShiftValues),max(scanPhShiftValues),720);
        fitPlotY = fitFunc(fitPlotX);
        maxPlot = plot(fitPlotX, fitPlotY,'r','LineWidth',2);

        [~,minSample] = min(squeeze(meanCalibrationFactors(mon,:,1)));
        legMin = sprintf('Min (samp%d, ampl=%.3f)',minSample,meanCalibrationFactors(mon,minSample,1));
        meanFrascatiOutput = nanmean(scanMixerSqrtDiode(mon,:,:,minSample),3);
        fitFunc = @(x) (meanCalibrationFactors(mon,minSample,1)).*sin( (meanCalibrationFactors(mon,minSample,2)).*x + (meanCalibrationFactors(mon,minSample,3)) ) + meanCalibrationFactors(mon,minSample,4);
        plot(scanPhShiftValues,meanFrascatiOutput,'.g');
        fitPlotX = linspace(min(scanPhShiftValues),max(scanPhShiftValues),720);
        fitPlotY = fitFunc(fitPlotX);
        minPlot = plot(fitPlotX, fitPlotY,'g','LineWidth',2);

        xlabel('Phase shifter setting [1GHz degrees]');
        ylabel('Mixer/sqrt(Diode)');
        legend([minPlot meanPlot maxPlot], {legMin,legMean,legMax},'Location','best');
        title(sprintf('%s Min, Mean, Max Fit Amplitude',monitorNames{mon}));
        grid on;

        %savePngEpsFig(saveDir,sprintf('MaxMeanMinFit_%s',monitorNames{mon}),1);

        % max, mean and min fits, linear
        if (doLinearCalibrations)
            figure;
            meanFrascatiOutput = nanmean(nanmean(scanMixerSqrtDiode(mon,:,:,calStartSample:calEndSample),4),3);
            legMean = sprintf('Mean (grad=%.4f)', meanLinearCalibrationFactors(mon,1));
            isInRange = (scanPhShiftValues >= meanZeroCrossings(mon)-rangeShifterOffset) & (scanPhShiftValues <= meanZeroCrossings(mon)+rangeShifterOffset);
            rangePhaseValues = scanPhShiftValues(isInRange)-meanZeroCrossings(mon);
            rangePhaseValues = (rangePhaseValues.*12)./3; % phase offset from zero crossing in 12GHz degrees
            rangeOutput = meanFrascatiOutput(isInRange);
            plot(rangePhaseValues,rangeOutput,'.b');
            hold all;
            fitPlotX = linspace(-linearRange,linearRange,720);
            fitPlotY = meanLinearCalibrationFactors(mon,1).*fitPlotX + meanLinearCalibrationFactors(mon,2);
            meanPlot = plot(fitPlotX, fitPlotY,'b','LineWidth',2);

            [~,maxSample] = max(squeeze(linearCalibrationFactors(mon,:,1)));
            legMax = sprintf('Max (samp%d, grad=%.4f)',maxSample,linearCalibrationFactors(mon,maxSample,1));
            meanFrascatiOutput = nanmean(scanMixerSqrtDiode(mon,:,:,maxSample),3);
            isInRange = (scanPhShiftValues >= meanZeroCrossings(mon,maxSample)-rangeShifterOffset) & (scanPhShiftValues <= meanZeroCrossings(mon,maxSample)+rangeShifterOffset);
            rangePhaseValues = scanPhShiftValues(isInRange)-meanZeroCrossings(mon,maxSample);
            rangePhaseValues = (rangePhaseValues.*12)./3; % phase offset from zero crossing in 12GHz degrees
            rangeOutput = meanFrascatiOutput(isInRange);
            plot(rangePhaseValues,rangeOutput,'.r');
            hold all;
            fitPlotX = linspace(-linearRange,linearRange,720);
            fitPlotY = linearCalibrationFactors(mon,maxSample,1).*fitPlotX + linearCalibrationFactors(mon,maxSample,2);
            maxPlot = plot(fitPlotX, fitPlotY,'r','LineWidth',2);

            [~,minSample] = min(squeeze(linearCalibrationFactors(mon,:,1)));
            legMin = sprintf('Min (samp%d, grad=%.4f)',minSample,linearCalibrationFactors(mon,minSample,1));
            meanFrascatiOutput = nanmean(scanMixerSqrtDiode(mon,:,:,minSample),3);
            isInRange = (scanPhShiftValues >= meanZeroCrossings(mon,minSample)-rangeShifterOffset) & (scanPhShiftValues <= meanZeroCrossings(mon,minSample)+rangeShifterOffset);
            rangePhaseValues = scanPhShiftValues(isInRange)-meanZeroCrossings(mon,minSample);
            rangePhaseValues = (rangePhaseValues.*12)./3; % phase offset from zero crossing in 12GHz degrees
            rangeOutput = meanFrascatiOutput(isInRange);
            plot(rangePhaseValues,rangeOutput,'.g');
            hold all;
            fitPlotX = linspace(-linearRange,linearRange,720);
            fitPlotY = linearCalibrationFactors(mon,minSample,1).*fitPlotX + linearCalibrationFactors(mon,minSample,2);
            minPlot = plot(fitPlotX, fitPlotY,'g','LineWidth',2);

            xlabel('Offset From Zero Crossing [12GHz degrees]');
            ylabel('Mixer/sqrt(Diode)');
            grid on;
            legend([minPlot meanPlot maxPlot], {legMin,legMean,legMax},'Location','best');
            title(sprintf('%s Min, Mean, Max Fit Gradient',monitorNames{mon}));
            xlim([-linearRange linearRange]);
            ylim([-2 2]);

            %savePngEpsFig(saveDir,sprintf('linearMaxMeanMinFit_%s',monitorNames{mon}),1);
        end
    end
end