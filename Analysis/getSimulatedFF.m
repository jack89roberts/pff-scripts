function [ simFFResults ] = getSimulatedFF( processedData, gain, limit, offset )
%getSimulatedFF Takes processed data and uses phases to calculate simulated
%feedforward results. OUTPUT IS ALWAYS MEAN SUBTRACTED!!
%
%   simFFResults - output struct containing calculated statistics etc.
%
%   processedData - struct containing processed data
%
%   gain - if given uses this gain in calculation. Can either be one value
%   or one value for each sample. If not given/empty calculates
%   optimal gain on mean and uses that. If 'sample' calculates optimal
%   sample by sample gains and uses them.
%
%   limit - if given limits the correction to +/-limit degrees. If not
%   given/empty uses unlimited correction.
%
%   offset - if given applies an offset to the upstream phase so it is not
%   centred in the correction range. Has no effect if limit not given/infinite.
%   Offset is then removed again in output.

%     fprintf('Simulating FF results...\n');
    
    % extract relevant variables from data struct
    if (processedData.subtractPhase(3)==0 && processedData.subtractPhase(2)==0)
        % need to calculate and apply mean subtraction
        subtractMon2 = nanmean(processedData.meanPulsePhase(2,:));
        subtractMon3 = nanmean(processedData.meanPulsePhase(3,:));
        mon3 = squeeze(processedData.phases(3,:,:))-subtractMon3;
        mon2 = squeeze(processedData.phases(2,:,:))-subtractMon2;
        meanMon2 = processedData.meanPulsePhase(2,:)-subtractMon2;
    else
        % data already mean subtracted (but need to apply to raw phases)
        mon3 = squeeze(processedData.phases(3,:,:))-processedData.subtractPhase(3);
        mon2 = squeeze(processedData.phases(2,:,:))-processedData.subtractPhase(2);
        meanMon2 = processedData.meanPulsePhase(2,:);
    end
    sampleRange = processedData.sampleRange;
    
    % calculate gain
    if (nargin<2 || isempty(gain))
        %  use optimal gain based on mean statistics
        corr = processedData.corrMeanMix2Mix3;
        jitMon3 = processedData.stdMeanPulsePhase(3);
        jitMon2 = processedData.stdMeanPulsePhase(2);
        gain = corr*(jitMon3/jitMon2);
        
    elseif (strcmp(gain,'sample'))
        % use optimal gains based on sample by sample statistics
        corr = processedData.corrSamplesMix2Mix3;
        jitMon3 = processedData.stdPhaseAlongPulse(3,:);
        jitMon2 = processedData.stdPhaseAlongPulse(2,:);
        gain = corr.*(jitMon3./jitMon2);
    
    elseif (~isnumeric(gain))
        % some string other than 'sample' was given
        error('invalid gain value')
        
    elseif (length(gain)~=1 && length(gain)~=processedData.nSamples)
        % number of gain values not 1 or number of samples
        error('invalid gain length')
    end
    
    % apply offset if given
    if (nargin>3 && ~isempty(offset))
        mon2 = mon2+offset;
    end
    
    % calculate correction
    if length(gain)==1
        simCorrection = gain.*mon2;
    else
        simCorrection = bsxfun(@times,mon2,gain); % multiplies each mon2 pulse by gain vector
    end
    
    % apply limit to correction
    if (nargin>2 && ~isempty(limit))
%         fprintf('Limiting correction...\n');
        simCorrection(simCorrection > limit) = limit;
        simCorrection(simCorrection < -limit) = -limit;
    end
    
    % apply correction, calculate stats
%     fprintf('Applying simulated correction...\n');
    simFFPhases = mon3-simCorrection;
    [meanSimFFPhase,flatnessSimFF,meanSimFFPhase_err,flatnessSimFF_err] = nanMeanStdErr(simFFPhases(:,sampleRange),2);
    
    [~,stdSimFFPhase,~,stdSimFFPhase_err] = nanMeanStdErr(meanSimFFPhase);

    if (nargin>3 && ~isempty(offset))
        % remove offset in simulation output
        offsetMeanSim = nanmean(meanSimFFPhase);
        simFFPhases = simFFPhases-offsetMeanSim;
        meanSimFFPhase = meanSimFFPhase-offsetMeanSim;
    end
    
    [meanSimFFAlongPulse,stdSimFFAlongPulse,meanSimFFAlongPulse_err,stdSimFFAlongPulse_err] = nanMeanStdErr(simFFPhases,1);
    [corrMeanSimFF,corrMeanSimFF_err] = nancorrcoef(meanMon2,meanSimFFPhase);
    
    % save results to struct
    simFFResults = struct();
    simFFResults.gain = gain;
    if (nargin>2 && ~isempty(limit))
        simFFResults.limit = limit;
    else
        simFFResults.limit = Inf;
    end
    if (nargin>3 && ~isempty(offset))
        simFFResults.offset = offset;
    else
        simFFResults.offset= 0;
    end
    simFFResults.simFFPhases = simFFPhases;
    simFFResults.meanSimFFPhase = meanSimFFPhase;
    simFFResults.stdSimFFPhase = stdSimFFPhase;
    simFFResults.meanSimFFPhase_err = meanSimFFPhase_err;
    simFFResults.stdSimFFPhase_err = stdSimFFPhase_err;
    simFFResults.meanSimFFAlongPulse = meanSimFFAlongPulse;
    simFFResults.stdSimFFAlongPulse = stdSimFFAlongPulse;
    simFFResults.meanSimFFAlongPulse_err = meanSimFFAlongPulse_err;
    simFFResults.stdSimFFAlongPulse_err = stdSimFFAlongPulse_err;
    simFFResults.corrMeanSimFF = corrMeanSimFF;
    simFFResults.corrMeanSimFF_err = corrMeanSimFF_err;
    simFFResults.flatnessSimFF = flatnessSimFF;
    simFFResults.flatnessSimFF_err = flatnessSimFF_err;
    simFFResults.simCorrection = simCorrection;
    
end

