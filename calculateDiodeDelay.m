function [ frascatiStaticDelay ] = calculateDiodeDelay( diodes, showResult )

    if (nargin<2)
        showResult = 0;
    end

    mon1Diode = squeeze(diodes(1,:,:));
    mon2Diode = squeeze(diodes(2,:,:));
    mon3Diode = squeeze(diodes(3,:,:));

%     mon1Diode = double(mon1Diode);
%     mon2Diode = double(mon2Diode);
%     mon3Diode = double(mon3Diode);
% 
    mon1Diode = removeBadPulses(mon1Diode);
    mon2Diode = removeBadPulses(mon2Diode);
    mon3Diode = removeBadPulses(mon3Diode);

    mon1Diode = getAlignedXCorr(mon1Diode,'first');
    mon2Diode = getAlignedXCorr(mon2Diode,'first');
    mon3Diode = getAlignedXCorr(mon3Diode,'first');

    [nPulses,nSamples] = size(mon2Diode);

    bestLags21 = NaN(1,nPulses);
    bestLags23 = NaN(1,nPulses);

    for p=1:nPulses

        if (sum(isnan(mon1Diode(p,:)))==0 && sum(isnan(mon2Diode(p,:)))==0)
            [corrs21,lags21] = xcorr(mon2Diode(p,:),mon1Diode(p,:));
        elseif (sum(isnan(mon1Diode(p,:)))<nSamples && sum(isnan(mon2Diode(p,:)))<nSamples)
            [corrs21,lags21] = nanxcorr(mon2Diode(p,:),mon1Diode(p,:));
        else
            corrs21 = NaN; lags21 = NaN;
        end
        [~,bestInd] = max(abs(corrs21));
        bestLags21(p) = lags21(bestInd);

        if (sum(isnan(mon2Diode(p,:)))==0 && sum(isnan(mon3Diode(p,:)))==0)
            [corrs23,lags23] = xcorr(mon2Diode(p,:),mon3Diode(p,:));
        elseif (sum(isnan(mon2Diode(p,:)))<nSamples && sum(isnan(mon3Diode(p,:)))<nSamples)
            [corrs23,lags23] = nanxcorr(mon2Diode(p,:),mon3Diode(p,:),300);
        else
            corrs23 = NaN; lags23 = NaN;
        end
        [~,bestInd] = max(abs(corrs23));
        bestLags23(p) = lags23(bestInd);   
    end

    [meanLag21,errLag21] =nanMeanStdErr(bestLags21);
    [meanLag23,errLag23] =nanMeanStdErr(bestLags23);
    frascatiStaticDelay = [meanLag21 0 meanLag23];
    fprintf('Delay 21: Mean %.2f, Error %.2f\n',meanLag21,errLag21);
    fprintf('Delay 23: Mean %.2f, Error %.2f\n',meanLag23,errLag23);
    
    if (showResult)
        figure;
        plot(bestLags21);
        hold all;
        plot(bestLags23);
        legend('Delay 2-1','Delay 3-2');
    end


end