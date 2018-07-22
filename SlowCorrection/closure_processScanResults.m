function scanResults = closure_processScanResults(corr465Offsets,corr765Offsets,bpmClosureIndices,initBPMH,initBPMS, scanBPMH,scanBPMS )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    scanCorrRatios = corr765Offsets(1,:)./corr465Offsets(1);
    n465Offsets = length(corr465Offsets);
    n765Offsets = length(scanCorrRatios);
    nBPMs = size(initBPMH,2);
    
    % align BPM output
    for indBPM=1:nBPMs
        [tmpBPMS,tmpBPMH] = getAlignedXCorr(squeeze(initBPMS(:,indBPM,:)),'end',{squeeze(initBPMH(:,indBPM,:))});
        tmpBPMH = tmpBPMH{1};

        initBPMS(:,indBPM,:) = tmpBPMS;
        initBPMH(:,indBPM,:) = tmpBPMH;
        
        for ind465 = 1:n465Offsets
            for ind765= 1:n765Offsets
                tmpBPMS = squeeze(scanBPMS(ind465,ind765,:,indBPM,:));
                tmpBPMH = squeeze(scanBPMH(ind465,ind765,:,indBPM,:));
                [tmpBPMS,tmpBPMH] = getAlignedXCorr(tmpBPMS,'end',{tmpBPMH});
                tmpBPMH = tmpBPMH{1};
                
                scanBPMS(ind465,ind765,:,indBPM,:) = tmpBPMS;
                scanBPMH(ind465,ind765,:,indBPM,:) = tmpBPMH;
            end
        end
            
    end
    
    % get bpm sample range
    bpmSampleRange = getBPMSampleRange(squeeze(initBPMS(1,:,:)),squeeze(initBPMH(1,:,:)));
    
    % remove bad pulses - to add
    % check transmission - to add
    
    % calculate init orbit
    meanInitH = nanmean(nanmean(initBPMH(:,:,bpmSampleRange)));
    meanInitS = nanmean(nanmean(initBPMS(:,:,bpmSampleRange)));
    stdInitH = nanstd(nanmean(initBPMH(:,:,bpmSampleRange)));
    stdInitS = nanstd(nanmean(initBPMH(:,:,bpmSampleRange)));
            
    % calculate mean orbit at each corrector setting
    meanScanBPMH = nanmean(nanmean(scanBPMH(:,:,:,:,bpmSampleRange),5),3);
    meanScanBPMS = nanmean(nanmean(scanBPMS(:,:,:,:,bpmSampleRange),5),3);
    stdScanBPMH = nanstd(nanmean(scanBPMH(:,:,:,:,bpmSampleRange),5),0,3);
    stdScanBPMS = nanstd(nanmean(scanBPMS(:,:,:,:,bpmSampleRange),5),0,3);

    % calculate difference orbit at each corrector setting
    meanDiffBPMH = NaN(n465Offsets,n765Offsets,nBPMs);
    meanDiffBPMS = NaN(n465Offsets,n765Offsets,nBPMs);
    stdDiffBPMH = NaN(n465Offsets,n765Offsets,nBPMs);
    stdDiffBPMS = NaN(n465Offsets,n765Offsets,nBPMs);
    
    for ind465Val=1:n465Offsets
        for ind765Val=1:n765Offsets
            for bpm=1:nBPMs
                meanDiffBPMH(ind465Val,ind765Val,bpm) = meanScanBPMH(ind465Val, ind765Val, bpm) - meanInitH(bpm);
                meanDiffBPMS(ind465Val,ind765Val,bpm) = meanScanBPMS(ind465Val, ind765Val, bpm) - meanInitS(bpm);

                stdDiffBPMH(ind465Val,ind765Val,bpm) = sqrt( stdScanBPMH(ind465Val,ind765Val,bpm).^2 + stdInitH(bpm).^2 );
                stdDiffBPMS(ind465Val,ind765Val,bpm) = sqrt( stdScanBPMS(ind465Val,ind765Val,bpm).^2 + stdInitS(bpm).^2 );
            end
        end
    end
    
    % calculate sum differences squared
    sumDiffSq = zeros(n465Offsets,n765Offsets);
    stdDiffSq = zeros(n465Offsets,n765Offsets); % d(a^2) = a^2 * sqrt((da/a)^2 + (da/a)^2), d(sum(a^2)) = sqrt(sum(d(a^2)^2))
    for ind465Val=1:n465Offsets
        for ind765Val=1:n765Offsets
            for bpm=bpmClosureIndices
                tmpDiffSq = meanDiffBPMH(ind465Val,ind765Val,bpm).^2;
                sumDiffSq(ind465Val,ind765Val) = sumDiffSq(ind465Val,ind765Val) + tmpDiffSq;
                stdDiffSq(ind465Val,ind765Val) = stdDiffSq(ind465Val,ind765Val) + (tmpDiffSq.^2).*( 2.* (stdDiffBPMH(ind465Val,ind765Val,bpm)./meanDiffBPMH(ind465Val,ind765Val,bpm)).^2) ;
            end
        end
    end
    stdDiffSq = sqrt(stdDiffSq);
    
    % fit sum differences squared vs. corr ratio
    fitCoeff = NaN(n465Offsets,3);
    fitCoeffStd = NaN(n465Offsets,3);
    fitGoF = NaN(1,n465Offsets);
    fit765Ratio = NaN(1,n465Offsets);
    fitRatioStd = NaN(1,n465Offsets);

    for ind465Val=1:n465Offsets
        weights = squeeze(1./(stdDiffSq(ind465Val,:).^2))';
        [fitResult, gof] = quadraticFit(scanCorrRatios', squeeze(sumDiffSq(ind465Val, :))', weights);

        fitCoeff(ind465Val,1) = fitResult.p1;
        fitCoeff(ind465Val,2) = fitResult.p2;
        fitCoeff(ind465Val,3) = fitResult.p3;
        fitGoF(ind465Val) = gof.rsquare;

        tmpConf = confint(fitResult);
        for c=1:3
            fitCoeffStd(ind465Val,c) = (fitCoeff(ind465Val,c) - tmpConf(1,c))./1.959964; % 1.96 to convert 95% conf. int. to std. error
        end

        fit765Ratio(ind465Val) = -fitCoeff(ind465Val,2)./(2.*fitCoeff(ind465Val,1));
        tmpDenomMean = 2.*fitCoeff(ind465Val,1);
        tmpNumerMean = fitCoeff(ind465Val,2);
        tmpDenomStd = 2.*fitCoeffStd(ind465Val,1);
        tmpNumerStd = fitCoeffStd(ind465Val,2);
        fitRatioStd(ind465Val) = fit765Ratio(ind465Val).*sqrt( (tmpNumerMean./tmpNumerStd).^2 + (tmpDenomMean./tmpDenomStd).^2 );
    end

    % get best estimate of ratio for orbit closure from mean of each 465
    % setting
    bestEstRatio = 0;
    sumWeights = 0;
    for ind465Val = 1:n465Offsets
        tmpWeight = 1./(fitRatioStd(ind465Val).^2);
        bestEstRatio = bestEstRatio + fit765Ratio(ind465Val).*tmpWeight;
        sumWeights = sumWeights + tmpWeight;
    end
    bestEstRatio = bestEstRatio./sumWeights

    % plot data
    for ind465Val = 1:n465Offsets
        fitFunc = @(x) fitCoeff(ind465Val,1).*x.^2 + fitCoeff(ind465Val,2).*x + fitCoeff(ind465Val,3);
        x = linspace(min(scanCorrRatios),max(scanCorrRatios),100);
        y = fitFunc(x);
        plot(x,y);
        hold all;
        plot(scanCorrRatios,squeeze(sumDiffSq(ind465Val,:)),'o');
    end
    
    % create output struct
    scanResults = struct();
    scanResults.scanCorrRatios=scanCorrRatios;
    scanResults.corr465Offsets=corr465Offsets;
    scanResults.corr765Offsets=corr765Offsets;
    scanResults.bpmClosureIndices=bpmClosureIndices;
    scanResults.initBPMH=initBPMH;
    scanResults.initBPMS=initBPMS;
    scanResults.scanBPMH=scanBPMH;
    scanResults.scanBPMS=scanBPMS;
    scanResults.bpmSampleRange=bpmSampleRange;
    scanResults.meanInitH=meanInitH;
    scanResults.meanInitS=meanInitS;
    scanResults.stdInitH=stdInitH;
    scanResults.stdInitS=stdInitS;
    scanResults.meanScanBPMH=meanScanBPMH;
    scanResults.meanScanBPMS=meanScanBPMS;
    scanResults.stdScanBPMH=stdScanBPMH;
    scanResults.stdScanBPMS=stdScanBPMS;
    scanResults.meanDiffBPMH=meanDiffBPMH;
    scanResults.meanDiffBPMS=meanDiffBPMS;
    scanResults.stdDiffBPMH=stdDiffBPMH;
    scanResults.stdDiffBPMS=stdDiffBPMS;
    scanResults.sumDiffSq=sumDiffSq;
    scanResults.stdDiffSq=stdDiffSq;
    scanResults.fitCoeff=fitCoeff;
    scanResults.fitCoeffStd=fitCoeffStd;
    scanResults.fitGoF=fitGoF;
    scanResults.fit765Ratio=fit765Ratio;
    scanResults.fitRatioStd=fitRatioStd;
    scanResults.bestEstRatio=bestEstRatio;

end

