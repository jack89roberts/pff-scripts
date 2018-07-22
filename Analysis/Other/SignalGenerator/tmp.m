m=1;
d=12;
stdResidVolt = NaN(nMixers,nPowerLevels);
stdResidPhas = NaN(nMixers,nPowerLevels);
for m=1:nMixers
    for d=1:nPowerLevels
        fitFunc = @(x) (mixerCalConsts(m,d,1).*sin( (mixerCalConsts(m,d,2).*x) + mixerCalConsts(m,d,3) )) + mixerCalConsts(m,d,4);
        fitX = 1:nSamples;
        fitY = fitFunc(fitX);
        stdResidVolt(m,d) = nanstd(squeeze(mixers(m,d,:))'-fitY);
        
        stdResidPhas(m,d) = asin(stdResidVolt(m,d)/mixerCalConsts(m,d,1));

%         figure;
%         plot(fitX,fitY,'k','LineWidth',1.5);
%         hold all;
%         plot(squeeze(mixers(m,d,:)),'o','Color',mixerColours(m,:));
%         title(sprintf('Mixer %d, Power %d dBm (Amplitude: %.2d V, Rsquare: %.3f)',m,powerLevels(d),abs(mixerCalConsts(m,d,1)),mixerCalRSquare(m,d)));
%         xlabel('Sample No.');
%         ylabel('Output [V]');
%         grid on;
%         legend('Fit','Data');
%         hold off;
%         

    end
end

figure;
y=squeeze(mixers(m,d,:))'-fitY;
y=asind(y./mixerCalConsts(m,d,1));
plot(1:nSamples,y,'LineWidth',2);
nanstd(squeeze(mixers(m,d,:))'-fitY)

figure;
plot(stdResidVolt');

figure;
plot(stdResidPhas');

figure;
plot((squeeze(mixers(m,d,:))'-fitY));

figure;
m=1;
d=10;
fitFunc = @(x) (mixerCalConsts(m,d,1).*sin( (mixerCalConsts(m,d,2).*x) + mixerCalConsts(m,d,3) )) + mixerCalConsts(m,d,4);
restFitFunc = @(x) (0.5654.*sin( (0.09355.*x) + 6.151 )) - 0.07204;
fitX = 1:nSamples;
fitY = fitFunc(fitX);
restFitY = restFitFunc(fitX);
plot(squeeze(mixers(m,d,:)),'o','Color',mixerColours(m,:));
hold all;
plot(fitX,fitY,'k','LineWidth',1.5);
plot(fitX,restFitY,'m--','LineWidth',1.5)
title(sprintf('Mixer %d, Power %d dBm (Amplitude: %.2d V, Rsquare: %.3f)',m,powerLevels(d),abs(mixerCalConsts(m,d,1)),mixerCalRSquare(m,d)));
xlabel('Sample No.');
ylabel('Output [V]');
grid on;
legend('Fit','Data');
hold off;

figure;
d=7;
for m=1:nMixers
plot(squeeze(mixers(m,d,:)),'o','Color',mixerColours(m,:));
hold all;
end
legend('Mixer 1','Mixer 2','Mixer 3');
for m=1:nMixers
   fitFunc = @(x) (mixerCalConsts(m,d,1).*sin( (mixerCalConsts(m,d,2).*x) + mixerCalConsts(m,d,3) )) + mixerCalConsts(m,d,4);
fitX = 1:nSamples;
fitY = fitFunc(fitX);
plot(fitX,fitY,'LineWidth',1.5,'Color',mixerColours(m,:));
hold all;
end
title(sprintf('Mixer Output at %d dBm ',powerLevels(d)));
xlabel('Sample No.');
ylabel('Output [V]');
xlim([92 187])
format_plots
hold off;

figure;
d=10;
for m=1:nMixers
fitX = 1:nSamples;
fitFunc = @(x) (mixerCalConsts(m,d,1).*sin( (mixerCalConsts(m,d,2).*x) + mixerCalConsts(m,d,3) )) + mixerCalConsts(m,d,4);
fitY = fitFunc(fitX);
plot(squeeze(mixers(m,d,:))-fitY','Color',mixerColours(m,:));
hold all;
end
legend('Mixer 1','Mixer 2','Mixer 3');
title(sprintf('Mixer Fit Residuals at %d dBm ',powerLevels(d)));
xlabel('Sample No.');
ylabel('Output [V]');
xlim([92 187])
format_plots
hold off;


fitPhas = asind((fitY-mixerCalConsts(m,d,4))./mixerCalConsts(m,d,1));
datPhase = real(asind((squeeze(mixers(m,d,:))'-mixerCalConsts(m,d,4))./mixerCalConsts(m,d,1)));
linPhae = real(asind((squeeze(mixers(m,d,:))'-initD)./0.2204));
figure;
plot(fitX,fitPhas);
hold all;
plot(fitX,datPhase);
plot(fitX,linPhae)
figure;
plot(datPhase-fitPhas);
hold all;
plot(linPhae-fitPhas);

tmpY=squeeze(mixers(m,d,134:146))';
tmpX=134:146;
polyfit(tmpX,tmpY,1)
figure;
plot(tmpX,tmpY,'o')
hold all;
plot(tmpX,polyval(polyfit(tmpX,tmpY,1),tmpX));
plot(fitX,fitY);

%%
ampM = (max(mixers,[],3)-min(mixers,[],3))/2;
diffM = NaN(nMixers,nPowerLevels);
relDiffM = NaN(nMixers,nPowerLevels);
for m=1:nMixers
    for d=1:nPowerLevels
        diffM(m,d) = ampM(m,d)-mixerCalConsts(m,d,1);
        relDiffM(m,d) = diffM(m,d)./mixerCalConsts(m,d,1);
    end
end
figure;
plot(powerLevels,(diffM'))
figure;
plot(powerLevels,(relDiffM'))


%%
%% Calibrations - mixer
% restricting fits to within +/- 0.5 A
fprintf('Calculating Mixer calibrations...\n');

% estimate frequency (number of samples per 360 degrees)
estFreq = NaN(nMixers,nPowerLevels);
for m=1:nMixers
    for d=1:nPowerLevels
        [pks,locs] = findpeaks(squeeze(mixers(m,d,:)),'MinPeakHeight',0,'NPeaks',7,'SortStr','descend'); % returns 7 largest peaks sorted from largest to smallest. 7 picked to match no. of peaks in data.
        estFreq(m,d) = nanmean(diff(sort(locs))); % average of the difference in peak locations. Sort needed as peaks are returned in order from largest to smallest but we want difference between neighbouring peaks.
    end
end
estFreq = nanmean(removeBadPulses(estFreq(:))); % overall best estimate, mean removing any outliers.

restMixerCalConsts = NaN(nMixers,nPowerLevels,4);
mixerCalRSquare = NaN(nMixers,nPowerLevels);

for m=1:nMixers
    fprintf('Mixer %d...\n',m);

    for d=1:nPowerLevels        
        tmpMixer = squeeze(mixers(m,d,:));
        isGood = ~isnan(tmpMixer);
        tmpMixer = tmpMixer(isGood);
        tmpSamples = 1:nSamples;
        tmpSamples = tmpSamples(isGood)';
        
        initA = (maxMixers(m,d) - minMixers(m,d))./2; % est amplitude
        initB = 2*pi/estFreq; % est frequency, s = no. samples between neighbouring max and min (180 degrees)
        initD = (maxMixers(m,d) + minMixers(m,d))./2; % est vertical offset
        
        sinArg = (tmpMixer-initD)./initA;
        sinArg(sinArg>1)= NaN; % avoid imaginary component if data point above amplitude guess
        sinArg(sinArg<-1)= NaN;
        sinArg = asin(sinArg);
        sinArg(diff(tmpMixer)<0) = NaN;
        initC = sinArg-(initB.*tmpSamples);
        initC = mod(initC,2*pi); % make everything between 0 and 2pi
        initC(initC>pi) = initC(initC>pi)-2*pi;
        initC(initC<-pi) = initC(initC<-pi)+2*pi;
        initC = nanmean(initC);

        initParams = [initA, initB, initC, initD];


            % need to select 1st quadrant where asin valid - positive
            % gradient in difference
            
            initC = mod(initC,2*pi); % make everything between 0 and 2pi
            initC(initC>pi) = initC(initC>pi)-2*pi;
            initC(initC<-pi) = initC(initC<-pi)+2*pi;
            initC = nanmean(initC);

            
        [~,maxInd] = findpeaks(smooth(tmpMixer));
        [~,minInd] = findpeaks(-smooth(tmpMixer));
        fitStart = maxInd(2);
        fitEnd = find(minInd-maxInd(2)>0,1);
        fitEnd = minInd(fitEnd);
%         fitEnd = maxInd(3);
        tmpSamples = fitStart:fitEnd;
        tmpMixer = tmpMixer(tmpSamples);
%         tmpSamples = tmpSamples(abs(tmpMixer-initD)<0.75.*initA);
%         tmpMixer = tmpMixer(abs(tmpMixer-initD)<0.75.*initA);
        
        [calibrationFactors, fitRSquare, fitConfInt] = offsetSinFit(tmpSamples', tmpMixer, initParams);
        
        % correct negative amplitude fit
        restMixerCalConsts(m,d,:) = calibrationFactors;           
        mixerCalRSquare(m,d) = fitRSquare;
    end
end

figure;
m=1;
d=7;
restFitFunc = @(x) (restMixerCalConsts(m,d,1).*sin( (restMixerCalConsts(m,d,2).*x) + restMixerCalConsts(m,d,3) )) + restMixerCalConsts(m,d,4);
fitX = 1:nSamples;
restFitY = restFitFunc(fitX);
fitFunc = @(x) (mixerCalConsts(m,d,1).*sin( (mixerCalConsts(m,d,2).*x) + mixerCalConsts(m,d,3) )) + mixerCalConsts(m,d,4);
fitY = fitFunc(fitX);

plot(fitX,restFitY,'k','LineWidth',1.5);
hold all;
plot(fitX,fitY,'m','LineWidth',1.5);
plot(squeeze(mixers(m,d,:)),'o','Color',mixerColours(m,:));
title(sprintf('Mixer %d, Power %d dBm (Amplitude: %.2d V, Rsquare: %.3f)',m,powerLevels(d),abs(restMixerCalConsts(m,d,1)),restMixerCalConsts(m,d)));
xlabel('Sample No.');
ylabel('Output [V]');
grid on;
legend('Rest Fit','All Fit','Data');
hold off;

figure
plot(squeeze(mixerCalConsts(1,:,1))')
hold all
plot(squeeze(restMixerCalConsts(1,:,1))','--')
legend('all fit','restricted fit')

figure
plot(squeeze(abs(mixerCalConsts(1,:,1)-restMixerCalConsts(1,:,1))./restMixerCalConsts(1,:,1)),'b')


%%
m=1;
d=7;
fitX = 1:nSamples;
fitFunc = @(x) (restMixerCalConsts(m,d,1).*sin( (restMixerCalConsts(m,d,2).*x) + restMixerCalConsts(m,d,3) )) + restMixerCalConsts(m,d,4);
fitY = fitFunc(fitX);
figure;
plot(abs(squeeze(mixers(m,d,:))-fitY'),'Color',mixerColours(m,:));
hold all;
% legend('Mixer 1','Mixer 2','Mixer 3');
title(sprintf('Mixer Fit Residuals at %d dBm ',powerLevels(d)));
xlabel('Sample No.');
ylabel('Output [V]');
xlim([93 125])
format_plots


figure
d=7;

m=1;
fitX = 1:nSamples;
fitFunc = @(x) (restMixerCalConsts(m,d,1).*sin( (restMixerCalConsts(m,d,2).*x) + restMixerCalConsts(m,d,3) )) + restMixerCalConsts(m,d,4);
fitY = fitFunc(fitX);
plot((squeeze(mixers(m,d,:))-delaySignal(fitY,-2)'),'Color',mixerColours(m,:),'LineWidth',2);
hold all;

m=2;
fitX = 1:nSamples;
fitFunc = @(x) (restMixerCalConsts(m,d,1).*sin( (restMixerCalConsts(m,d,2).*x) + restMixerCalConsts(m,d,3) )) + restMixerCalConsts(m,d,4);
fitY = fitFunc(fitX);
plot((squeeze(mixers(m,d,:))-fitY'),'Color',mixerColours(m,:),'LineWidth',2);

m=3;
fitX = 1:nSamples;
fitFunc = @(x) (restMixerCalConsts(m,d,1).*sin( (restMixerCalConsts(m,d,2).*x) + restMixerCalConsts(m,d,3) )) + restMixerCalConsts(m,d,4);
fitY = fitFunc(fitX);
plot((squeeze(mixers(m,d,:))-fitY'),'Color',mixerColours(m,:),'LineWidth',2);
hold all;


legend('Mixer 1','Mixer 2','Mixer 3');
title(sprintf('Mixer Fit Residuals at %d dBm ',powerLevels(d)));
xlabel('Sample No.');
ylabel('Output [V]');
format_plots
hold off;

%%
figure
d=10;

m=1;
fitX = 1:nSamples;
fitFunc = @(x) (restMixerCalConsts(m,d,1).*sin( (restMixerCalConsts(m,d,2).*x) + restMixerCalConsts(m,d,3) )) + restMixerCalConsts(m,d,4);
fitY = fitFunc(fitX);
plot(asind((squeeze(mixers(m,d,93:125))-fitY(93:125)')./restMixerCalConsts(m,d,1)),'Color',mixerColours(m,:),'LineWidth',2);
hold all;

m=2;
fitX = 1:nSamples;
fitFunc = @(x) (restMixerCalConsts(m,d,1).*sin( (restMixerCalConsts(m,d,2).*x) + restMixerCalConsts(m,d,3) )) + restMixerCalConsts(m,d,4);
fitY = fitFunc(fitX);
plot(asind((squeeze(mixers(m,d,79:111))-fitY(79:111)')./restMixerCalConsts(m,d,1)),'Color',mixerColours(m,:),'LineWidth',2);

m=3;
fitX = 1:nSamples;
fitFunc = @(x) (restMixerCalConsts(m,d,1).*sin( (restMixerCalConsts(m,d,2).*x) + restMixerCalConsts(m,d,3) )) + restMixerCalConsts(m,d,4);
fitY = fitFunc(fitX);
plot(asind((squeeze(mixers(m,d,57:89))-fitY(57:89)')./restMixerCalConsts(m,d,1)),'Color',mixerColours(m,:),'LineWidth',2);
hold all;

legend('Mixer 1','Mixer 2','Mixer 3');
title(sprintf('Mixer Fit Residuals at %d dBm ',powerLevels(d)));
xlabel('Sample No.');
ylabel('Residual [degrees]');
format_plots
hold off;


%%
[x,y,n]=getDataFromFig;

m1=y{3};
m2=y{2};
m3=y{1};

fit1=polyfit(powerInVolts(diodePowerLevelsToFit),m1(diodePowerLevelsToFit),2);
fit2=polyfit(powerInVolts(diodePowerLevelsToFit),m2(diodePowerLevelsToFit),2);
fit3=polyfit(powerInVolts(diodePowerLevelsToFit),m3(diodePowerLevelsToFit),2);

fitX = linspace(powerInVolts(1),powerInVolts(end),1000);
hold all;
plot(fitX,polyval(fit1,fitX),'--b','LineWidth',1.5);
plot(fitX,polyval(fit2,fitX),'--r','LineWidth',1.5);
plot(fitX,polyval(fit3,fitX),'--g','LineWidth',1.5);
legend('Mixer 1 (data)','Mixer 2 (data)','Mixer 3 (data)','Mixer 1 (quadratic)','Mixer 2 (quadratic)','Mixer 3 (quadratic)')
xlim([0.2890 1.2842])