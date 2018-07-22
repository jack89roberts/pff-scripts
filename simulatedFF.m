%%
saveDir = '/home/jack/PhaseFeedforward/Analysis/Simulation/EffectOfLimitedCorrection';
addpath('/home/jack/Documents/MATLAB/PhaseFeedforward/ctfmod/PhaseFeedforward/Analysis');
addpath('/home/jack/Documents/MATLAB/PhaseFeedforward/ctfmod/PhaseFeedforward');

%% make correlated distribution for upstream and downstream phase
jitterUp = 1;
jitterDown = 1;
correl = 0.95;
nPoints = 1000000;

phaseUp = random('norm',0,jitterUp,1,nPoints);
phaseDown = random('norm',0,jitterUp,1,nPoints);
phaseDown = correl*phaseUp + sqrt(1-(correl^2)).*phaseDown;
jitterRatio = jitterDown/jitterUp;
phaseDown = phaseDown*jitterRatio;

distJitterUp = std(phaseUp)
distJitterDown = std(phaseDown)
distCorrel = corrcoef(phaseUp,phaseDown);
distCorrel = distCorrel(1,2)
eqnCorrJitter = distJitterDown.*sqrt(1-(distCorrel.^2))

%% phase along pulse - just quadratic shape added to distribution above

phaseAlongPulse = (-(-50:50).^2)*0.002 + 2;

pulsePhasesUp = NaN(nPoints,length(phaseAlongPulse));
pulsePhasesDown = NaN(nPoints,length(phaseAlongPulse));
for i=1:nPoints
    pulsePhasesUp(i,:) = phaseAlongPulse + phaseUp(i);
    pulsePhasesDown(i,:) = phaseAlongPulse + phaseDown(i);
end

%% FF results with different limits on max correction
corrLimitVals = 0:0.05:3;%0:0.01:4;
%gainMultipliers = 0.5:0.01:1.5;

gain = (distCorrel*distJitterDown)/distJitterUp;

fracOutLimit = NaN(1,length(corrLimitVals));
distCorrJitter = NaN(1,length(corrLimitVals));
distCorrCorrel = NaN(1,length(corrLimitVals));
% fracOutLimit = NaN(1,length(gainMultipliers));
% distCorrJitter = NaN(1,length(gainMultipliers));
% distCorrCorrel = NaN(1,length(gainMultipliers));
figure;
load('/home/jack/Documents/MATLAB/cloudPlot/myColorMap.mat');
for i=1:length(corrLimitVals)
    corrLimit = corrLimitVals(i);
% for i=1:length(gainMultipliers) 
%     corrLimit = corrLimitVals(1);
%    gain = (distCorrel*distJitterDown*gainMultipliers(i))/distJitterUp;

    % calculate corrected phase
    subtractPhase = gain*phaseUp;
    fracOutLimit(i) = sum(subtractPhase>corrLimit | subtractPhase<-corrLimit)/nPoints;
    subtractPhase(subtractPhase>corrLimit) = corrLimit;
    subtractPhase(subtractPhase<-corrLimit) = -corrLimit;
    corrPhase = phaseDown - subtractPhase;

    distCorrJitter(i) = std(corrPhase);
    tmpCorr = corrcoef(phaseUp,corrPhase);
    distCorrCorrel(i) = tmpCorr(1,2);
    
    cloudPlot(phaseUp,corrPhase,[],[],[250 250])
    title(sprintf('Max Correction = %.2f^o, Correlation = %.2f, Downstream Jitter = %.2f^o',corrLimitVals(i),distCorrCorrel(i),distCorrJitter(i)));
    xlabel('Upstream Phase [degrees]');
    ylabel('Corrected Downstream Phase [degrees]');
    xlim([-4 4])
    ylim([-4 4])
    colormap(myColorMap)
    %caxis([0 30])
    format_plots;
    savePlot(saveDir,sprintf('scatterLim%.2f',corrLimitVals(i)));

end

% figure;
% plot(corrLimitVals,fracOutLimit,'LineWidth',2);
% title('Fraction of Points Outside Limit')
% ylabel('Fraction')
% xlabel('Maximum Correction [degrees]');
% format_plots;
% savePlot(saveDir,'fracPointsVsMaxCorr');
% 
% figure;
% plot(corrLimitVals,distCorrJitter,'LineWidth',2);
% title('Corrected Downstream Jitter')
% ylabel('Jitter [degrees]')
% xlabel('Maximum Correction [degrees]');
% format_plots;
% savePlot(saveDir,'jitterVsMaxCorr');
% 
% figure;
% plot(corrLimitVals,distCorrCorrel,'LineWidth',2);
% title('Correlation: Upstream Phase and Corrected Downstream Phase')
% ylabel('Correlation')
% xlabel('Maximum Correction [degrees]');
% format_plots;
% savePlot(saveDir,'correlationVsMaxCorr');

%% FF results with different gains
corrLimit = 1.5;
gainMultipliers = 0:0.1:3;

optGain = (distCorrel*distJitterDown)/distJitterUp;
gainValues = optGain.*gainMultipliers;
fracOutLimit = NaN(1,length(gainMultipliers));
distCorrJitter = NaN(1,length(gainMultipliers));
distCorrCorrel = NaN(1,length(gainMultipliers));
figure;
load('/home/jack/Documents/MATLAB/cloudPlot/myColorMap.mat');
for i=1:length(gainMultipliers) 
    gain = gainValues(i);

    % calculate corrected phase
    subtractPhase = gain*phaseUp;
    fracOutLimit(i) = sum(subtractPhase>corrLimit | subtractPhase<-corrLimit)/nPoints;
    
    subtractPhase(subtractPhase>corrLimit) = corrLimit;
    subtractPhase(subtractPhase<-corrLimit) = -corrLimit;
    corrPhase = phaseDown - subtractPhase;

    distCorrJitter(i) = std(corrPhase);
    tmpCorr = corrcoef(phaseUp,corrPhase);
    distCorrCorrel(i) = tmpCorr(1,2);
    
    cloudPlot(phaseUp,corrPhase,[],[],[250 250])
    title(sprintf('Gain = %.2f, Correlation = %.2f, Downstream Jitter = %.2f^o',gain,distCorrCorrel(i),distCorrJitter(i)));
    xlabel('Upstream Phase [degrees]');
    ylabel('Corrected Downstream Phase [degrees]');
    xlim([-4 4])
    ylim([-4 4])
    hold all;
    plot([corrLimit/gain corrLimit/gain],[-4 4],'k','LineWidth',2);
    plot([-corrLimit/gain -corrLimit/gain],[-4 4],'k','LineWidth',2);
    %caxis([0 30])
    colormap(myColorMap)
    format_plots;
    savePlot(saveDir,sprintf('scatterLim%.2fGain%.2f',corrLimit,gain));
    hold off;
end

figure;
plot(gainValues,fracOutLimit,'LineWidth',2);
hold all;
format_plots;
plot([optGain optGain],get(gca,'YLim'),'k--','LineWidth',2)
title('Fraction of Points Outside Limit')
ylabel('Fraction')
xlabel('Gain');
savePlot(saveDir,'fracPointsVsMaxCorr');

figure;
plot(gainValues,distCorrJitter,'LineWidth',2);
hold all;
title('Corrected Downstream Jitter')
ylabel('Jitter [degrees]')
xlabel('Gain');
format_plots;
plot([optGain optGain],get(gca,'YLim'),'k--','LineWidth',2)
savePlot(saveDir,'jitterVsMaxCorr');

figure;
plot(gainValues,distCorrCorrel,'LineWidth',2);
hold all;
format_plots;
plot([optGain optGain],get(gca,'YLim'),'k--','LineWidth',2)
title('Correlation: Upstream Phase and Corrected Downstream Phase')
ylabel('Correlation')
xlabel('Gain');
savePlot(saveDir,'correlationVsMaxCorr');

%% FF results ALONG PULSE with different limits on max correction
corrLimit = 2.5;
pulsePhasesCorr = NaN(nPoints,length(phaseAlongPulse));
gain = (distCorrel*distJitterDown)/distJitterUp;

% calculate corrected phase
for i=1:nPoints
    subtractPhase = gain*pulsePhasesUp(i,:);
    fracOutLimit = sum(subtractPhase>corrLimit | subtractPhase<-corrLimit)/nPoints;
    subtractPhase(subtractPhase>corrLimit) = corrLimit;
    subtractPhase(subtractPhase<-corrLimit) = -corrLimit;
    pulsePhasesCorr(i,:) = pulsePhasesDown(i,:) - subtractPhase;
end

figure;
plot(mean(pulsePhasesUp));
hold all;
plot(mean(pulsePhasesDown))
plot(mean(pulsePhasesCorr))
legend('up','down','corr')
title('mean along pulse')

figure;
plot(std(pulsePhasesUp));
hold all;
plot(std(pulsePhasesDown))
plot(std(pulsePhasesCorr))
legend('up','down','corr')
title('std along pulse')

%% equation: jitter vs correlation
% 
% jitterUp = 1;
% rho = 0:0.001:1;
% 
% jitterOptiGain = jitterUp.*sqrt(1 - (rho.^2)); % jitter after ff if input is correlation*upstream phase
% jitterFullGain = sqrt(2).*jitterUp.*sqrt(1-rho); % jitter after ff if input is upstream phase
% 
% figure;
% plot(rho,jitterFullGain);
% hold all;
% plot(rho,jitterOptiGain);
% legend('PhiD - PhiU','PhiD - a*PhiU')
% xlabel('CORRELATION')
% ylabel('JITTER')
