%% errorFromFONTPhaseReconstruction
% Script to look at the effect that using the phase=Mixer/A approximation
% in the feedforward algorithm has on the measured upstream, downstream and
% corrected phase jitters.

% "FONT" method: as used in PFF algorithm (Mixer/A).
% "SiS" method: as used in data analysis (asin((Mixer-d)/A))

% Resolution of monitors not considered (i.e. assumes resolution = zero).

%% Inputs: 
% jitters/correlations and monitor calibration constants in second section.

jitterUp = 0.8; %0.9 % upstream phase jitter (degrees)
jitterDown = 0.8; %1.4 % downstream phase jitter
offsets = -80:80;% mean offset in upstream phase (mean upstream phase = offset, mean downstream phase = jitRatio*offset)
correl = 0.97; %0.89 % upstream-downstream phase correlation
nPoints = 100000; % no. points in generated distribution

A_up = 1.05; % upstream calibration amplitude
d_up = 0.046; % upstream calibration min-max offset
A_down = 1.013; % downstream calibration amplitude
d_down = 0.016; % downstream calibration min-max offset

%%
addpath('../Analysis/');
addpath('../');
%% basic difference across full 360 degrees phase range

phase = 0:0.01:360;

% convert generated phases in distribution to equivalent mixer outputs
Mixer = A_up*sind(phase) + d_up;

% reconstruct measured phase using FONT/SiS method
phaseFONT = (Mixer./A_up)*(180/pi);
phaseSiS = asind((Mixer-d_up)./A_up);

figure;
plot(phase-180,phaseFONT,'LineWidth',3);
hold all;
plot(phase-180,phaseSiS,'LineWidth',3);
plot(phase-180,phaseFONT-phaseSiS,'LineWidth',3);
xlabel('Real Beam Phase [degrees]')
ylabel('Measured Phase [degrees]')
legend('FONT','SiS','Diff')
title('Difference Between Phase Reconstruction Methods')
format_plots;
% savePlot('/home/jack/PhaseFeedforward/Analysis/Simulation/FONTPhaseReconstruction','fullResponse_360degrees');
xlim([-10 10]);
% savePlot('/home/jack/PhaseFeedforward/Analysis/Simulation/FONTPhaseReconstruction','fullResponse_pm10degrees');
xlim([-90 90]);
% savePlot('/home/jack/PhaseFeedforward/Analysis/Simulation/FONTPhaseReconstruction','fullResponse_pm90degrees');

%% make correlated distribution for upstream and downstream phase

nOffsets = length(offsets);
stdCorrFONT = NaN(1,nOffsets);
stdCorrSiS = NaN(1,nOffsets);
for i=1:nOffsets
    
% make the correlated distributions
phaseUp = random('norm',0,jitterUp,1,nPoints);

phaseDown = random('norm',0,jitterUp,1,nPoints);
phaseDown = correl*phaseUp + sqrt(1-(correl^2)).*phaseDown;
jitterRatio = jitterDown/jitterUp;
phaseDown = phaseDown*jitterRatio;

% add the mean offset
offset = offsets(i);
phaseUp = phaseUp+offset;
phaseDown = phaseDown+(offset*jitterRatio);

% phase = 0:0.01:2*pi;
% phase=random('norm',0,1.4,1,100000)+30;
% phase = phase*(pi/180);
% figure;
% histogram(phase);

% some basic statistics from generated distribution
distJitterUp = std(phaseUp);
distJitterDown = std(phaseDown);
distCorrel = nancorrcoef(phaseUp,phaseDown);

optGain = (distJitterDown/distJitterUp)*distCorrel;


%% convert phase in to mixer output

% convert generated phases in distribution to equivalent mixer outputs
MixerUp = A_up*sind(phaseUp) + d_up;
MixerDown = A_down*sind(phaseDown) + d_down;

% reconstruct measured phase using FONT/SiS method
phaseUpFONT = (MixerUp./A_up)*(180/pi);
phaseUpSiS = asind((MixerUp-d_up)./A_up);

stdUpFONT = std(phaseUpFONT); % upstream jtiter on FONT measurement
stdUpSiS = std(phaseUpSiS); % upstream jitter on SiS measurement
stdDiffUpFONT = std(phaseUpFONT-phaseUp); % jitter of difference between FONT and real phase
stdDiffUpSiS = std(phaseUpSiS-phaseUp); % jitter of difference between SiS and real phase (0)

% calculate optimal corrected phase using FONT and SiS measurements
realCorrPhaseFONT = phaseDown - optGain*phaseUpFONT;
realCorrPhaseSiS  = phaseDown - optGain*phaseUpSiS;

% convert the real corrected phases in to the measured phase in both
% methods
mixerCorrPhaseFONT = A_down*sind(realCorrPhaseFONT) + d_down;
mixerCorrPhaseSiS = A_down*sind(realCorrPhaseSiS) + d_down;

measCorrPhaseFONT = (mixerCorrPhaseFONT./A_down)*(180/pi);
measCorrPhaseSiS = asind((mixerCorrPhaseSiS-d_down)./A_down);

% some statistics
stdRealCorrPhaseFONT = std(realCorrPhaseFONT);
stdRealCorrPhaseSiS = std(realCorrPhaseSiS);
stdCorrFONT(i) = stdRealCorrPhaseFONT;
stdCorrSiS(i) = stdRealCorrPhaseSiS;

stdMeasCorrPhaseFONT = std(measCorrPhaseFONT);
stdMeasCorrPhaseSiS = std(measCorrPhaseSiS);

end

%%
figure;
plot(offsets,stdCorrFONT,'LineWidth',2)
hold all;
plot(offsets,stdCorrSiS,'LineWidth',2)
legend('FONT','SiS')
xlabel('Mean Phase Offset [degrees]')
ylabel('Downstream Jitter [degrees]')
title(sprintf('init up jit = %.2f^o, init down jit = %.2f^o, correl =%.3f',jitterUp,jitterDown,correl));
format_plots;
xlim([-80 80])