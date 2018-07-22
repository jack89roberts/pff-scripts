%% Noisy Signal
% Use Fourier transforms to find the frequency components of a signal buried 
% in noise.
% 
% Specify the parameters of a signal with a sampling frequency of 1 kHz and
% a signal duration of 1 second.

% Copyright 2015 The MathWorks, Inc.

datFile = '/home/jack/PhaseFeedforward/Analysis/201612/20161215_2340_Straight_Gain_1200_K1Del_7_151216_FFOff.mat';

datLabels = {
    'Upstream'
    'Downstream'
};
datFormats = {
   '-^'
   '-v'
};
datCols = [0 0.8 0; 0 0 1];

monInds = [1 3];
myRange = 600:755;
Fs = 357e6;            % Sampling frequency 

%%
T = 1/Fs;             % Sampling period       
L = length(myRange);             % Length of signal
t = (0:L-1)*T;        % Time vector

%%
load(datFile);
lineRefs = cell(1,length(monInds));
figure;

for i=1:length(monInds)

    X = squeeze(phases(monInds(i),:,myRange));
    for p=1:nPulses % subtract average of each pulse separately
        X(p,:) = X(p,:)-nanmean(X(p,:));
    end
    X = X';

    %% 
    % Compute the Fourier transform of the signal. 

    Y = fft(X);
    % Y = nanmean(Y,2);
    %%
    % Compute the two-sided spectrum |P2|.  Then compute the single-sided
    % spectrum |P1| based on |P2| and the even-valued signal length |L|.

    P2 = abs(Y/L);
    P1 = P2(1:L/2+1,:);
    P1(2:end-1,:) = 2*P1(2:end-1,:);

    %% Calculate mean and std spectrum
    addpath('/home/jack/cernbox/DocsLinux/MATLAB/PhaseFeedforward/ctfmod/PhaseFeedforward');
    [meanP1,stdP1] = nanMeanStdErr(P1,2);

    %% 
    % Define the frequency domain |f| and plot the single-sided amplitude
    % spectrum |P1|.  The amplitudes are not exactly at 0.7 and 1, as expected, because of the added 
    % noise. On average, longer signals produce better frequency approximations.

    f = Fs*(0:(L/2))/L;
    f = f/1e6;

    %lineRefs{i}=shadedErrorBar(f,meanP1,stdP1,{datFormats{i},'Color',datCols{i},'MarkerFaceColor',datCols{i}});
    plot(f,meanP1,datFormats{i},'Color',datCols(i,:),'MarkerFaceColor',datCols(i,:));
    hold all;
end
%%
title('FFT of Upstream and Downstream Phase')
xlabel('f [MHz]')
ylabel('Amplitude [a.u.]')
%legend([lineRefs{1}.mainLine lineRefs{2}.mainLine],datLabels);
legend(datLabels);
xlim([0 40])