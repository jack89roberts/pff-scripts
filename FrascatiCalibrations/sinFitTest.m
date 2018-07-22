close all;
%%
x = 0:0.2:2*pi;
y = 1*sin(1*x + 1) +1 + random('norm',0,0.1,1,length(x));
z = cos(x);
% figure;
% plot(x,y,'o');
%%

% sensible guesses for fit parameters
initA = (max(y)-min(y))/2 % amplitude
initB = 1 % phase units, 4x for 3GHz -> 12GHz
initD = (max(y)+min(y))/2 % vertical offset
% calculate phase offset using other guesses for initial
% parameters and phase/output at first point.
sinArg = (y-initD)./initA;
sinArg(sinArg>1)= NaN; % avoid imaginary component if data point above amplitude guess
sinArg(sinArg<-1)= NaN;
sinArg = asin(sinArg);
% need to select 1st quadrant where asin valid - positive
% gradient in difference
sinArg(diff(y)<0) = NaN;
initC = (sinArg-(initB.*x));
initC = mod(initC,2*pi);
initC(initC>pi) = initC(initC>pi)-2*pi;
initC(initC<-pi) = initC(initC<-pi)+2*pi;
initC = nanmean(initC)
initParams = [initA initB initC initD];
%%
[fitConst, fitRSquare, fitConfInt]=offsetSinFit(x',y',initParams);
fitConst
fitConfInt/2

figure;
plot(x,y,'o');
hold all;
fitY = fitConst(1)*sin(fitConst(2).*x + fitConst(3)) + fitConst(4);
plot(x, fitY,'r','LineWidth',2);
initY = initA*sin(initB.*x + initC) + initD;
plot(x, initY,'r--');
xlabel('x');
ylabel('y');
grid on;
legend('Data','Fit','Initial Guess');
