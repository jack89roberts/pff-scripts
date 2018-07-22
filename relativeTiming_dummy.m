avg=0;
sigma=1;
ampFX2 = 0.8;

x_st=-3:0.01:0;
fx_st=1/sqrt(2*pi)/sigma*exp(-(x_st-avg).^2/2/sigma/sigma);
x_en=0:0.01:3;
fx_en=1/sqrt(2*pi)/sigma*exp(-(x_en-avg).^2/2/sigma/sigma);
fx_mid = max(fx_st).*ones(1,3*length(fx_st));
fx1 = [zeros(1,500) fx_st fx_mid fx_en zeros(1,500)];

x_st=-3:0.01:0;
fx_st=1/sqrt(2*pi)/sigma*exp(-(x_st-avg).^2/2/sigma/sigma);
x_en=0:0.01:3;
fx_en=1/sqrt(2*pi)/sigma*exp(-(x_en-avg).^2/2/sigma/sigma);
fx_mid = max(fx_st).*ones(1,3*length(fx_st));
fx2 = -[zeros(1,500) ampFX2*fx_st ampFX2*fx_mid ampFX2*fx_en zeros(1,500)];

delay = -200:50:200;
nDelays = length(delay);

delayFX2 = NaN(nDelays,length(fx1));
diffFX2 = NaN(nDelays,length(fx1));
for i=1:nDelays
    delayFX2(i,:) = delaySignal(fx2,delay(i));
    diffFX2(i,:) = fx1+delayFX2(i,:);
end

figure;
plot(max(diffFX2(:,300:1000),[],2)-min(diffFX2(:,300:1000),[],2),'o');
hold all;
plot(max(diffFX2(:,1500:2200),[],2)-min(diffFX2(:,1500:2200),[],2),'o');
figure;
plot(diffFX2');

firstSlope =max(diffFX2(:,300:1000),[],2)-min(diffFX2(:,300:1000),[],2);
firstSlopeA = firstSlope(1:5);
firstSlopeB = firstSlope(5:end);

figure;
plot(fx1,'b','LineWidth',2);
hold all;
plot(delayFX2(end,:),'r','LineWidth',2);
plot(diffFX2(end,:),'k','LineWidth',2);
xlabel('Time [a.u.]')
ylabel('Residual Kick [a.u.]')
legend('K1','K2','SUM')
title('Simulated Response to Offset Kicks')
xlim([200 2500])
plot([788 788],[-0.4 0.4],'b--','LineWidth',1);
plot([1719 1719],[-0.4 0.4],'b--','LineWidth',1);
plot([987 987],[-0.4 0.4],'r--','LineWidth',1);
plot([1926 1926],[-0.4 0.4],'r--','LineWidth',1);
format_plots;
savePlot(saveDir,'simulation')

% hold all;
% plot(fx1,'b');
% plot(delayFX2','r')
% figure;
% plot(fx1)
% hold all
% plot(fx2)
% plot(fx1+fx2)
% 
