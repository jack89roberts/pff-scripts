close all

gain=624;
zeroSamp = 421;

phas = squeeze(nanmean(FONTData.ADCs,2));
phas = phas(2,:);
phas = averageSamples(phas,5);
phas = phas-phas(zeroSamp);


offUp = 200;
offDown = -200;

phasOffUp = phas+offUp;
phasOffDown = phas+offDown;


dacOpt = -phas*(gain/64);
dacUp = -phasOffUp*(gain/64);
dacDown = -phasOffDown*(gain/64);
dacOpt(dacOpt>4095) = 4095;
dacOpt(dacOpt<-4096) = -4096;
dacUp(dacUp>4095) = 4095;
dacUp(dacUp<-4096) = -4096;
dacDown(dacDown>4095) = 4095;
dacDown(dacDown<-4096) = -4096;


myx = [212 616];

figure;
plot(phas,'b','LineWidth',2);
hold all;
plot(phasOffUp,'r','LineWidth',2);
% plot(phasOffDown);
xlim(myx);
ylim([-500 3500])
plot([zeroSamp zeroSamp],get(gca,'YLim'),'k--');
xlabel('Sample No.')
ylabel('ADC2 Input [counts]')
format_plots;
grid off;

figure;
plot(dacOpt,'b','LineWidth',2);
hold all;
plot(dacUp,'r','LineWidth',2);
% plot(dacDown);
xlim(myx);
ylim([-4500 500])
plot([zeroSamp zeroSamp],get(gca,'YLim'),'k--');
xlabel('Sample No.')
ylabel('DAC Output [counts]')
format_plots;
grid off;