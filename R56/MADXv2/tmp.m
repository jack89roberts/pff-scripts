load('/home/jack/Documents/MATLAB/PhaseFeedforward/ctfmod/PhaseFeedforward/R56/MADXv2/t566.mat');
phase = myResults.opticsPhase;
deltaP = myResults.deltaP;
r56Vals = myResults.opticsR56;
r56TL1 = myResults.r56Values;
t566Vals =  myResults.opticsT566;

diffP = NaN(length(r56Vals),length(deltaP));
for r=1:length(r56Vals)
    diffP(r,1) = phase(r,2)-phase(r,1);
end
for r=1:length(r56Vals)
    for p=2:(length(deltaP)-1)
        diffP(r,p) = phase(r,p+1)-phase(r,p-1);
    end
end
for r=1:length(r56Vals)
    diffP(r,end) = phase(r,end)-phase(r,end-1);
end


[~,bestR56Inds] = min(abs(diffP));
bestR56 = r56Vals(bestR56Inds);

figure;
plot(deltaP,bestR56+(r56TL1(1)-r56Vals(1)),'o','MarkerFaceColor','b','LineWidth',2)
hold all;
plot(deltaP,(-2*t566Vals(bestR56Inds).*deltaP)+(r56TL1(1)-r56Vals(1)),'k','LineWidth',2);
xlabel('\Deltap / p')
ylabel('R56 in TL1 [m]')
ylim([-0.35 0.65])
title('Simultated Optimal TL1 Optics vs. Energy Offset')
legend('MADX','Equation')
format_plots