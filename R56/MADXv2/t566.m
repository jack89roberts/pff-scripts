myPath ='/home/jack/Documents/MATLAB/PhaseFeedforward/ctfmod/PhaseFeedforward/R56/MADXv2';
addpath('/home/jack/Documents/MATLAB/madx2Matlab');
addpath('../');

cd(myPath);

% Create model
machineModel = madx2matlab(fullfile(myPath,'generalScript.madx'),...
    '',...
    fullfile('/usr/local/bin/madx'));
...%fullfile(myPath,'./madx'));


% assign a twiss command line not as the default one
machineModel.twissCommandLine=[...
    'select, flag=twiss, clear;',char(10),...
    'select, flag=twiss, column=name, s, BETX, ALFX, L, ANGLE, x, px, dx, BETY, ALFY, y, py, dy,mux,muy, t, pt,re56, t566;',char(10),...
    'twiss, RMATRIX,',char(10),...
    'BETX=initialBETX,ALFX=initialALFX,',char(10),...
    'BETY=initialBETY,ALFY=initialALFY,',char(10),...
    'X=0,Y=0,',char(10),...
    'DELTAP=initialDeltaP;',char(10),...
    ];

% sets of TL1 optics
[quadNames, r56Values, quadCurrents] = loadR56CurrentsTable(140);
quadNames = strrep(quadNames,'.','.I');

deltaP = -0.02:0.0025:0.02;


% initial values for twiss
initialBETX =    10;
initialBETY =    1.5;
initialALFX =    0;
initialALFY =    0;
initialX = 0;
initialY = 0;
initialDeltaP = 0.0;

%%

opticsR56 = NaN(1,length(r56Values));
opticsT566 = NaN(1,length(r56Values));

for r=1:length(r56Values)
    fprintf('-------------------------------\n');
    fprintf('R56 = %.3f...\n',r56Values(r));
    fprintf('-------------------------------\n');
    
    machineModel.setModelParameters(quadNames,quadCurrents(r,:));

    % set initial settings and compute optics
    machineModel.setModelParameters('initialBETX', initialBETX);
    machineModel.setModelParameters('initialBETY', initialBETY);
    machineModel.setModelParameters('initialALFX', initialALFX);
    machineModel.setModelParameters('initialALFY', initialALFY);
    machineModel.setModelParameters('initialX', initialX);
    machineModel.setModelParameters('initialY', initialY);
    machineModel.setModelParameters('initialDELTAP', initialDeltaP);

    nominalOptic = machineModel.computeOptics();

    if (r==1)
        opticsPhase = NaN(length(r56Values),length(deltaP));
    end
    
    finalT = NaN(length(deltaP),1);

    for p=1:length(initialDeltaPs)
        fprintf('deltaP/P = %.3f...\n',deltaP(p));
        machineModel.setModelParameters('initialDELTAP', deltaP(p));
        auxOptics = machineModel.computeOptics();

        finalT(p) = auxOptics.DATA(end).T;
    end
    
    % 1st order effect not included in normal twiss output, must be added
    % afterwards.
    finalT = finalT' + deltaP.*nominalOptic.DATA(end).RE56;
    opticsPhase(r,:) = finalT;
    
    % fit
    tFit = polyfit(initialDeltaPs,finalT,2);
    opticsR56(r) = tFit(2);
    opticsT566(r) = tFit(1);

end

%%
figure;
plot(r56Values,opticsR56,'LineWidth',2);
xlabel('R56 in TL1 [m]')
ylabel('R56 Upstream-Downstream [m]')
title('R56 Between Phase Monitors vs. R56 in TL1');
format_plots;

figure;
plot(r56Values,opticsT566);
xlabel('R56 in TL1 [m]')
ylabel('T566 Upstream-Downstream [m]')
title('T566 Between Phase Monitors vs. R56 in TL1');
format_plots;

myCols = varycolor(length(r56Values));
figure;
for i=1:length(r56Values)
    plot(deltaP,opticsPhase(i,:)*(360/0.025),'Color',myCols(i,:));
    hold all;
end
set(gcf, 'Colormap', myCols);
figColBar = colorbar;
allTicks = linspace(0,1,length(r56Values)+1);
set(figColBar,'YTick',allTicks(1:10:length(r56Values)));
set(figColBar,'YTickLabel',r56Values(1:10:end));
ylabel(figColBar,'R56 in TL1 [m]');
xlabel('\Delta p / p')
ylabel('Phase [degrees]')
title('Phase vs. Beam Energy: Including Higher Orders')
format_plots;

figure;
for i=1:length(r56Values)
    plot(deltaP,opticsR56(i).*deltaP.*(360/0.025),'Color',myCols(i,:))
    hold all
end
set(gcf, 'Colormap', myCols);
figColBar = colorbar;
allTicks = linspace(0,1,length(r56Values)+1);
set(figColBar,'YTick',allTicks(1:10:length(r56Values)));
set(figColBar,'YTickLabel',r56Values(1:10:end));
ylabel(figColBar,'R56 in TL1 [m]');
xlabel('\Delta p / p')
ylabel('Phase [degrees]')
title('Phase vs. Beam Energy: R56 only')
format_plots;

figure;
plot(deltaP,abs(deltaP.*opticsR56(61).*(360/0.025)),'LineWidth',2);
hold all;
plot(deltaP,abs((deltaP.^2).*opticsT566(61).*(360/0.025)),'LineWidth',2)
xlabel('\Deltap / p')
ylabel('Phase [degrees]')
title({'Energy Dependent Phase Shift' 'R56 = 0m in TL1'})
legend('R56','T566')
format_plots

myResults = struct();
myResults.opticsR56 = opticsR56;
myResults.opticsT566 = opticsT566;
myResults.opticsPhase = opticsPhase;
myResults.deltaP = deltaP;
myResults.r56Values = r56Values;
save('t566.mat','myResults')
