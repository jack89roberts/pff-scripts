% to run from my laptop
%myPath='/Users/davideg/CERN/work/scripts/particleTracking/testPhaseAdvance';
%addpath(fullfile(myPath,'../../../madx2MatlabLib'));

% to run on CERN technical network
%myPath='/nfs/cs-ccr-nfs6/vol29/share/ctf/data/ctfmod/MatLab/PhaseFeedforward/R56/MADX';
myPath ='/home/jack/Documents/MATLAB/PhaseFeedforward/ctfmod/PhaseFeedforward/R56/MADXv2';
%addpath('/acc/oper/Linux/mcr/matlab/madx2Matlab');
addpath('/home/jack/Documents/MATLAB/madx2Matlab');
addpath('../');

cd(myPath);

% Create model
machineModel = madx2matlab(fullfile(myPath,'generalScript.madx'),...
    '',...
    fullfile('/usr/local/bin/madx'));
...%fullfile(myPath,'./madx'));


%% assign a twiss command line not as the default one
machineModel.twissCommandLine=[...
    'select, flag=twiss, clear;',char(10),...
    'select, flag=twiss, column=name, s, BETX, ALFX, L, ANGLE, x, px, dx, BETY, ALFY, y, py, dy,mux,muy, t, pt,re56, t566;',char(10),...
    'twiss, RMATRIX,',char(10),...
    'BETX=initialBETX,ALFX=initialALFX,',char(10),...
    'BETY=initialBETY,ALFY=initialALFY,',char(10),...
    'X=0,Y=0,',char(10),...
    'DELTAP=initialDeltaP;',char(10),...
    ];


%% initial values for twiss
initialBETX =    10;
initialBETY =    1.5;
initialALFX =    0;
initialALFY =    0;
initialX = 0;
initialY = 0;
initialDeltaP = 0.0;


% set initial settings and compute optics
machineModel.setModelParameters('initialBETX', initialBETX);
machineModel.setModelParameters('initialBETY', initialBETY);
machineModel.setModelParameters('initialALFX', initialALFX);
machineModel.setModelParameters('initialALFY', initialALFY);
machineModel.setModelParameters('initialX', initialX);
machineModel.setModelParameters('initialY', initialY);
machineModel.setModelParameters('initialDELTAP', initialDeltaP);

nominalOptic = machineModel.computeOptics();


% plot phase advance
subplot(1,2,1)
plot([nominalOptic.DATA.S], [nominalOptic.DATA.DX])
subplot(1,2,2)
plot([nominalOptic.DATA.S], [nominalOptic.DATA.RE56]);


%% Nominal currents, energy dependence
initialDeltaPs = -0.02:0.0025:0.02;

finalT = NaN(length(initialDeltaPs),1);

for i=1:length(initialDeltaPs)
    machineModel.setModelParameters('initialBETX', initialBETX);
    machineModel.setModelParameters('initialBETY', initialBETY);
    machineModel.setModelParameters('initialALFX', initialALFX);
    machineModel.setModelParameters('initialALFY', initialALFY);
    machineModel.setModelParameters('initialX', 0);
    machineModel.setModelParameters('initialY', 0);
    machineModel.setModelParameters('initialDELTAP', initialDeltaPs(i));
    auxOptics = machineModel.computeOptics();

    finalT(i) = auxOptics.DATA(end).T;
end

% 1st order effect not included in normal twiss output, must be added
% afterwards.
finalT = finalT' + initialDeltaPs.*nominalOptic.DATA(end).RE56;

% fit
tFit = polyfit(initialDeltaPs,finalT,2);
tFitVals = polyval(tFit,initialDeltaPs);
optimalDeltaP = -tFit(2)./(2*tFit(1));

%%
%plot(initialDeltaPs, finalT*(360/0.025));
figure;
plot(initialDeltaPs, (finalT/0.025)*360,'LineWidth',2);
title(['from ', nominalOptic.DATA(3).NAME , ' to ', nominalOptic.DATA(end-2).NAME, ' (S = ',num2str(nominalOptic.DATA(end-2).S),' m)' ])
xlabel('initial \Delta p / p_0')
%ylabel('T [12GHz degrees]')
ylabel('Phase [degrees]')
grid
hold all;
plot(initialDeltaPs,(tFitVals/0.025)*360,'LineWidth',2);
legend(...
    sprintf('TWISS (R56=%.2f)',nominalOptic.DATA(end).RE56),...
    sprintf('FIT (R56=%.2f, T566=%.2f)',tFit(2),tFit(1)),...
    'Location','SouthWest');
plot([optimalDeltaP optimalDeltaP],get(gca,'YLim'),'k');
title('MADX: Optics Phase Dependence on Energy')
format_plots;
% savePlot('/home/jack/PhaseFeedforward/Analysis/201512/20151221_FONT','phaseVsEnergy');

%% Nominal currents, X dependence
initialDeltaXs = -0.001:0.0001:0.001;

finalTX = NaN(length(initialDeltaXs),1);

for i=1:length(initialDeltaXs*(360/0.025))
    machineModel.setModelParameters('initialDELTAP', 0);
    machineModel.setModelParameters('initialX', initialDeltaXs(i));
    machineModel.setModelParameters('initialY', 0);
    auxOptics = machineModel.computeOptics();

    finalTX(i) = auxOptics.DATA(end).T;
end
%%
plot(initialDeltaXs, finalTX*(360/0.025));
title(['from ', nominalOptic.DATA(3).NAME , ' to ', nominalOptic.DATA(end-2).NAME, ' (S = ',num2str(nominalOptic.DATA(end-2).S),' m)' ])
xlabel('initial x_0 [m]')
ylabel('T [12GHz degrees]')
grid

%% Nominal currents, Y dependence
initialDeltaYs = -0.001:0.0001:0.001;

finalTY = NaN(length(initialDeltaYs),1);

for i=1:length(initialDeltaYs)
    machineModel.setModelParameters('initialDELTAP', 0);
    machineModel.setModelParameters('initialX', 0);
    machineModel.setModelParameters('initialY', initialDeltaYs(i));
    auxOptics = machineModel.computeOptics();

    finalTY(i) = auxOptics.DATA(end).T;
end
%%
plot(initialDeltaYs, finalTY*(360/0.025));
title(['from ', nominalOptic.DATA(3).NAME , ' to ', nominalOptic.DATA(end-2).NAME, ' (S = ',num2str(nominalOptic.DATA(end-2).S),' m)' ])
xlabel('initial y_0 [m]')
ylabel('T [12GHz degrees]')
grid

%%
plot([nominalOptic.DATA.S],[nominalOptic.DATA.DX])

machineModel.generalExecMadxScript('','value, cr.iqdh0840;')

%% R56
[quadNames, r56Values, quadCurrents] = loadR56CurrentsTable(140);
quadNames = strrep(quadNames,'.','.I');

deltaP = 0.0002;
finalPhaseR56 = NaN(length(deltaP),length(r56Values));
finalR56 = NaN(length(deltaP),length(r56Values));
for p=1:length(deltaP)
    fprintf('-----------------------\n');
    fprintf('p = %.3f\n',deltaP(p));
    fprintf('-----------------------\n');
    
    for i=1:length(r56Values)
        tic       
        fprintf('R56 = %.3f... ',r56Values(i));
        
        machineModel.setModelParameters(quadNames,quadCurrents(i,:));
        machineModel.setModelParameters('initialX',0);
        machineModel.setModelParameters('initialY',0);
        machineModel.setModelParameters('initialDELTAP',deltaP(p));
        auxOptics = machineModel.computeOptics();

        finalPhaseR56(p,i) = auxOptics.DATA(end).T;
        finalR56(p,i) = auxOptics.DATA(end).RE56;
        
        fprintf('%.3f seconds\n',toc);
    end
end

%%
plot(r56Values, finalPhaseR56*(360/0.025));
title(['from ', nominalOptic.DATA(3).NAME , ' to ', nominalOptic.DATA(end-2).NAME, ' (S = ',num2str(nominalOptic.DATA(end-2).S),' m)' ])
xlabel('R56 [m]')
ylabel('T [12GHz degrees]')
grid

%%
figure;
legLabels = cell(1,length(deltaP));
for p=1:length(deltaP)
    plot(r56Values,finalPhaseR56(p,:)*(360/0.025));
    hold on;
    legLabels{p} = ['\Delta p / p_0 = ' num2str(deltaP(p))];
end
title(['from ', nominalOptic.DATA(3).NAME , ' to ', nominalOptic.DATA(end-2).NAME, ' (S = ',num2str(nominalOptic.DATA(end-2).S),' m)' ])
xlabel('R56 [m]')
ylabel('T [12GHz degrees]')
legend(legLabels);
%%
plot(r56Values,finalR56);
title(['from ', nominalOptic.DATA(3).NAME , ' to ', nominalOptic.DATA(end-2).NAME, ' (S = ',num2str(nominalOptic.DATA(end-2).S),' m)' ])
xlabel('Set R56 [m]')
ylabel('Final R56 [m]')
grid

%% T566 vs R56 comparison using sector map


% new twiss command
mycommand=[...
    'select, flag=twiss, clear;',char(10),...
    'select, flag=twiss, column=name, s, BETX, ALFX, x, px, dx, BETY, ALFY, y, py, dy,mux,muy, t, pt,re56, t566;',char(10),...
    'twiss, RMATRIX,SECTORMAP,',char(10),...
    'BETX=initialBETX,ALFX=initialALFX,',char(10),...
    'BETY=initialBETY,ALFY=initialALFY,',char(10),...
    'X=initialX,Y=initialY,',char(10),...
    'DELTAP=initialDELTAP;',char(10),...
    ];
% command to get sectortable
myoutput = 'write, table=sectortable;';

% matlab interface to sectortable
myTFSlines=machineModel.generalExecMadxScript(mycommand, myoutput);
myTFStable = madx2matlab.parseTFSTable(myTFSlines);


finalPhaseR56 = NaN(length(r56Values),1);
finalT566 = NaN(length(r56Values),1);
finalR56 = NaN(length(r56Values),1);
for i=1:length(r56Values)
    fprintf('R56 = %.3f\n',r56Values(i));
    machineModel.setModelParameters(quadNames,quadCurrents(i,:));
    machineModel.setModelParameters('initialX',0);
    machineModel.setModelParameters('initialY',0);
    machineModel.setModelParameters('initialDELTAP',0);
    tic
    myTFSlines=machineModel.generalExecMadxScript(mycommand, myoutput);
    fprintf('To execute MADX script: %.3f\n',toc);
    tic
    myTFStable = madx2matlab.parseTFSTable(myTFSlines);
    fprintf('To parse TFS table: %.3f\n',toc);
    
    tic
    finalRMatrix = eye(6,6);
    finalTMatrix = eye(6,6,6);
    for d=1:length(myTFStable.DATA)
        newRMatrix = NaN(6,6);
        newTMatrix = NaN(6,6,6);
        for r=1:6
            for c=1:6
                eval(sprintf('newRMatrix(r,c) = myTFStable.DATA(d).R%d%d;',r,c));
                for k=1:6
                    eval(sprintf('newTMatrix(r,c,k) = myTFStable.DATA(d).T%d%d%d;',r,c,k));
                end
            end
        end
        finalRMatrix = newRMatrix*finalRMatrix;
        finalTMatrix = newTMaxtrix*finalTMatrix;
    end
    fprintf('To construct matrices: %.3f\n',toc);
    
    finalR56(i) = finalRMatrix(5,6);
    finalT566(i) = finalTMatrix(5,6,6);
end

figure;
plot(r56Values,finalR56);
title(['from ', nominalOptic.DATA(3).NAME , ' to ', nominalOptic.DATA(end-2).NAME, ' (S = ',num2str(nominalOptic.DATA(end-2).S),' m)' ])
xlabel('Set R56 [m]')
ylabel('Final R56 [m]')
grid

figure;
plot(r56Values,finalT566);
title(['from ', nominalOptic.DATA(3).NAME , ' to ', nominalOptic.DATA(end-2).NAME, ' (S = ',num2str(nominalOptic.DATA(end-2).S),' m)' ])
xlabel('Set R56 [m]')
ylabel('Final T566 [m]')
grid

%%
finalRMatrix = eye(6,6);

for i=1:length(myTFStable.DATA)
    newRMatrix = NaN(6,6);
    for r=1:6
        for c=1:6
            eval(sprintf('newRMatrix(r,c) = myTFStable.DATA(i).R%d%d;',r,c));
        end
    end
    newRMatrix
    finalRMatrix = newRMatrix*finalRMatrix;
end

finalRMatrix(5,6)