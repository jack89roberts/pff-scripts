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
    'twiss, RMATRIX,',char(10),...
    'BETX=initialBETX,ALFX=initialALFX,',char(10),...
    'BETY=initialBETY,ALFY=initialALFY,',char(10),...
    'X=initialX,Y=initialY,',char(10),...
    'PX=initialPX,',char(10),...
    'DELTAP=initialDeltaP;',char(10),...
    ];

kcktot = 0.001;
kckratio = 0:0.1:4;

% initial values for twiss
initialBETX =    10;
initialBETY =    1.5;
initialALFX =    0;
initialALFY =    0;
initialX = 0;
initialPX = 0;
initialY = 0;
initialDeltaP = 0.0;

% twiss table indices for kicker entrances/exits and chicane exit
ind480Ent = 88;
ind480Exit = 103;
ind780Ent = 155;
ind780Exit = 171;
ind800Exit = 175;

%%



for r=1:length(kckratio)
    fprintf('%d out of %d: %.5f\n',r,length(kckratio),kckratio(r))
    
    % set initial settings and compute optics
    machineModel.setModelParameters('initialBETX', initialBETX);
    machineModel.setModelParameters('initialBETY', initialBETY);
    machineModel.setModelParameters('initialALFX', initialALFX);
    machineModel.setModelParameters('initialALFY', initialALFY);
    machineModel.setModelParameters('initialX', initialX);
    machineModel.setModelParameters('initialY', initialY);
    machineModel.setModelParameters('initialDELTAP', initialDeltaP);
    machineModel.setModelParameters('initialPX', initialPX);


    machineModel.setModelParameters('cc.kcktot', kcktot);
    machineModel.setModelParameters('cc.kckratio', kckratio(r));
    
    nominalOptic = machineModel.computeOptics();
    tmpX = [nominalOptic.DATA.X];
    tmpPX = [nominalOptic.DATA.PX];

    if (r==1)
        opticsX = NaN(length(kcktot),length(tmpX));
        opticsPX = NaN(length(kcktot),length(tmpPX));
    end
    
    opticsX(r,:) = tmpX;
    opticsPX(r,:) = tmpPX;
    
end

x480Ent = opticsX(:,ind480Ent);
x480Exit = opticsX(:,ind480Exit);
x780Ent = opticsX(:,ind780Ent);
x780Exit = opticsX(:,ind780Exit);
x800Exit = opticsX(:,ind800Exit);

px480Ent = opticsPX(:,ind480Ent);
px480Exit = opticsPX(:,ind480Exit);
px780Ent = opticsPX(:,ind780Ent);
px780Exit = opticsPX(:,ind780Exit);
px800Exit = opticsPX(:,ind800Exit);

figure;
plot(kckratio,x800Exit);
title('x800')

figure;
plot(kckratio,px800Exit);
title('px800')
