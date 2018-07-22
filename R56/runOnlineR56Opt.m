clear all; close all;
%%

r56Obj = OnlineR56Opt();

r56Obj.saveData = true;
r56Obj.beamEnergy = 121;
r56Obj.nAvg = 80; % no. pulses to average over when calculating stats
r56Obj.initR56 = 0.0;
r56Obj.r56Step = 0.1;
r56Obj.phaseSampleRange = 589:689;
r56Obj.bpmSampleRange = 589:689;
r56Obj.wiggleGun = true;


%%

r56Obj.startOptimisation();