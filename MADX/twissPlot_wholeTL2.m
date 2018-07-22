addpath('/home/jack/Documents/MATLAB/madx2Matlab');

twissTable=madx2matlab.parseTFSTableFromFile('/home/jack/Documents/PFF Chicane Design/TL2CleanScripts/output/WholeTL2Line.twiss');
% madx2matlab.plotMachineTwissFunctions(twissTable)	
% madx2matlab.plotMachineTwissElements(twissTable)

%% plot twiss functions

% extract functions
auxS = [twissTable.DATA.S];

auxBetaX = [twissTable.DATA.BETX];
auxBetaY = [twissTable.DATA.BETY];
auxDispX = [twissTable.DATA.DX];
auxDispY = [twissTable.DATA.DY];
auxR56 = [twissTable.DATA.RE56];
auxX = [twissTable.DATA.X];
auxT = [twissTable.DATA.T];

figure
output_axes = [];
%subplot('Position',positionVector) positionVector := [left,bottom,width,height]

%% plot the machine elements

output_axes(1) = subplot('Position',[0.1,0.1,0.85,0.2]);

yscale = 1;

% some baseline
minS = min([twissTable.DATA.S]);
maxS = max([twissTable.DATA.S]);
line([minS maxS], [0 0],'Color','k');

% SBENDs
myIndex=madx2matlab.findElement(twissTable,'KEYWORD','SBEND');
madx2matlab.plotTwissSingleElements(myIndex, twissTable, 'bend', yscale);

% RBENDs
myIndex=madx2matlab.findElement(twissTable,'KEYWORD','RBEND');
madx2matlab.plotTwissSingleElements(myIndex, twissTable, 'bend', yscale);

% % HKICKERs
% myIndex=madx2matlab.findElement(twissTable,'KEYWORD','HKICKER');
% madx2matlab.plotTwissSingleElements(myIndex, twissTable, 'kicker', yscale);
% 
% % VKICKERs
% myIndex=madx2matlab.findElement(twissTable,'KEYWORD','VKICKER');
% madx2matlab.plotTwissSingleElements(myIndex, twissTable, 'kicker', yscale);
% 
% % KICKERs
% myIndex=madx2matlab.findElement(twissTable,'KEYWORD','KICKER');
% madx2matlab.plotTwissSingleElements(myIndex, twissTable, 'kicker', yscale);

% QUADRUPOLEs
myIndex=madx2matlab.findElement(twissTable,'KEYWORD','QUADRUPOLE');
% remove 490 650 790 which are split in to multiple elements (done manually
% below)
myQuadNames = {twissTable.DATA(myIndex).NAME};
indsNot490 = ~cellfun('length',strfind(myQuadNames,'CC.QDH0490'));
indsNot650 = ~cellfun('length',strfind(myQuadNames,'CC.QDL0650'));
indsNot790 = ~cellfun('length',strfind(myQuadNames,'CC.QDH0790'));
myIndex = myIndex(indsNot490 & indsNot650 & indsNot790);
myQuadsK1 = [twissTable.DATA(myIndex).K1L];
myFocusingIdxes = myIndex(myQuadsK1 > 0);
myDeFocusingIdxes = myIndex(myQuadsK1 < 0);
madx2matlab.plotTwissSingleElements(myFocusingIdxes, twissTable, 'focusing', yscale);
madx2matlab.plotTwissSingleElements(myDeFocusingIdxes, twissTable, 'defocusing', yscale);

% 
% % SEXTUPOLEs
% myIndex=madx2matlab.findElement(twissTable,'KEYWORD','SEXTUPOLE');
% madx2matlab.plotTwissSingleElements(myIndex, twissTable, 'sextupole', yscale);
% 
% % MONITORs
% myIndex=madx2matlab.findElement(twissTable,'KEYWORD','MONITOR');
% madx2matlab.plotTwissSingleElements(myIndex, twissTable, 'monitor', yscale);

%% manually add elements (kickers/quads) that are split in to multiple types in MADX
%
kickerXCoordinates = [0 1 1 0];
kickerYCoordinates = [0.25 0.25 -0.25 -0.25];
%
focusingXCoordinates = [0 0.2 0.5 0.8 1  0.8 0.5  0.2];
focusingYCoordinates = [0 0.7   1 0.7 0 -0.7  -1 -0.7];
%
defocusingXCoordinates = [0 0.3 0.45  0.3  0  1  0.7 0.55 0.7 1];
defocusingYCoordinates = [1 0.7   0 -0.7 -1 -1 -0.7   0 0.7 1];

%QDL650
tmpLength = 0.2950;
tmpInitialPosition  = 24.10521369;
tmpLineStyle = '-';
tmpFaceColor = 'b';
tmpXdata = focusingXCoordinates*tmpLength + tmpInitialPosition;
tmpYdata = focusingYCoordinates*yscale;
patch(tmpXdata, tmpYdata,'k','FaceColor',tmpFaceColor,'LineStyle',tmpLineStyle)

%QDH490
tmpLength = 0.38;
tmpInitialPosition  = 17.37496271;
tmpLineStyle = '-';
tmpFaceColor = 'b';
tmpXdata = focusingXCoordinates*tmpLength + tmpInitialPosition;
tmpYdata = focusingYCoordinates*yscale;
patch(tmpXdata, tmpYdata,'k','FaceColor',tmpFaceColor,'LineStyle',tmpLineStyle)

%QDH790
tmpLength = 0.38;
tmpInitialPosition  = 28.23648083;
tmpLineStyle = '-';
tmpFaceColor = 'b';
tmpXdata = defocusingXCoordinates*tmpLength + tmpInitialPosition;
tmpYdata = defocusingYCoordinates*yscale;
patch(tmpXdata, tmpYdata,'k','FaceColor',tmpFaceColor,'LineStyle',tmpLineStyle)

%KHD0480
tmpLength =  0.7880;
tmpInitialPosition  = 17.1141271;
tmpLineStyle = '-';
tmpFaceColor = 'g';
tmpXdata = kickerXCoordinates*tmpLength + tmpInitialPosition;
tmpYdata = kickerYCoordinates*yscale;
patch(tmpXdata, tmpYdata,'k','FaceColor',tmpFaceColor,'LineStyle',tmpLineStyle)

%KHD0780
tmpLength =  0.7880;
tmpInitialPosition  = 27.97563083;
tmpLineStyle = '-';
tmpFaceColor = 'g';
tmpXdata = kickerXCoordinates*tmpLength + tmpInitialPosition;
tmpYdata = kickerYCoordinates*yscale;
patch(tmpXdata, tmpYdata,'k','FaceColor',tmpFaceColor,'LineStyle',tmpLineStyle)


%%
ylim([-1 1])
set(output_axes(1), 'YTick', [], 'Box', 'on', 'FontSize',16)
xlabel('s [m]','FontSize',16)
grid on

%% plot the twiss functions

% prepare axes for twiss functions
output_axes(2) = subplot('Position',[0.1,0.3,0.85,0.65]);

% plot(auxS, auxBetaX,'Color','b','LineWidth',2)
% hold on
% plot(auxS, auxBetaY,'r','LineWidth',2)
% legend('\beta_x','\beta_y')
% ylabel('Beta [m]','FontSize',16)

% plot(auxS, auxDispX,'Color','b','LineWidth',2)
% hold on
% plot(auxS, auxDispY,'Color','r','LineWidth',2)
% legend('D_x','D_y')
% ylabel('Dispersion [m]','FontSize',16)

% plot(auxS, auxR56,'Color','b','LineWidth',2)
% legend('R56')
% ylabel('R56 [m]','FontSize',16)

% plot(auxS, auxX,'Color','b','LineWidth',2)
% legend('x')
% ylabel('x [m]','FontSize',16)

plot(auxS, (auxT/0.025)*360,'Color','b','LineWidth',2)
legend('Phase')
ylabel('Phase [degrees]','FontSize',16)


% hold off
grid on
% xlabel('s [m]','FontSize',16)
xlim([-0.1, max(auxS)+0.1])
set(gca,'FontSize',16)

linkaxes(output_axes,'x')
% set some paramaters
set(output_axes(2), 'XTickLabel', [], 'Box', 'on' ,'FontSize',16)

% fix xlim
% xlim(output_axes(2), [-0.1, max(auxS)+0.1])
xlim([0 36.8])

%% lines to highlight sections

% vertical chicane
line([10.1994215 10.1994215],get(gca,'YLim'),'Color','r','LineStyle','--','LineWidth',2)
line([15.45258002 15.45258002],get(gca,'YLim'),'Color','r','LineStyle','--','LineWidth',2)
t=text(12.8260, 0.95*max(get(gca,'YLim')), 'VERTICAL CHICANE','Color','r','FontSize',14);
t.HorizontalAlignment = 'Center';

% horizontal chicane
line([18.56406302 18.56406302],get(gca,'YLim'),'Color','r','LineStyle','--','LineWidth',2)
line([29.94136479 29.94136479],get(gca,'YLim'),'Color','r','LineStyle','--','LineWidth',2)
t=text(24.2527, 0.95*max(get(gca,'YLim')), 'HORIZONTAL CHICANE','Color','r','FontSize',14);
t.HorizontalAlignment = 'Center';

% 480 kicker
line([17.1141271 17.1141271],get(gca,'YLim'),'Color','g','LineStyle','--','LineWidth',2)
line([17.90211271 17.90211271],get(gca,'YLim'),'Color','g','LineStyle','--','LineWidth',2)
t=text(17.4631, 0.95*max(get(gca,'YLim')), 'K1','Color',[0 0.75 0],'FontSize',14);
t.HorizontalAlignment = 'Center';

% 780 kicker
line([27.97563083 27.97563083],get(gca,'YLim'),'Color','g','LineStyle','--','LineWidth',2)
line([28.76363083 28.76363083],get(gca,'YLim'),'Color','g','LineStyle','--','LineWidth',2)
t=text(28.3696, 0.95*max(get(gca,'YLim')), 'K2','Color',[0 0.75 0],'FontSize',14);
t.HorizontalAlignment = 'Center';



