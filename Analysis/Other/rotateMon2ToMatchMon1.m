load('/home/jack/PhaseFeedforward/Analysis/201511/20151116_1400_Resolution_0dBPhMonInputs.mat');
%%
figure;
plot(timeAxis,meanPhaseAlongPulse(1,:),'b','LineWidth',2)
hold all
plot(timeAxis,meanPhaseAlongPulse(2,:),'r','LineWidth',2)
title('Upstream Phase Along Pulse')
xlabel('Time [ns]')
ylabel('Phase [degrees]')
xlim([80 1230])
legend('Mon1','Mon2')
format_plots

%% rotation about x axis (origin manually set where traces cross ~middle pulse)
x=(1:nSamples)-647;
y= meanPhaseAlongPulse(2,:);
calcRange = 530:695;
rotAng = 0:0.001:4;
sumDiffSq = NaN(1,length(rotAng));

for i=1:length(rotAng)
    R = [cosd(rotAng(i)) -sind(rotAng(i)); sind(rotAng(i)) cosd(rotAng(i))];
    myY = R*[x;y];
    myY = myY(2,:);

    sumDiffSq(i) = sum((meanPhaseAlongPulse(1,calcRange)-myY(calcRange)).^2);
end

figure;
plot(rotAng,sumDiffSq);

[a,b]=min(sumDiffSq);
fprintf('best rotation = %.3f\n',rotAng(b))

figure;   
plot(timeAxis,meanPhaseAlongPulse(1,:),'b','LineWidth',2);
hold all;
R = [cosd(rotAng(b)) -sind(rotAng(b)); sind(rotAng(b)) cosd(rotAng(b))];
myY = R*[x;y];
myY = myY(2,:);
plot(timeAxis,myY,'r','LineWidth',2)
title([sprintf('Upstream Phase Along Pulse: Mon2 Rotated %.1f',rotAng(b)) '\circ'])
xlabel('Time [ns]')
ylabel('Phase [degrees]')
xlim([80 1230])
legend('Mon1','Mon2')
format_plots

figure;
plot(timeAxis,meanPhaseAlongPulse(1,:)-meanPhaseAlongPulse(2,:),'b','LineWidth',2)
hold all;
plot(timeAxis,meanPhaseAlongPulse(1,:)-myY,'r','LineWidth',2);
title('Difference Between Mon1 and Mon2')
xlabel('Time [ns]')
ylabel('Phase [degrees]')
xlim([80 1200])
legend('Original','Mon2 Rotated')
format_plots

ang=2.127;
R = [cosd(ang) -sind(ang); sind(ang) cosd(ang)];
x=(1:nSamples)-647;
y=meanPhaseAlongPulse(2,:);
myY = R*[x;y];
myY = myY(2,:);
