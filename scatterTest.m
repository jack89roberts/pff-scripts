function [ ] = scatterTest( x,y,z)
%scatterTest 
    % colour and size of markers depends on given z values
    % sine (blue -> green -> red)

    zScale = z-min(z);
    zScale = zScale./max(zScale);
    
    colScale = (1.2*pi)*zScale;
    colR = sin(colScale+(3*pi/2));
    colG = sin(colScale);
    colB = sin(colScale+(pi/2));
    colours = zeros(length(y),3);
    colours(:,1) = colR;
    colours(:,2) = colG;
    colours(:,3) = colB;
    
    sizes = 64*((zScale/2)+0.5);
    sizes = (64*zScale)+8;
    figure;
    subplot(4,1,4);
    scatter(z,ones(1,length(y)),sizes,colours,'o','filled');
    xlabel('z');
    set(gca,'YTick',[]);
    subplot(4,1,1:3);
    scatter(x,y,sizes,colours,'o','filled');
    xlabel('x');
    ylabel('y');
    
end

