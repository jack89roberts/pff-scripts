function [ ] = scatterColourProgression( x,y,markerType,markerSize,showColours)
%scatterColourProgression Makes a scatter plot where colour of each marker
%varies from blue to green to red depending on it's position in the array
%(first points blue, middle green, end red).

%   Detailed explanation goes here
    
    % linear blue to red
%     colours = zeros(length(y),3);
%     colours(:,3) = 1;
%     colourStep = 1./length(y);
%     for i=2:length(y)
%         colours(i,:) = [colours(i-1,1)+colourStep, 0, colours(i-1,3)-colourStep];
% %         colours(i,:) = [colours(i-1,1)+colourStep, colours(i-1,2)+colourStep, colours(i-1,3)+colourStep];
%     end
        
    % linear blue to green to red
%     colours = zeros(length(y),3);
%     colours(1,:) = [0 0 1];
%     mid = ceil(length(y)/2);
%     stepStart = 1./(mid-1);
%     for i=2:mid
%         colours(i,:) = [0 colours(i-1,2)+stepStart colours(i-1,3)-stepStart];
%     end
%     stepEnd = 1./(length(y)-mid);
%     for i=mid+1:length(y)
%         colours(i,:) = [colours(i-1,1)+stepEnd colours(i-1,2)-stepEnd 0];
%     end
    
    % sine (blue -> green -> red)
    colX = linspace(0,1.2*pi,length(y));
    colR = sin(colX+(3*pi/2));
    colG = sin(colX);
    colB = sin(colX+(pi/2));
    colours = zeros(length(y),3);
    colours(:,1) = colR;
    colours(:,2) = colG;
    colours(:,3) = colB;

    if (nargin<4)
        markerSize = 40;
    end
    if (nargin<3)
        markerType = 'o';
    end

    if (nargin<5 || showColours)
        figure;
        subplot(5,1,5);
%         plot(colR,'r');
%         hold all;
%         plot(colG,'g');
%         plot(colB,'b');
        scatter(1:length(y),ones(1,length(y)),64,colours,markerType,'filled');
%         title('Colour progression');
        xlabel('Colour of Each Data Point No.');
        set(gca,'YTick',[]);
        subplot(5,1,1:4);
    end
    
    scatter(x,y,markerSize,colours,markerType,'filled');
    
end

