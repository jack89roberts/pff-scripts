%% Variables used for line formatting in plot commands
% this script must be run first to create the variables if you want to use them 
lineWidthBig = 2;
lineWidthSmall = 1.5;
markerSizeBig = 7;
markerType = 'o';
smallLegendLocation = 'Best';
%% Call this script after making the plot (including all labels, legends etc.) to make some changes
%% to font sizes, legend etc.

% forces figures to be saved with same aspect ratio
set(gcf,'PaperPositionMode','auto');

% change font sizes
set(gca,'FontSize',14)
set(findall(gcf,'type','text'),'FontSize',14)

% add grid lines
grid on;

% format legend
figLegend = findobj(gcf,'Type','axes','Tag','legend');
if (~isempty(figLegend))
    set(figLegend,'FontSize',10);
    % If the legend has more than 3 entries, move it outside the plot area.
    if (length(get(figLegend,'String')) > 3) 
        set(figLegend,'Location','EastOutside')
    % Otherwise, try to find the best place for the legend inside the plot area.
    else
        set(figLegend,'Location',smallLegendLocation)
    end
    
    % fix bug where legend moves in saved images
    set(figLegend,'Location','none');
end