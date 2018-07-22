set(gcf,'Colormap',binColours)
figColBar = colorbar;
%set(figColBar,'YTickLabel',legendLabels)
set(figColBar,'YTick',(1:nWindows)/nWindows);
set(figColBar,'YTickLabel',legendLabels)
ylabel(figColBar, 'Cut on upstream phase [no. sigma]')