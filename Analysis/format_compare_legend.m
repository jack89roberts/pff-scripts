if (nDataSets > 6)
    set(gcf, 'Colormap', datSetColors);
    figColBar = colorbar;
    
    allTicks = linspace(0,1,nDataSets+1);
    if (nDataSets>8)
        set(figColBar,'YTick',allTicks(1:2:nDataSets));
%         set(figColBar,'YTick',1:2:nDataSets);
        set(figColBar,'YTickLabel',dataSetValues(1:2:end));
    else
        set(figColBar,'YTick',allTicks(1:nDataSets));
%         set(figColBar,'YTick',1:nDataSets);
        set(figColBar,'YTickLabel',dataSetValues);   
    end
    ylabel(figColBar,dataSetValueLabel);
else
    legend(legLabels,'Location','NorthWest');
end