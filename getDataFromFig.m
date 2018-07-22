function [ XData,YData,Names, D] = getDataFromFig(ax)
%getDataFromFig Extracts data points from axis ax.
%   [ XData,YData,Names ] = getDataFromFig(ax)
%   XData: x values for each data series.
%   YData: y values for each data series.
%   Names: name of each data series.

    if (nargin<1)
        ax = gca;
    end

    D = get(ax,'Children');
    XData = get(D,'XData');
    YData = get(D,'YData');
    Names = get(D,'DisplayName');

end

