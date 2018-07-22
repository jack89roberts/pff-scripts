function [ fitConsts, fitErr ] = linFit( x, y, errY, plotResult )
%linFit Linear fit with errors, removing any NaN values in data.
% errY is the standard error of each y value, used to weight fit. If not
% given all set to 1.

    if (iscolumn(x) ~= iscolumn(y))
        y = y';
    end
    
    if (nargin==2 || isempty(errY))
%         errY = ones(size(y));
        useErr = 0;
    else
        useErr = 1;    
        if (iscolumn(errY) ~= iscolumn(y))
            errY = errY';
        end
    end
    
    
    if (nargin<4 || isempty(plotResult))
        plotResult = 0;
    end
        
    notNaNX = ~isnan(x);
    notNaNY = ~isnan(y);
    if (useErr)
        notNaNerrY = ~isnan(errY);
        notNaNAll = notNaNX & notNaNY & notNaNerrY;
        errY = errY(notNaNAll);
    else
       notNaNAll = notNaNX & notNaNY; 
    end
    x = x(notNaNAll);
    y = y(notNaNAll);
    N = length(x);
    
    if (N>2 && useErr)
        sumXoverESq = sum(x./(errY.^2));
        sumYoverESq = sum(y./(errY.^2));
        sumXYoverESq = sum((x.*y)./(errY.^2));
        sum1OverEsq = sum(1./(errY.^2));
        sumXSqOverESq = sum((x.^2)./(errY.^2));
        
        slope = (sumXoverESq.*sumYoverESq)-(sumXYoverESq.*sum1OverEsq);
        slope = slope./((sumXoverESq.^2)-(sumXSqOverESq.*sum1OverEsq));
        
        intercept = (sumXYoverESq - slope.*sumXSqOverESq)./sumXoverESq;
        
        slope_err = sum1OverEsq./((sumXSqOverESq.*sum1OverEsq)-(sumXoverESq.^2));
        slope_err = sqrt(slope_err);
        
        intercept_err = sumXSqOverESq./((sumXSqOverESq.*sum1OverEsq)-(sumXoverESq.^2));
        intercept_err = sqrt(intercept_err);
        
%         sigE = sum((1./(errY.^2)).*((y-intercept-slope.*x).^2))./(N-2);
%         slope_err = slope_err.*sigE;
%         intercept_err = intercept_err.*sigE;
%         
    elseif (N>2 && ~useErr)
        sumXY = sum(x.*y);
        sumX = sum(x);
        sumY = sum(y);
        sumXSq = sum(x.^2);
        
        slope = ((N*sumXY)-(sumX.*sumY))./((N.*sumXSq)-(sumX.^2));
        intercept = ((sumXSq.*sumY)-(sumX.*sumXY))./((N.*sumXSq)-(sumX.^2));
        
        S = sqrt(sum((y-(slope.*x)-intercept).^2)./(N-2));
        slope_err = S.*sqrt(N./((N.*sumXSq)-(sumX.^2)));
        intercept_err = S.*sqrt(sumXSq./((N.*sumXSq)-(sumX.^2)));
        
    else
        slope = NaN;
        intercept = NaN;
        slope_err = NaN;
        intercept_err = NaN;
    end
    
    if (plotResult)
        fitResult = slope.*x+intercept;
        
        figure;
        subplot(1,2,1)
        if (useErr)
            errorbar(x,y,errY,'bo','MarkerFaceColor','b')
        else
            plot(x,y,'bo','MarkerFaceColor','b')
        end
        hold all;
        plot(x,fitResult,'-r');
        xlabel('x')
        ylabel('y')
        legend('DATA','FIT')
        title(sprintf('y = (%.2f%c%.2f)*x + (%.2f%c%.2f)',slope,char(177),slope_err,intercept,char(177),intercept_err));
        
        subplot(1,2,2)
        plot(x,y-fitResult,'bo','MarkerFaceColor','b');
        title(sprintf('RESIDUALS (std=%.4f)',std(y-fitResult)))
        xlabel('x')
        ylabel('y-fitY')
        
    end
    
    fitConsts = [slope, intercept];
    fitErr = [slope_err, intercept_err];


end

