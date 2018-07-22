function [ fitConsts, fitErr ] = linFit2( x, y, w)
%linFit Linear fit with errors, removing any NaN values in data.
% errY is the standard error of each y value, used to weight fit. If not
% given all set to 1.

    sumW = sum(w);
    sumWX = sum(w.*x);
    sumWY = sum(w.*y);
    sumWXY = sum(w.*x.*y);
    sumWXSq = sum(w.*x.^2);
    
    delta = sumW.*sumWXSq - (sumWX.^2);
    
    intercept = ((sumWXSq.*sumWY)-(sumWX.*sumWXY))./delta
    slope = ((sumW.*sumWXY)-(sumWX.*sumWY))./delta
    
    sigE = sum(w.*((y-intercept-(slope.*x)).^2))./(length(x)-2);
    sigE = sqrt(sigE);
    
    intercept_err = sigE.*sqrt(sumWXSq./delta)
    slope_err = sigE.*sqrt(sumW./delta)
  
    figure;
    errorbar(x,y,w,'o');
    hold all;
    plot(x,slope.*x + intercept);

end

