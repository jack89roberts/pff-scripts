function [coeffs,rsquare,confInts] = nanpolyfit( x, y, n, weights )
%NANPOLYFIT Polyfit ignoring nan values
% [coeffs,rsquare,confInts] = nanpolyfit( x, y, n, weights )

    if (iscolumn(x) ~= iscolumn(y))
        y = y';
    end
    
    if (nargin==3 || isempty(weights))
        weights = ones(size(y));
    end
    
    if (iscolumn(weights) ~= iscolumn(y))
        weights = weights';
    end

    notNaNX = ~isnan(x);
    notNaNY = ~isnan(y);
    notNaNWeights = ~isnan(weights);
    notNaNAll = notNaNX & notNaNY & notNaNWeights;
    validX = x(notNaNAll);
    validY = y(notNaNAll);
    validWeights = weights(notNaNAll);

    if (sum(notNaNAll)>2)
        [xData, yData, weightsData] = prepareCurveData( validX, validY, validWeights );
        ft = fittype( sprintf('poly%d',n) );

        opts = fitoptions( 'Method', 'LinearLeastSquares' );
        opts.Weights = weightsData;
        
        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts );

        coeffs = coeffvalues(fitresult);
        rsquare = gof.rsquare;
        confInts = confint(fitresult);
        
    else
        coeffs = NaN(1,n+1);
        rsquare = NaN;
        confInts = NaN(2,n+1);     
    end

% OLD FUNCTION
%function [p,S] = nanpolyfit( x, y, n )
%     if (sum(notNaNBoth)>2)
%         [p,S] = polyfit(validX, validY, n);
%     else
%         %error('Less than 2 data pairs with non-NaN values');
%         p = NaN(1,n);
%         S = {};
%         S.error = 'Less than 2 data pairs with non-NaN values';
%     end
    
end

