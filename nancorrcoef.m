function [ R,err,P,RLO,RUP ] = nancorrcoef( x, y )
%NANCORRCOEF Corrcoef ignoring nan values
%   

    if (iscolumn(x) ~= iscolumn(y))
        y = y';
    end

    notNaNX = ~isnan(x);
    notNaNY = ~isnan(y);
    notNaNBoth = notNaNX & notNaNY;
    
    if (sum(notNaNBoth)>2)
        validX = x(notNaNBoth);
        validY = y(notNaNBoth); 
        [ R,P,RLO,RUP ] = corrcoef(validX, validY);    
        
        R = R(1,2);
        P = P(1,2);
        N = length(validX);
        err = sqrt((1-(R.^2))./(N-2));

    else
        R = NaN;
        err = NaN;
        P = NaN;
        RLO = NaN(2,2);
        RUP = NaN(2,2);
        %error('Less than 2 data pairs with non-NaN values');
    end
    
end

