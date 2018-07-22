function [ chiSqVal ] = chisq( observed, expected)
% calculates chi squared value from observed and expected values
    
    chiSqVal = sum( ((observed-expected).^2)./abs(expected));
    
end

