function [ phase ] = getPhaseNewNormalisation(mixer,diode,a,b,c,d,e)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

    numerator = mixer-c;
    %denominator = 10.^(((a./d)*(sqrt(abs(meanDiode))-e))+b);
    denominator = sqrt(abs(diode))-e;
    denominator = denominator.*(a/d);
    denominator = denominator + b;
    denominator = 10.^denominator;
    
    phase = numerator./denominator;
    phase(phase>1) = 1;
    phase(phase<-1) = -1;
        
    phase = asind(phase);
    
end

