function phase = getPhaseSin( output, calibA, calibB, calibC )
% BPR phase calculation given input calibration constants
%   output = Asin(bx+c)
%   phase = [asin(output/A)-c]/b

    phase =  (asin((output./calibA))-calibC)./calibB ;
    

end

