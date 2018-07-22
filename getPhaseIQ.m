function [phase,power] = getPhaseIQ( inpI, inpQ )
% Returns array of phase (degrees) from arrays of I and Q values
%   for PETS phase measurement
    
    % subtract offset
    I = inpI - nanmean(inpI(1:(round(length(inpI)/10))));
    Q = inpQ - nanmean(inpQ(1:(round(length(inpQ)/10))));

    % calculate phase
    phase = atan2(double(Q),double(I));
    phase = phase.*(180./pi);
%     for i=1:length(phase)
%         if(phase(i)<0) 
%             phase(i) = phase(i) + 180;
%         end
%     end
    power = sqrt(I.^2 + Q.^2);
    

end

