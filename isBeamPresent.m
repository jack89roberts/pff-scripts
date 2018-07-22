function [ beamIsPresent ] = isBeamPresent()
%ISBEAMPRESENT Returns 1 if gun authorised and enabled, 0 otherwise
%   correctStatus is the value returned by GUN-AUTH/Acquisition/status when
%   the gun is on.
% BUG - Correct status when beam enabled isn't one fixed value. Need to
% figure out what GUN-AUTH status actually means.
    nAttempts = 0;
    exitCount = 10;
    while (nAttempts < exitCount)
        try
            enabled=JGetCoValueFesa('CX.SGUN-CTTIM','EnableStatus','enableStatus');

            if (enabled==1) && (min(JGetCoArrayFesa('CL.STBPM0402S', 'Samples', 'samples')) < -2)
                beamIsPresent=1;
            else
                beamIsPresent=0;
            end
            
            return;
        
        catch
            nAttempts = nAttempts + 1;
        end    
    end
    
    disp('All attempts to identify beam failed');
    beamIsPresent = 0;
    
    
    
    
    %beamIsPresent = 1;
end

