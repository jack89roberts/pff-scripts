function [ volts ] = dBmToVolts( dBm, impedance )
%dBmToVolts Converts power in dBm to volts

    if nargin<2
        impedance = 50;
    end
    
    volts = sqrt(impedance./1000)*10.^(dBm./20);

end

