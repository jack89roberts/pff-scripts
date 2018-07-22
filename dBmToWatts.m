function [ watts ] = dBmToWatts( dBm )
%dBmToWatts Converts power in dBm to Watts

    watts = (10.^(dBm./10))./1000;

end

