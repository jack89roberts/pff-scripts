function [quadNames, tableR56Values, tableQuadCurrents] = loadR56CurrentsTable(beamEnergy)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    % load table of currents
    if (strcmp(getenv('USER'),'jack') == 1)
        r56Table = importdata('/home/jack/Documents/MATLAB/PhaseFeedforward/ctfmod/PhaseFeedforward/devices/currentsVsR56TL1.dat');
    else
        r56Table = importdata('/nfs/cs-ccr-nfs6/vol29/share/ctf/data/ctfmod/MatLab/PhaseFeedforward/devices/currentsVsR56TL1.dat');
    end
    quadNames = strtrim(r56Table.colheaders(2:end));
    tableR56Values = r56Table.data(:,1);
    tableQuadCurrents = r56Table.data(:,2:end);
    tableQuadCurrents = tableQuadCurrents.*(beamEnergy./100); % scale currents to beam energy

end

