function [currents, quadNames ] = getCurrentsTL1()
%getCurrentsTL1 Get quadrupole currents in TL1 (same format as setR56 script).

    % load R56 table to get quad names in order
    r56Table = importdata('/nfs/cs-ccr-nfs6/vol29/share/ctf/data/ctfmod/MatLab/PhaseFeedforward/devices/currentsVsR56TL1.dat');
    quadNames = strtrim(r56Table.colheaders(2:end));

    % get quad currents in the machine
    selector = 'SCT.USER.SETUP';
    currents = NaN(1,length(quadNames));
    fprintf('-------------------------------------\n')
    fprintf('Getting current Quad settings...\n')
    fprintf('-------------------------------------\n')

    for i=1:length(quadNames)
        currents(i) = matlabJapc.staticGetSignal(selector, sprintf('%s/SettingPPM#current',quadNames{i}));
        fprintf('%s: %.2f A\n',quadNames{i},currents(i)); 
    end

end

