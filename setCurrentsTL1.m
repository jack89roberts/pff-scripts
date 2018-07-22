function [] = setCurrentsTL1( currents )
%setCurrentsTL1 Changes quadrupole currents in TL1 to given values (assuming currents comes from output of setR56TL1).
%   R56 vaule to set must exist in the devices/currentsVsR56TL1.dat table.

    % load table of currents to get quad names
    r56Table = importdata('devices/currentsVsR56TL1.dat');
    quadNames = strtrim(r56Table.colheaders(2:end));
    selector = 'SCT.USER.SETUP';

    % set new currents (twice for safety)
    fprintf('-------------------------------------\n')
    fprintf('Setting new Quad currents...\n')
    fprintf('-------------------------------------\n')

    for i=1:length(quadNames)
        matlabJapc.staticSetSignal(selector, sprintf('%s/SettingPPM#current',quadNames{i}), currents(i));
    end
    pause(0.5);
    for i=1:length(quadNames)        
        fprintf('%s: %.2f A\n',quadNames{i},currents(i)); 
        matlabJapc.staticSetSignal(selector, sprintf('%s/SettingPPM#current',quadNames{i}), currents(i));
    end
end

