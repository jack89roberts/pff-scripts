function [ newCurrents, oldCurrents, quadNames ] = setR56TL1( r56ValueToSet, beamEnergy )
%setR56TL1 Changes quadrupole currents in TL1 to vary R56.
%   R56 vaule to set must exist in the devices/currentsVsR56TL1.dat table.

    % load table of currents
    r56Table = importdata('/nfs/cs-ccr-nfs6/vol29/share/ctf/data/ctfmod/MatLab/PhaseFeedforward/devices/currentsVsR56TL1.dat');
    quadNames = strtrim(r56Table.colheaders(2:end));
    tableR56Values = r56Table.data(:,1);
    tableQuadCurrents = r56Table.data(:,2:end);
    tableQuadCurrents = tableQuadCurrents.*(beamEnergy./100); % scale currents to beam energy
    
    r56ValueIndex = find(ismember(tableR56Values,r56ValueToSet)); 
    if (isempty(r56ValueIndex))
        % if couldn't find the R56 value in the table, use the closest
        % value to the one the user gave.
        oldR56ValueToSet = r56ValueToSet;
        diffR56Values = tableR56Values-r56ValueToSet;
        r56ValueToSet = min(abs(diffR56Values));
        warning('R56 value of %.3f does not exist in devices/currentsVsR56TL1.dat. Will be set to %.3f instead.',oldR56ValueToSet,r56ValueToSet);
    end
    
    % get quad currents in the machine
    selector = 'SCT.USER.SETUP';
    oldCurrents = NaN(1,length(quadNames));
    fprintf('-------------------------------------\n')
    fprintf('Getting current Quad settings...\n')
    fprintf('-------------------------------------\n')

    for i=1:length(quadNames)
        oldCurrents(i) = matlabJapc.staticGetSignal(selector, sprintf('%s/SettingPPM#current',quadNames{i}));
        fprintf('%s: %.2f A\n',quadNames{i},oldCurrents(i)); 
    end

    % set new currents (twice for safety)
    fprintf('-------------------------------------\n')
    fprintf('Setting new Quad currents...\n')
    fprintf('-------------------------------------\n')

    newCurrents = tableQuadCurrents(r56ValueIndex,:);
    for i=1:length(quadNames)
        matlabJapc.staticSetSignal(selector, sprintf('%s/SettingPPM#current',quadNames{i}), newCurrents(i));
    end
    pause(0.5);
    for i=1:length(quadNames)        
        fprintf('%s: %.2f A\n',quadNames{i},newCurrents(i)); 
        matlabJapc.staticSetSignal(selector, sprintf('%s/SettingPPM#current',quadNames{i}), newCurrents(i));
    end
    
    fprintf('-------------------------------------\n')
    fprintf('R56 in TL1 set to %.3f m\n',r56ValueToSet);
end

