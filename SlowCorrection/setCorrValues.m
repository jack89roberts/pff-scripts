function setCorrValues(setVal465,setVal765)
    corr465Name = 'CC.DHF0465';
    corr765Name = 'CC.DHF0765';
    corrSetName = 'SettingPPM#current';
    cycleName = 'SCT.USER.SETUP';

    fprintf('Changing correctors...\n');
    setName465 = [corr465Name '/' corrSetName];
    setName765 = [corr765Name '/' corrSetName];

    changed465 = 0;
    changed765 = 0;           
    while (~changed465 || ~changed765)
        matlabJapc.staticSetSignal(cycleName,setName465,setVal465);
        matlabJapc.staticSetSignal(cycleName,setName765,setVal765);
        pause(4);

        % check that new acquired corrector currents agree with set
        % value to +/- 0.02. If not try again.
        [new465,new765] = getCorrValues();
        changed465 = abs(new465-setVal465) <= 0.02;
        changed765 = abs(new765-setVal765) <= 0.02;

        fprintf('Set 465: %.2f, Acq 465: %.2f\n',setVal465,new465);
        fprintf('Set 765: %.2f, Acq 765: %.2f\n',setVal765,new765);
    end

end