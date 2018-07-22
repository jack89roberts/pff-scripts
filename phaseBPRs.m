    % need to do a bpr calibration for CL.STBPR0290S
    phaseCL475 = real(getPhaseSin(lastDataStruct.CL.STBPR0475S.Samples.samples.value, 11.719413, 0.017954, 0.582959));
    phaseCT532 = real(getPhaseSin(lastDataStruct.CT.STBPR0532S.Samples.samples.value, 123.533712, 0.017031, 0.255416));
    phaseCR505 = real(getPhaseSin(lastDataStruct.CR.STBPR0505S.Samples.samples.value, 0.770575, 0.018457, -1.813096));

    subplot(4,2,1); 
    plot(lastDataStruct.CL.STBPR0290S.Samples.samples.value);
    title('CL.STBPR0290');
    
    subplot(4,2,2);
    plot(phaseCL475);
    title('CL.STBPR0475');
    
    subplot(4,2,4); 
    plot(phaseCT532);
    title('CT.STBPR0532');
    
    subplot(4,2,5);
    plot(phaseCR505);
    title('CR.STBPR0505');