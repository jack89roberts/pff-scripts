function [val465,val765] = getCorrValues()
    corr465Name = 'CC.DHF0465';
    corr765Name = 'CC.DHF0765';
    corrAcqName = 'Acquisition#currentAverage';
    cycleName = 'SCT.USER.SETUP';
    
    devName465 = [corr465Name '/' corrAcqName];
    devName765 = [corr765Name '/' corrAcqName];

    val465 = matlabJapc.staticGetSignal(cycleName,devName465);
    val765 = matlabJapc.staticGetSignal(cycleName,devName765);
end
