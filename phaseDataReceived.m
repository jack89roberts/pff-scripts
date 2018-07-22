
function phaseDataReceived(lastDataStruct)
    % Check a device from each monitor, print a message if data inside is
    % empty or NaN.
    try 
        tmp = lastDataStruct.CC_SVBPI0535S.Samples.samples.value;
        ccBPMIsBad = (isnan(max(tmp)) || isempty(tmp));
    catch
        ccBPMIsBad = 1;
    end
    
    try
        tmp = lastDataStruct.CB_SVBPS0210S.Samples.samples.value;
        cbBPMIsBad = (isnan(max(tmp)) || isempty(tmp));
    catch
       cbBPMIsBad = 1; 
    end
    
    try
        tmp = lastDataStruct.CT_SVBPM0285S.Samples.samples.value;
        ctBPMIsBad = (isnan(max(tmp)) || isempty(tmp));
    catch
       ctBPMIsBad = 1;  
    end
    
        
    try
        tmp = lastDataStruct.CT_SCOPE01_CH01.Acquisition.value.value;
        frascatiScopeIsBad = (isnan(max(tmp)) || isempty(tmp));
    catch
       frascatiScopeIsBad = 1;   
    end
    
    try
        tmp = lastDataStruct.CE_SCOPE03_CH01.Acquisition.value.value;
        petsScopeIsBad = (isnan(max(tmp)) || isempty(tmp));
    catch
        petsScopeIsBad = 1;   
    end
    
    try
        tmp = lastDataStruct.CC_STBPR0915S.Samples.samples.value;
        bprIsBad = (isnan(max(tmp)) || isempty(tmp));
    catch
        bprIsBad = 1;
    end
    
    try    
        tmp = lastDataStruct.CK_SVPKI02P.Samples.samples.value;
        pkiIsBad =  (isnan(max(tmp)) || isempty(tmp));
    catch
        pkiIsBad = 1;
    end
      
    try
        tmp = lastDataStruct.CK_SVPSI03P.Samples.samples.value;
        psiIsBad = (isnan(max(tmp)) || isempty(tmp));
    catch
        psiIsBad = 1;
    end
       
    try
        tmp = lastDataStruct.CL_SVBPM0402S.Samples.samples.value;
        clBPMIsBad = (isnan(max(tmp)) || isempty(tmp));
    catch
        clBPMIsBad = 1;
    end
       
    try
        tmp = lastDataStruct.CR_SVBPI0275S.Samples.samples.value;
        crBPMIsBad = (isnan(max(tmp)) || isempty(tmp));   
    catch
        crBPMIsBad = 1;
    end
    
    dataIsBad = [frascatiScopeIsBad,petsScopeIsBad,ccBPMIsBad,cbBPMIsBad,ctBPMIsBad, clBPMIsBad,crBPMIsBad, bprIsBad,pkiIsBad,psiIsBad];
    
    if sum(dataIsBad) > 0
        dataNames = {'Frascati', 'PETS', 'CC BPM', 'CB BPM', 'CT BPM', 'CL BPM', 'CR BPM', 'BPR', 'PKI', 'PSI'};
        warningStr = 'Bad data:';
        for i=1:length(dataIsBad)
            if (dataIsBad(i))
                warningStr = [warningStr ', ' dataNames{i}];
            end
        end
        warning('%s',warningStr);
    end
    
    
    %phaseBPRs;
    %phaseFrascati;
    %phasePETS;

    %orbitBPMs;
end
