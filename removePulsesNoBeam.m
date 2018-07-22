function CTFData = removePulsesNoBeam( CTFData, refName, refLimit )
%removePulsesNoBeam Sets pulses with no beam to NaN for all signals in
%CTFData
%   CTFData: Array of matlabJapc monitor data structs
%   refName: Name of signal to be used for checking whether beam is present
%   (default: CL_SVBPM0402S.Samples.samples.value)
%   refLimit: Beam defined as present if max(abs(refName))>=refLimit.
%   Default 2 if refName not given. If refName given must be user defined.
    if (nargin==2)
        error('If refName is given, refLimit must be defined.');
    end
    if (nargin<2)
        refName = 'CL_SVBPM0402S.Samples.samples.value';
        refLimit = 2;
    end
    
    try
        refSignal = extractCTFSignalFromMergedData(refName,CTFData);
    catch
        fprintf('!!!Ref signal not found in CTFData. No pulses removed.\n');
        return;
    end
    
    maxSignal = max(abs(refSignal),[],2);
    pulseHasBeam = maxSignal>=refLimit;
       
    signalList = CTFData(1).parameters;   
    signalList = strrep(signalList,'.','_');
    signalList = strrep(signalList,'-','_');
    signalList = strrep(signalList,'/','.');
    signalList = strrep(signalList,'#','.');
    signalList = strcat(signalList,'.value'); 
    
    nPulses = length(CTFData);
    nSignals = length(signalList);

    for p=1:nPulses
        if (~pulseHasBeam(p))
            for s=1:nSignals
                sigName = signalList{s};
                evalStr = sprintf('CTFData(p).%s = NaN(size(CTFData(p).%s));',sigName,sigName);
                eval(evalStr);
            end
        end
    end
    
end

