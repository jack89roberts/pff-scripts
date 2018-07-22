function [  corrH, corrV ] = extractCorrFromCTFData( corrName, CTFData)
%extractCorrFromCTFData Returns corrector currents from CTFData
%   corr name must not include H or V - added here.
    nPulses = length(CTFData);    
  
    corrAcqNameH = [corrName(1:4) 'H' corrName(5:end)];
    corrAcqNameV = [corrName(1:4) 'V' corrName(5:end)];
    
    corrH = NaN(1,nPulses);
    corrV = NaN(1,nPulses);
    
    nBadH = 0;
    nBadV = 0;
    
    for t=1:nPulses
        try
            tmpH = eval(['CTFData(t).' corrAcqNameH '.Acquisition.current.value;']);
            tmpH = tmpH(1);
            corrH(t) = tmpH;
        catch
            corrH(t) = NaN;
            nBadH = nBadH+1;
        end
        try
            tmpV = eval(['CTFData(t).' corrAcqNameV '.Acquisition.current.value;']);
            tmpV = tmpV(1);
            corrV(t) = tmpV;
        catch
            corrV(t) = NaN;
            nBadV = nBadV+1;
        end
    end
    
    if (nBadH>0 || nBadV>0)
        fprintf('%s bad data: H - %d pulses, V - %d pulses\n',corrName,nBadH,nBadV);
    end

end

