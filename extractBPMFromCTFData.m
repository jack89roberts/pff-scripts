function [  bpmH, bpmS, bpmV ] = extractBPMFromCTFData( bpmName, CTFData)
%extractBPMFromCTFData Returns array of BPM signals (samples) for BPM with
%name bpmName in merged data CTFData.
%   BPM name must have same format as in the data struct, e.g. CT_SVBPM0285
%   with underscore instead of dot.
    
    nPulses = length(CTFData);    
    try
        nSamples = length(eval(['CTFData(1).' bpmName 'S.Samples.samples.value']));
    catch
         nSamples = length(eval(['CTFData(1).' strrep(bpmName,'_','.') 'S.Samples.samples.value']));       
    end
    
    if (nSamples < 100 || isempty(nSamples) || isnan(nSamples))
        fprintf('%s bad number of samples (%d).\n',bpmName,nSamples);
    end
    
    %fprintf('%s %d samples \n',bpmName,nSamples);
    
    bpmH = NaN(nPulses,nSamples);
    bpmS = NaN(nPulses,nSamples);
    bpmV = NaN(nPulses,nSamples);

    
    try
        eval(['CTFData(1).' bpmName 'V.Samples.samples.value;']);
        useVertical = 1;
    catch
        try
            eval(['CTFData(1).' strrep(bpmName,'_','.') 'V.Samples.samples.value;']);
            useVertical = 1;           
        catch
            useVertical = 0;
            bpmV = {};
        end
    end
    
    nBadH = 0;
    nBadS = 0;
    nBadV = 0;
    
    for t=1:nPulses
        try
            bpmH(t,:) = eval(['CTFData(t).' bpmName 'H.Samples.samples.value;']);
        catch
            try
                bpmH(t,:) = eval(['CTFData(t).' strrep(bpmName,'_','.') 'H.Samples.samples.value;']);
            catch
                bpmH(t,:) = NaN(1,nSamples);
                nBadH = nBadH+1;
            end
        end
        try
            bpmS(t,:) = eval(['CTFData(t).' bpmName 'S.Samples.samples.value;']);
        catch
            try
                bpmS(t,:) = eval(['CTFData(t).' strrep(bpmName,'_','.') 'S.Samples.samples.value;']);
            catch
                bpmS(t,:) = NaN(1,nSamples);
                nBadS = nBadS+1;
            end
        end
        if (useVertical)
            try
                bpmV(t,:) = eval(['CTFData(t).' bpmName 'V.Samples.samples.value;']);
            catch
                try
                    bpmV(t,:) = eval(['CTFData(t).' strrep(bpmName,'_','.') 'V.Samples.samples.value;']);
                catch
                    bpmV(t,:) = NaN(1,nSamples);
                    nBadV = nBadV+1;
                end
            end
        end
    end
    
    if (nBadH>0 || nBadS>0 || nBadV>0)
        fprintf('%s bad data: H - %d pulses, S - %d pulses, V - %d pulses\n',bpmName,nBadH,nBadS,nBadV);
    end

end

