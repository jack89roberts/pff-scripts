function [ bpmNames, bpmH_MeansAtCursor, bpmH_Samples, bpmS_MeansAtCursor, bpmS_Samples ] = extractSectionBPMsFromCTFData( sectionPrefixes, CTFData )
%extractBPMsFromCTFData Get BPM signals from merged CTFData
%
%   sectionPrefixes - MUST BE CELL ARRAY (even if only one prefix).
%   Possible values CC,CB,CM and CT. e.g. sectionPrefixes = {'CC'} or
%   {'CC','CB'}. Data from all BPMs with given prefixes will be extracted.
%
    for p=1:length(sectionPrefixes)

        devsFile = ['/nfs/cs-ccr-nfs6/vol29/share/ctf/data/ctfmod/MatLab/PhaseFeedforward/devices/bpm' sectionPrefixes{p} '.devs'];

        fileIDDevs = fopen(devsFile);
        devices = textscan(fileIDDevs,'%s');
        devices = devices{1};
        fclose(fileIDDevs);
        
        if p==1
            bpmNames = devices;
        else
            bpmNames = [bpmNames; devices];
        end
        
    end
    
    nBPMs = length(bpmNames);
    nPulses = length(CTFData);
    nSamples = 768;
    
    bpmH_MeansAtCursor = NaN*ones(nBPMs,nPulses);
    bpmH_Samples = NaN*ones(nBPMs,nPulses,nSamples);
    bpmS_MeansAtCursor = NaN*ones(nBPMs,nPulses);
    bpmS_Samples = NaN*ones(nBPMs,nPulses,nSamples);
    for bpm=1:length(bpmNames)
         [  tmpH_MeanAtCursor, tmpH_Samples, tmpS_MeanAtCursor, tmpS_Samples ] = extractBPMFromCTFData( bpmNames{bpm}, CTFData );

         bpmH_MeansAtCursor(bpm,:) = tmpH_MeanAtCursor;
         bpmH_Samples(bpm,:,1:length(tmpH_Samples)) = tmpH_Samples;
         bpmS_MeansAtCursor(bpm,:) = tmpS_MeanAtCursor;
         bpmS_Samples(bpm,:,1:length(tmpS_Samples)) = tmpS_Samples;

    end
    
end

