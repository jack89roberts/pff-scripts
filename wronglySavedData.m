dataDir = '/nfs/cs-ccr-nfs6/vol29/share/ctf/data/ctfmod/MatLab/PhaseFeedforward/data/20141204_1617_FF_NewFMW_Gate365to465_GainK1_-0_GainK2_0';

tmpList = dir([dataDir '/*.mat']);
fileNames = cell(1,length(tmpList));
fileDateStamps = NaN*ones(1,length(tmpList));
beamEnabled = NaN*ones(1,length(tmpList));
for i=1:length(tmpList)
    fprintf(1, '%d out of %d\n', i, length(tmpList));
    
    fileNames{i} = tmpList(i).name;
    tmpNums = regexp(fileNames{i}, '\.', 'split');
    tmpNumString = '';
    for j=1:(length(tmpNums)-1)
        tmpNumString = [tmpNumString tmpNums{j}];
    end
    
    fileDateStamps(i) = str2double(tmpNumString);
    
    try     
        load([dataDir '/' tmpList(i).name])
        tmpTransmission = min(tmpDataStruct.CT.SVBPI0105S.Samples.samples.value);
        if min(tmpTransmission) < -2
            beamEnabled(i) = 1;
        else
            beamEnabled(i) = 0;
        end
    catch
        beamEnabled(i) = 0;
    end
end

subplot(1,2,1)
plot(fileDateStamps,'b.');
title('time stamps');
subplot(1,2,2)
plot(beamEnabled,'r');
title('beam enabled');

plot(diff(fileDateStamps))
plot(fileDateStamps,ones(1,length(fileDateStamps)),'.')