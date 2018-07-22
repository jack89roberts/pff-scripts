baseDir = '/home/jack/PhaseFeedforward/CTFData/';
fileList = fopen([baseDir 'frascatiCalibrationFiles.txt']); % list of frascati calibration constant files
frascatiCalFiles = textscan(fileList,'%s');
fclose(fileList);

frascatiCalFiles = frascatiCalFiles{1};
nCalFiles = length(frascatiCalFiles);
dates = cell(nCalFiles,1);
times = cell(nCalFiles,1);

%  generate list of expected .mat files from names of frascati calibration
%  constant files. Also extract dates and times.
frascatiMatFiles = strrep(frascatiCalFiles,'Constants','');
for i=1:nCalFiles
    frascatiMatFiles{i} = [frascatiMatFiles{i} '.mat'];
    
    tmpName = regexp(frascatiCalFiles{i},'/','split');
    if length(tmpName)>1 % found a '/', last element should give name without leading directory paths
        tmpName = tmpName{length(tmpName)};
    end
    tmpName = regexp(tmpName,'_','split');
    dates{i} = tmpName{2};
    times{i} = tmpName{3};
end
dateTime = strcat(dates,times); % join together dates and times for plot axis
dateTime = str2double(dateTime);

calibrationConstants = NaN*ones(nCalFiles,3,4);
useMixerOverSqrtDiode = NaN*ones(nCalFiles,1);

comments = cell(nCalFiles,1);

for i=1:nCalFiles
    [tmpCalConst,tmpUseMixDio] = loadFrascatiCalibrationConstants(sprintf('%s%s',baseDir,frascatiCalFiles{i}));
    calibrationConstants(i,:,:) = tmpCalConst;
    useMixerOverSqrtDiode(i) = str2double(tmpUseMixDio);
    
    tmpMatFile = load(sprintf('%s%s',baseDir,frascatiMatFiles{i}));
    
    comments{i} = tmpMatFile.dataSetComment;
    
    fprintf('FILE %d: %s, %s. Comment: %s\n',i,dates{i},times{i},comments{i});
end

% figure;
% plot(dateTime,squeeze(calibrationConstants(:,1,1)),'-o');
% hold all;
% plot(dateTime,squeeze(calibrationConstants(:,2,1)),'-o');
% plot(dateTime,squeeze(calibrationConstants(:,3,1)),'-o');
% legend('Mon1','Mon2','Mon3');
% title('amplitude')
% 
% figure;
% plot(dateTime,squeeze(calibrationConstants(:,1,4)),'-o');
% hold all;
% plot(dateTime,squeeze(calibrationConstants(:,2,4)),'-o');
% plot(dateTime,squeeze(calibrationConstants(:,3,4)),'-o');
% legend('Mon1','Mon2','Mon3');
% title('offset')
% 
% figure;
% plot(dateTime,useMixerOverSqrtDiode,'-o');
% title('useMixerOverSqrtDiode');

figure;
plot(squeeze(calibrationConstants(:,1,1)),'-o');
hold all;
plot(squeeze(calibrationConstants(:,2,1)),'-o');
plot(squeeze(calibrationConstants(:,3,1)),'-o');
legend('Mon1','Mon2','Mon3');
title('amplitude')

figure;
plot(squeeze(calibrationConstants(:,1,4)),'-o');
hold all;
plot(squeeze(calibrationConstants(:,2,4)),'-o');
plot(squeeze(calibrationConstants(:,3,4)),'-o');
legend('Mon1','Mon2','Mon3');
title('offset')

figure;
plot(useMixerOverSqrtDiode,'-o');
title('useMixerOverSqrtDiode');