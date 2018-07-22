baseDir = '/home/jack/PhaseFeedforward/CTFData/201412/';
dataSetListFile = '201412_DataSetCalibrations.txt'; % list of data sets with calibration file to use (comma separated)

%%
datFile = fopen([baseDir dataSetListFile]);
delimiter = ',';
formatSpec = '%s%s%[^\n\r]';
importedData = textscan(datFile, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(datFile);

dataDirs = importedData{:, 1};
calFiles = importedData{:, 2};

for i=1:length(dataDirs)
    fprintf(1,'-----------------------------------\n');
    fprintf('DATA SET %d OUT OF %d: %s\n',i,length(dataDirs), dataDirs{i});
    fprintf(1,'-----------------------------------\n');
    
    tmpDir = [baseDir dataDirs{i}];
    tmpCal = [baseDir 'frascatiCalibrations/frascatiCalibrationConstants_' calFiles{i}];
    
    processedData = processData(tmpDir,tmpCal);
    
    if (~isempty(processedData))
        fprintf(1,'Saving data from %s\n', dataDirs{i});
        saveDir = [baseDir dataDirs{i} '/processed/'];
        if (~exist(saveDir,'dir'))
            mkdir(saveDir);
        end

        save([saveDir dataDirs{i} '.mat'],'processedData');
    else
        fprintf(1,'All files in %s were bad, nothing saved.\n', dataDirs{i});
    end
    
end