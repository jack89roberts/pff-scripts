function [ processedData ] = loadProcessedData( dataSetNames )
%loadProcessedData Outputs a cell array of multiple processed data files
%   dataSetNames = cell array of names (not paths)
%   Assumes first 6 characters of dataSetNames are yyyyMM (e.g. 201412)
%   Assumes data sets are stored in yyyyMM folders in baseDir
%   Assumes processed data stored in processed folder, with same name as
%   data set, inside data set folder.
%   i.e. /home/jack/PhaseFeedforward/CTFData/<yyyyMM>/<dataSetName>/processed/<dataSetName>.mat
    baseDir = '/home/jack/PhaseFeedforward/CTFData/';
    processedDir = 'processed/';

    nDataSets = length(dataSetNames);
    processedData = cell(1,nDataSets);

    for i=1:nDataSets
        yearMonth = dataSetNames{i}(1:6);
        dataPath = [baseDir yearMonth '/' dataSetNames{i} '/' processedDir dataSetNames{i} '.mat'];
        tmpData = load(dataPath);
        processedData{i} = tmpData.processedData;
    end


end

