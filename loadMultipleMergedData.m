function [ CTFData, calStruct, dataDir ] = loadMultipleMergedData( dataSetNames, frascatiCalTimeStamps )
%loadMultipleMergedData Loads data from multiple data sets 
%   dataSetName and frascatiCalStamp must have same length.
%   outputs are cell arrays with each entry same as loadMergedData output
    tic;
    
    nDataSets = length(dataSetNames);
    
    CTFData = cell(1,nDataSets);
    calStruct = cell(1,nDataSets);
    dataDir = cell(1,nDataSets);
    
    fprintf('Loading data from %d data sets...\n',nDataSets);
    
    for i=1:nDataSets
        fprintf('----------------------------------------------------------------\n');
        fprintf('Data set %d of %d: %s\n',i,nDataSets,dataSetNames{i});
        fprintf('----------------------------------------------------------------\n');

        if (nargin<2 || isempty(frascatiCalTimeStamps))
            [tmpCTFData,tmpCalStruct,tmpDataDir] = loadMergedData(dataSetNames{i}); 
        else
            [tmpCTFData,tmpCalStruct,tmpDataDir] = loadMergedData(dataSetNames{i},frascatiCalTimeStamps{i});
        end
        
        CTFData{i} = tmpCTFData;
        calStruct{i} = tmpCalStruct;
        dataDir{i} = tmpDataDir;
        
    end
    
    fprintf('----------------------------------------------------------------\n');
    fprintf('Finished!\n');
    toc;
end

