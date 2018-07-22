function [CTFData] = saveMergedData( dataDir )
%saveMergedData Merges data from dataDir in to one array of structs and
%then saves it to the directory 'Merged' in dataDir.

    CTFData = mergeMatMonData(dataDir);

    mergedDir = [dataDir '/Merged'];

    fprintf('Saving data to %s...\n',mergedDir);
    
    if (~exist(mergedDir,'dir'))
        mkdir(mergedDir);
    end
    
    [dataBaseDir,~] = fileparts(dataDir);
    dataSetName = strrep(dataDir,[dataBaseDir '/'],'');
    
    saveName = [mergedDir '/' dataSetName '.mat'];
    
    
    varsInSpace = whos;
    sizeVars = sum([varsInSpace.bytes]);
    if (sizeVars>2e9) % split in to multiple files?
        error('file too big');
        %save(saveName,'CTFData','-v7.3');
    else % if smaller than 2GB use whatever the default is
        save(saveName,'CTFData');
    end

end

