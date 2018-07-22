function CTFData = mergeMatMonData( dataDir )
% mergeMatMonData Takes all the .mat files (output files from
% MatlabMonitor) and creates an array CTFData that contains all the data
% structs.
%   One trigger per file, name of struc in file is tmpDataStruct
%   File order in CTFData array is the same as returned by dir/ls
%   (alphanumerical order?)
%
%   If there's some error the function loops through all the files and
%   only loads those with the longest structures (hopefully skips bad
%   data). However, this function is slow in general and if there's bad
%   data it's even slower.
    tic
    homeDir = pwd;
    cd(dataDir);
    fileList = dir('*.mat');
    nFiles = length(fileList);
    if (nFiles < 1)
        cd(homeDir);
        error(['There are no mat files in ' dataDir]);
    else
        fprintf(1,'Loading data from %d files...\n', nFiles);
    end

    
    [~,dataSetName,~] = fileparts(dataDir);
    dataYear = str2double(dataSetName(1:4));
    if (isnan(dataYear) || dataYear > 2020 || dataYear < 2012)
        warning('Unable to determine year of data from data path. Assuming data uses new MatLab JAPC libraries (from 2015)');
        newJAPC = 1;
    elseif (dataYear >= 2015)
        newJAPC = 1;
    else
        newJAPC = 0;
    end
    
    startIndex = 1;
    mergeSuccess = 0;
    while(startIndex<11)
        clear CTFData;
        
        try % try without removing bad files to start with as it's faster and we might be lucky
            % load files in to CTFData array
            for i=startIndex:nFiles
                fprintf(1,'File %d of %d: %s\n',i,nFiles,fileList(i).name);
                tmp = open(fileList(i).name);
                
                if (i==startIndex)
                    if (~newJAPC)
                        CTFData(i-startIndex+1) = tmp.tmpDataStruct;
                    else
                        CTFData(i-startIndex+1) = tmp.myDataStruct;
                    end
                else
                    % sometimes have errors when fields are saved in a
                    % different order, force the same order here.
                    if (~newJAPC)
                        tmp = orderfields(tmp.tmpDataStruct,CTFData(1));     
                    else
                        tmp = orderfields(tmp.myDataStruct,CTFData(1));
                    end
                    CTFData(i-startIndex+1) = tmp;
                end
                
            end
            
            mergeSuccess = 1;
            break;
            
        catch err % some error, probably some bad files
            startIndex = i+1;
            warning('File %s gave an error. Retrying with a start index of %d',fileList(i).name,i+1);
        end
    end
    
    if (mergeSuccess)
        fprintf(1,'Successfully loaded %d files! %d files were skipped.\n',length(startIndex:nFiles),startIndex-1);
    else
        fprintf(1,'Data merging failed with this error:\n');
        fprintf(1,'%s\n',err.message);
    end
    cd(homeDir);
    toc
end
