function [ propertyNames, propertyValues ] = savePropValues( devsFile, propsFile, outputFile )
%savePropValues Gets
%   Saves properties in file propsFile for devices in file devsFile to
%   outputFile (single value properties only).
%    addpath('/nfs/cs-ccr-nfs6/vol29/share/ctf/data/ctfmod/MatLab/JavaCoInterface/');

    fileIDDevs = fopen(devsFile);
    devices = textscan(fileIDDevs,'%s');
    devices = devices{1};
    fclose(fileIDDevs);
    
    fileIDProps = fopen(propsFile);
    properties = textscan(fileIDProps,'%s');
    properties = properties{1};
    fclose(fileIDProps);

    nDevices = length(devices);
    nProperties = length(properties);
    nPropsToGet = nProperties*nDevices;
    
    splitProps = cell(nProperties,3);
    for p=1:nProperties
        %splitProps(p,:) = strsplit(properties{p},{'/','#',' ',','});      
        splitProps(p,:) = regexpi(properties{p},'(/|#|,)','split');
    end
    
    propertyNames = cell(nPropsToGet,1);
    propertyValues = NaN*ones(nPropsToGet,1);
    outFileID = fopen(outputFile,'w');
    for d=1:nDevices
        for p=1:nProperties
            tmpDevName = [devices{d} splitProps{p,1}];
            propertyNames{p + (nProperties*(d-1))} = [devices{d} properties{p}];
            %propertyValues(p + (nProperties*(d-1))) = JGetCoValueFesa(tmpDevName, splitProps{p,2}, splitProps{p,3});
            propertyValues(p + (nProperties*(d-1))) = matlabJapc.staticGet('SCT.USER.ALL',tmpDevName, splitProps{p,2}, splitProps{p,3});
            fprintf(outFileID,'%s, %f\n',[devices{d} splitProps{p,1} ', ' splitProps{p,2} ', ' splitProps{p,3}],propertyValues(p + (nProperties*(d-1))));
        end
    end
    
    fclose(outFileID);
    
end

