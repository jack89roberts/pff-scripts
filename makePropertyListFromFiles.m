function propList = makePropertyListFromFiles( devsFile, propsFile )
% makePropertyListFromFiles Takes a file containing a list of device names
% and a file containing a list of property names and generates a list of
% device properties to subscribe to.
%   devsFile and propsFile are file name strings
    
    % changed to use textscan instead of importdata which is very slow
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
    
    nDevicePropertiesToSubscribe = nDevices*nProperties;
    
    propList = cell(1,nDevicePropertiesToSubscribe);
    
    for d=1:nDevices
        for p=1:nProperties
            propList{p + (nProperties*(d-1))} = [devices{d} properties{p}];
        end
    end
        
end

