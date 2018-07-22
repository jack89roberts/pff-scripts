function [ signalList, monitorIdxs ] = combinePropertyLists( varargin )
%combinePropertyLists Takes cell arrays of property lists 
% (e.g. from makePropertyListFromFile) and outputs
% one combined list of properties and monitor idxs for matlab monitor.
%   varargin is the property lists to combine, e.g. bpmList, bprList
%   items from each list are assigned a different monitorIdx for the matlab monitor
    
    nLists = length(varargin);
    listLengths = zeros(1,nLists);
    for l=1:nLists
        listLengths(l) = length(varargin{l});
    end
    nPropertiesCombined = sum(listLengths);
    
    signalList = cat(2,varargin{:}); % combine cell arrays in to one
    monitorIdxs = zeros(1,nPropertiesCombined);
    
    monitorIdxs(1:listLengths(1)) = 1;
    startInd = listLengths(1) + 1;
    idx = 2;
    for l=2:nLists
        monitorIdxs( startInd : (startInd+listLengths(l)-1) ) = idx;
        startInd = startInd+listLengths(l);
        idx = idx + 1;
    end

end

