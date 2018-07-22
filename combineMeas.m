function [ combMean,combErr ] = combineMeas( values, errors )
%combineMeas Takes measurements of same quantity with their errors and
%outputs best estimate of mean and error after combining the measurements.
% [ combMean,combErr ] = combineMeas( values, errors )
    % ignore NaN values
    goodMeas = (~isnan(values))&(~isnan(errors));
    if (sum(goodMeas)<1)
        warning('all values NaN')
        combMean=NaN;
        combErr=NaN;
        return;
    end
    values = values(goodMeas);
    errors = errors(goodMeas);

    % calculate mean, error
    sumInvErrSq = sum(1./(errors.^2));
    sumValOverErrSq = sum(values./(errors.^2) );
    
    combMean = sumValOverErrSq./sumInvErrSq;
    combErr = 1/sqrt(sumInvErrSq);

end

