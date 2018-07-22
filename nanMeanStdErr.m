function [datMean, datStd, datMeanErr, datStdErr ] = nanMeanStdErr(X, dim, squeezeOut)
%nanMeanStdErr Returns mean, standard deviation and their errors for data
%in X along dimension dim ignoring NaN values. If squeezeOut is true or not
%given, all outputs are squeezed to remove singleton dimensions.
%[datMean, datStd, datMeanErr, datStdErr ] = nanMeanStdErr(X, dim, squeezeOut)    
    goodVals = ~isnan(X);

    if (nargin<2 || isempty(dim))
        datMean = nanmean(X);
        datStd = nanstd(X);
        N = sum(goodVals);   
    else
        datMean = nanmean(X,dim);
        datStd = nanstd(X,0,dim);
        N = sum(goodVals,dim);
    end
    
    datMeanErr =  datStd./sqrt(N);
    datStdErr = datStd./sqrt(2*N);

    if (nargin<3 || isempty(squeezeOut) || squeezeOut==1)
        datMean = squeeze(datMean);
        datStd = squeeze(datStd);
        datMeanErr = squeeze(datMeanErr);
        datStdErr = squeeze(datStdErr);
    end
    
end

