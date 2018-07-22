function [correlations,lags] = nanxcorr(A,B,maxLag,dim)
%NANXCORR Uses NANCOV to calculate the cross correlation between two
%signals when NaN is present.
%
% SYNOPSIS: correlations = nanxcorr(A,B,wt)
%           correlations = nanxcorr(A,B,wt,dim)
%
%
% Modiefied from acmartin. xies@mit.edu Jan 2012.
% Edited by Jack Roberts, June 2015.

if nargin<3
    maxLag = 10;
end
if nargin<4
    dim = 2;
end

if any(size(B) ~= size(A)), error('Input matrices must be the same size.'); end

lags = -maxLag:maxLag;

switch dim
    case 2
        [N,T] = size(A);
        correlations = nan(N,2*maxLag+1);
        
        for i = 1:N
            for t = lags
                signal = cat(1,A(i,max(1,1+t):min(T,T+t)), ...
                    B(i,max(1,1-t):min(T,T-t)));
                signal = signal';
                cov_mat = nancov(signal);
                variances = diag(cov_mat);
                corr = cov_mat./sqrt(variances*variances');
                correlations(i,t+maxLag+1) = corr(1,2);
                
            end
        end
    case 1
        [T,N] = size(A);
        correlations = nan(N,2*maxLag+1);
        
        for i = 1:N
            for t = lags
                signal = cat(2,A(max(1,1+t):min(T,T+t),i), ...
                    B(max(1,1-t):min(T,T-t),i));
                cov_mat = nancov(signal);
                variances = diag(cov_mat);
                corr = cov_mat./sqrt(variances*variances');
                correlations(i,t+maxLag+1) = corr(1,2);
                
            end
        end
        
        
    otherwise
        error('Unsupported dimension.');
end