function R = myautocorr(y,maxlag)
% MYAUTOCORR: Computes autocorrelation for single-column data, up to maxlag
% 
% Inputs:  y : single data series, in a [N 1] vector
%          maxlag : total number of lags, single integer (used as P)
% Output:  R : autocorrelation, [maxlag 1] vector, starting at lag 1
% ------------------------------------------------------------------------

% 2018 Ji Hyun Bak

%% unpack input

if(maxlag>numel(y))
    error('myautocorr: maxlag should be less than data length.');
end

dy = y(:) - mean(y(:)); % force column vector

%% set up auxiliary functions

% variance and covariance
getvar = @(dy) sum(dy.*dy)/numel(dy);
getcov = @(dy1,dy2) sum(dy1.*dy2)/numel(dy1);

% autocorrelation at fixed lag
getcorr = @(dy,lag) getcov(dy((lag+1):end),dy(1:(end-lag)))/getvar(dy);

%% compute autocorrelation

R = zeros(maxlag,1);
for lag = 1:maxlag
   R(lag) = getcorr(dy,lag);
end

end
