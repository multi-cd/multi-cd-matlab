function svec_new = renumber_clusters(svec)
% RENUMBER_CLUSTERS: Renumbers domain (cluster) indices in a state vector
% such that they constitute a consecutive range of integers from 1 to Kmax, 
% where Kmax is the total number of domains.
% 
% Input:    svec (original domain state) -- should be an N-vector
% Output:   svec_new -- of the same dimension as the input, 
%                       but possibly with re-numbered domain indices
%                       from 1 to Kmax (totla number of domains)
% ------------------------------------------------------------------------

% Copyright 2018 Min Hyeok Kim & Ji Hyun Bak

[~,~,svec_new] = unique(svec,'stable');

end
