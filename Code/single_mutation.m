function svec_mut = single_mutation(svec)
% SINGLE_MUTATION: Makes a single mutation in the domain state space.
% 
% Input:    svec (original domain state)
%           : should be an N-vector, that spans a consecutive integer range
% Output:   svec_mut (new domain state after a single mutation)
%           : of the same dimension as input
% -------------------------------------------------------------------------

% Copyright 2018 Min Hyeok Kim & Ji Hyun Bak


%% unpack input

% total number of loci
N = numel(svec);

% max index of existing domains (Kmax domains total)
Kmax = max(svec); 

% check that svec spans all domains 1:Kmax
if(numel(unique(svec))~=Kmax)
    error('domain indices should span a consecutive range of integers 1:K');
end


%% single mutation

% ===== choose a locus

% choose one locus at random (mutation will be made at this locus)
j = randsample(N,1);
s_j = svec(j); % domain index of this locus


% ===== make mutation

% choose the target domain
s_i = randsample(1:Kmax,1); % choose a target domain index
if(s_i==s_j)
    s_i = Kmax+1; % if self (s_j), add a new domain (Kmax+1)
end

% apply mutation
svec_mut = svec;
svec_mut(j) = s_i;

end
