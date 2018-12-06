function [myH,myE,myK] = HS_calculation_all(Cmat,s_set,lambda)
% HS_calculation_all: Calculates the cost function 
% given a domain solution and a data matrix,
% by summing over all domains in the solution.
% 
% INPUT:
%   - Cmat: input correlation matrix, [N N] symmetric matrix
%   - s_set: domain solution of interest, N-vector
%   - lambda: the Multi-CD parameter, a single number
% 
% OUTPUT:
%   - myH: the total cost function, or the "hamiltonian",
%          corresponding to the (negative log) posterior in Bayesian setup
%   - myE: cost function associated with the goodness of clustering,
%          or the "energy", corresponding to the (negative log) likelihood
%   - myK: cost function associated with the simplicity of solution,
%          corresponding to the (negative log) prior
% ------------------------------------------------------------------------

% Copyright 2018 Min Hyeok Kim & Ji Hyun Bak


% calculate the hamiltonian value for cluster state 
myE = getE(s_set,Cmat);

% calculate the regularization value 
myK = getK_effNumClust(s_set);
     
% total hamiltonian
myH = myE + lambda*myK;

end

% ========================================================================

function myE = getE(s_set,Cmat)
% The likelihood of a clustering state given correlation matrix
% [ref: Giada and Marsili (2001), PRE]
% and report the negative log likelihood, or the "energy".

Kmax = max(s_set); % number of clusters

Ek_list = zeros(Kmax,1);

for k = 1:Kmax % for each cluster
    
    ix = (s_set==k); % index elements in cluster

    Cmat_sub = Cmat(ix,ix); 
    ck = sum(Cmat_sub(:)); % sum of correlation matrix elements in cluster
    nk = sum(ix); % number of elements (size of cluster)
    
    if(nk<=1)
        Hk = 0;
    elseif(ck == nk^2) % in case of all elements just "1"
        ck = nk*1 + 2*nchoosek(nk,2)*0.99; % set large enough value 
        Hk = log(ck/nk) + (nk-1)*log((nk^2-ck)/(nk^2-nk));
    elseif(ck<=nk) % in case anti-correlation dominant 
        Hk = 0;
    else
        Hk = log(ck/nk) + (nk-1)*log((nk^2-ck)/(nk^2-nk));
    end

    Ek_list(k) = Hk;
end

myE = sum(Ek_list);

end

function myK = getK_effNumClust(s_set)
% The sparsity prior, penalizing fragmented solutions;
% the value is larger for more fragmented solutions.
% Quantifies effective number of clusters (exp(entropy))

% get cluster sizes
nk_list = zeros(max(s_set),1); % cluster number list 
for k = 1:max(s_set)
    nk_list(k) = sum(s_set==k); % number of elements in cluster s
end
N = sum(nk_list);

% effective number of clusters [*NOTE: log(number) ~ entropy ]
myK = exp(sum(-log(nk_list./N).*nk_list./N));

end
