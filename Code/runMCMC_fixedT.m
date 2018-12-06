function [s_list,HS_list,tau,term_status] = ...
    runMCMC_fixedT(costfun,T,s_init,mcmcOptions)
% runMCMC_fixedT: Runs MCMC sampling at a fixed temperature T.
% 
% INPUT:
%   - costfun: a function handle for the cost function, with one argument 
%              that is the domain state vector -- e.g., costfun(s)
%   - T : the tempering factor / temperature (single number)
%   - s_init: initial state (N-vector)
%   - mcmcOptions: a struct variable that contains MCMC options
% 
% OUTPUT:
%   - s_list: the trajectory of domain states sampled in the MCMC chain
%             [N M] array, where M is the final # samples in MCMC chain
%             each column is a domain state
%   - HS_list: the trajectory of cost function values in MCMC chain
%             [M 1] vector
%   - tau: the relaxation time of the chain
%   - term_status: reports by what termination condition the program exits
% ------------------------------------------------------------------------

% Copyright 2018 Min Hyeok Kim & Ji Hyun Bak

%% initialize

% initial state
s_set = s_init;
HS = costfun(s_set);

% unpack MCMC options
minIterMCMC = mcmcOptions.minIterMCMC;
maxIterMCMC = mcmcOptions.maxIterMCMC;
numSampAtEq = mcmcOptions.numSampAtEq; % final sampling after reaching equilibrium
if(isfield(mcmcOptions,'checkEvery'))
    checkEvery = mcmcOptions.checkEvery; % set frequency to calculate tau
else
    checkEvery = floor(maxIterMCMC/20); % default
end

% termination status
term_status = 0; % 0 means running


%% iterate

cnt = 0;  % step index within a MCMC chain
tau = Inf; % initialize just in case of early termination

% track and store the MCMC chain
s_list = zeros(numel(s_init),maxIterMCMC);
HS_list = zeros(maxIterMCMC,1);

while 1
    
    cnt=cnt+1;
    
    % propose next move
    s_set_cdd = single_mutation(s_set); % single mutation
    s_set_cdd = renumber_clusters(s_set_cdd); % renumber domain indices
    
    HSpos = costfun(s_set_cdd); % cost function at the proposed move
    delta_HS = HSpos-HS; % the "energy" difference
    acceptance_prob = exp(-delta_HS/T); % Metropolis-Hastings
    
    % accept or reject the move
    if(acceptance_prob>rand)
        s_set = s_set_cdd; % move accepted: update cluster state
        HS = HSpos;
    else
        % move rejected: keep previous state
    end
    s_list(:,cnt) = s_set; % store current sample
    HS_list(cnt) = HS; % store current H value
    
    
    % ==== Check for the stopping conditions: =========================
    
    % update the relaxation time every once in a while
    if(rem(cnt,checkEvery)==0 && cnt>=minIterMCMC)
        autoC = myautocorr(HS_list(1:cnt),cnt-1); % autocorrelation function
        diff1 = abs(autoC-exp(-1));
        tau = find(diff1<1e-3,1,'first'); % relaxation time tau
        if isempty(tau)
            tau=inf;
        end
        
        % Stopping condition 1: enough samples at equilibrium
        if(cnt>5*tau+numSampAtEq) 
            term_status = 1;
            break;
        end
    end
	
    % Stopping condition 2: iteration takes too long
    if(cnt>maxIterMCMC) 
        disp(' - premature end of MCMC. maxIter reached.')
        term_status = 2;
        break;
    end
    
    % Stopping condition -1: reports a single-domain state
    if max(s_set)==1 
        term_status = -1; % (negative value: passes a stop signal for the entire SA)
        break;
    end
    % =================================================================
    
end

% trim lists
idx_keep = (cnt-numSampAtEq)+(1:numSampAtEq);
s_list = s_list(:,idx_keep);
HS_list = HS_list(idx_keep,:);

end
