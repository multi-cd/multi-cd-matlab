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

    % Copyright 2018-2020 Min Hyeok Kim & Ji Hyun Bak

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


    %% iterate

    cnt = 0;  % step index within a MCMC chain
    tau = Inf; % initialize just in case of early termination

    % track and store the MCMC chain
    s_list = zeros(numel(s_init),maxIterMCMC);
    HS_list = zeros(maxIterMCMC,1);
    
    term_status = 0; % termination status (0 means running)
    while term_status == 0
        cnt=cnt+1;
        
        % propose next move
        s_set_cdd = single_mutation(s_set); % single mutation
        s_set_cdd = renumber_clusters(s_set_cdd); % renumber domain indices
        
        % Metropolis-Hastings
        HSpos = costfun(s_set_cdd); % cost function at the proposed move
        delta_HS = HSpos-HS; % the "energy" difference
        acceptance_prob = exp(-delta_HS/T);
        if (acceptance_prob>rand)
            % move accepted: update cluster state
            s_set = s_set_cdd;
            HS = HSpos;
        else
            % move rejected: keep previous state
        end
        s_list(:,cnt) = s_set; % store current sample
        HS_list(cnt) = HS; % store current H value
        
        % update the relaxation time every once in a while
        if (rem(cnt,checkEvery)==0 && cnt>=minIterMCMC)
            tau = update_relaxation_time(HS_list, cnt);
        end
        
        % check for the stopping conditions
        K = max(s_set); % number of clusters in the current state
        term_status = check_stopping_conditions(cnt, tau, K, ...
                                                maxIterMCMC, numSampAtEq);
        % if (term_status ~= 0)
        %     break;
        % end
        
    end

    % trim lists
    idx_keep = (cnt-numSampAtEq)+(1:numSampAtEq);
    s_list = s_list(:,idx_keep);
    HS_list = HS_list(idx_keep,:);

end

function tau = update_relaxation_time(HS_list, cnt)
    % thresholds
    C_thresh = exp(-1);
    epsilon = 1e-3;
    
    % calculate the relaxation time tau
    autoC = myautocorr(HS_list(1:cnt),cnt-1); % autocorrelation function
    diff1 = abs(autoC - C_thresh);
    tau = find(diff1 < epsilon, 1, 'first');
    if isempty(tau)
        tau = inf;
    end
end

function term_status = check_stopping_conditions(cnt, tau, K, maxIterMCMC, numSampAtEq)
    % check for stopping conditions.
    %
    % term_status (int): termination status
    % 0: running normally
    % positive value: end this MCMC chain but continue the simulated annealing
    %   1: normal termination (enough samples at equilibrium)
    %   2: potentially premature termination (reached maxIterMCMC)
    % negative value: passes a stop signal for the entire SA
    %   -1: stuck in a single-domain solution (wouldn't escape this state)
    % -------------------------------------------------------------------------

    % Stopping condition 1: enough samples at equilibrium
    % NOTE: this was previously checked when and only when tau was re-calculated
    % if(rem(cnt,checkEvery)==0 && cnt>=minIterMCMC)
    if (cnt > 5 * tau + numSampAtEq) 
        term_status = 1;
        return
    end
    % end
    
    % Stopping condition 2: iteration takes too long
    % NOTE: the inequality was previously a ">"; now fixed to a ">="
    if (cnt >= maxIterMCMC)
        % disp(' - premature end of MCMC. maxIter reached.')
        term_status = 2;
        return
    end
    
    % Stopping condition -1: reports a single-domain state
    if (K == 1)
        term_status = -1;
        return
    end
    
    % otherwise, keep running
    term_status = 0;
end
