function [s_set, HS, more_output] = runSimulatedAnnealing(costfun, N, opts)
    % Function to run simulated annealing at fixed lambda.
    % See demo2_MultiCD.m for a step-by-step demonstration.
    
    % === unpack options
    
    saOptions = opts.SA;
    mcmcOptions = opts.MCMC;
    quenchOptions = opts.quench;
    
    talkative = opts.talkative;
    
    
    % ==== initial state generation
    
    % randomly select the number of clusters
    % (to allow sampling of highly clustered or highly fragmented solutions)
    K_init = randperm(N,1);
    
    % randomly assign each loci to a cluster
    s_init = randsample(K_init,N,true);
    s_init = renumber_clusters(s_init);
    
    % K_init = max(s_init); % fix in case some domains were never assigned
    
    
    % ==== initial temperature decision
    
    T_init = initial_temp_decision(costfun,s_init);
    
    
    % ==== simulated annealing
    
    % initialize
    T = T_init;
    s_set = s_init;
    if talkative
        disp(['Initial T=',num2str(T,'%1.1f')]);
    end
    
    c_cool = saOptions.c_cool;
    maxIterSA = saOptions.maxIterSA;
    T_stop = saOptions.T_stop;
    
    HS_track = NaN(maxIterSA,3); % [best mean std]
    T_track = NaN(maxIterSA,3);
    
    for r = 1:maxIterSA
        
        % run MCMC sampling at this temperature
        [s_list,HS_list,tau,term_status_MCMC] = runMCMC_fixedT(costfun,T,s_set,mcmcOptions);

        quit_cond_SA = (term_status_MCMC<0);
        
        % choose best sample (cost minimizer)
        [HS,ipick] = min(HS_list);
        s_set = s_list(:,ipick);
        
        % store results from current step
        T_track(r) = T;
        HS_track(r,:) = [HS,mean(HS_list,'omitnan'),std(HS_list,'omitnan')];
        
        if talkative
            disp(['Step #',num2str(r),...
                ', T=',num2str(T,'%1.1f'),', H=',num2str(HS,'%1.1f')]);
        end
        
        if(quit_cond_SA)
            if talkative
                disp('* single-domain solution reached -- force quitting SA.');
            end
            break
        end
        
        % update temperature
        T = T*c_cool; % cooling by a constant factor
        if(T<T_stop)
            break
        end
        
    end
    
    if talkative
        disp('End of simulated annealing.');
    end
    
    
    % === Final quenching by gradient descent
    
    % zero temperature quenching (gradient descent)
    [s_set,HS,HS_list,term_status] = runQuench_zeroT(costfun,s_set,quenchOptions);
    if talkative
        if(term_status==1)
            disp(' - bad sign: descending too far. trial limit reached.');
            %warning('maybe SA was not enough. check temperature schedule.');
        elseif(term_status==2)
            disp('Minimum reached successfully.');
        end
    end
    
    more_output = struct('HS_list', HS_list);
    
end
