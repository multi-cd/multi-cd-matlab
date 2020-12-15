function [s_set, HS, more_output] = runSimulatedAnnealing(costfun, N, opts)
    % runSimulatedAnnealing: Runs simulated annealing at fixed lambda.
    % See demo2_MultiCD.m for a step-by-step demonstration.
    % --------------------------------------------------------------------

    % Copyright 2018-2020 Min Hyeok Kim & Ji Hyun Bak


    % === unpack options

    saOptions = opts.SA;
    mcmcOptions = opts.MCMC;
    quenchOptions = opts.quench;

    talkative = opts.talkative;
    save_output = opts.save_output;


    % ==== initial state generation

    % randomly select the number of clusters
    % (to allow sampling of highly clustered or highly fragmented solutions)
    K_init = randperm(N,1);

    % randomly assign each loci to a cluster
    s_init = randsample(K_init,N,true);
    s_init = renumber_clusters(s_init);

    K_init = max(s_init); % fix in case some domains were never assigned


    % ==== initial temperature decision

    T_init = initial_temp_decision(costfun,s_init);
    % disp(['initial temp. = ',num2str(T_init,'%1.1f')]);

    if save_output
        % save initial information
        save_initial_info(opts, s_init, costfun(s_init), T_init)
    end


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

    t_elap_total = 0;

    for iter = 1:maxIterSA

        % run MCMC sampling at this temperature
        tic;
        [s_list,HS_list,tau,term_status_MCMC] = runMCMC_fixedT(costfun,T,s_set,mcmcOptions);
        t_elap = toc;
        t_elap_total = t_elap_total + t_elap;

        if (term_status_MCMC == 2 && talkative)
            disp(' - premature end of MCMC. maxIter reached.')
        end
        quit_cond_SA = (term_status_MCMC<0);

        % choose best sample (cost minimizer)
        [HS,ipick] = min(HS_list);
        s_set = s_list(:,ipick);

        % store results from current step
        T_track(iter) = T;
        HS_track(iter,:) = [HS,mean(HS_list,'omitnan'),std(HS_list,'omitnan')];

        if (save_output && rem(iter,20) == 0)
            % save sampling output from this T
            save_samp_output_fixedT(opts, iter, HS_list, s_set, t_elap, T)
        end

        if talkative
            disp(['Step #',num2str(iter),...
                ', T=',num2str(T,'%1.1f'),', H=',num2str(HS,'%1.1f')]);
        end

        if quit_cond_SA
            if talkative
                disp('* single-domain solution reached -- force quitting SA.');
            end
            break
        end

        % update temperature
        T = T*c_cool; % cooling by a constant factor
        if (T<T_stop)
            break
        end

    end

    if talkative
        disp('End of simulated annealing.');
    end


    % === Final quenching by gradient descent

    % zero temperature quenching (gradient descent)
    tic;
    [s_set,HS,HS_list,term_status] = runQuench_zeroT(costfun,s_set,quenchOptions);
    if talkative
        if (term_status==1)
            disp(' - bad sign: descending too far. trial limit reached.');
            %warning('maybe SA was not enough. check temperature schedule.');
        elseif (term_status==2)
            disp('Minimum reached successfully.');
        end
    end

    t_elap = toc;
    t_elap_total = t_elap_total + t_elap;
    if save_output
        save_quench_output(opts, HS_list, s_set, t_elap_total);
    end

    more_output = struct('HS_list', HS_list);

end

function save_initial_info(opts, s_init, H_init, T_init)
    % save initial information
    save_folder_name = opts.out_path;
    save_name = [save_folder_name, '/initial_info.mat'];

    % match legacy namespace
    % eval(['save ',savename,' samplingtype r s_set HS_sum initial_temp comb_ratio '])
    initial_info.samplingtype = 'Single_mut';
    initial_info.r = opts.lambda;
    initial_info.s_set = s_init;
    initial_info.HS_sum = H_init;
    initial_info.initial_temp = T_init;
    initial_info.comb_ratio = 0.5; % no longer used
    save(save_name, '-struct', 'initial_info');

end

function save_samp_output_fixedT(opts, iter, HS_list, s_set, t_elap, T)
    save_folder_name = opts.out_path;
    save_name = [save_folder_name, '/results', num2str(iter), '.mat'];

    % match legacy namespace
    % eval(['save ',save_name,' HS_list s_set t_elap tmp'])
    output = [];
    output.HS_list = HS_list;
    output.s_set = s_set;
    output.t_elap = t_elap;
    output.tmp = T;
    save(save_name, '-struct', 'output');

end

function save_quench_output(opts, HS_list, s_set, t_elap_total)
    % NOTE: also save every 1000 trials as /ZeroT# ?
    save_folder_name = opts.out_path;
    save_name=[save_folder_name,'/ZeroTfinal'];

    % match legacy namespace
    % eval(['save ',save_name,' s_set time_cost HS_list']);
    output = [];
    output.s_set = s_set;
    output.HS_lisit = HS_list;
    output.time_cost = t_elap_total;
    save(save_name, '-struct', 'output');

end
