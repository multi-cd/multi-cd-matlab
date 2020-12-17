function [s_set, HS, more_output] = runSimulatedAnnealing(costfun, N, opts)
    % runSimulatedAnnealing: Runs simulated annealing at fixed lambda.
    % See demo2_MultiCD.m for a step-by-step demonstration.
    % --------------------------------------------------------------------

    % Copyright 2018-2020 Min Hyeok Kim & Ji Hyun Bak


    % === unpack options

    saOptions = opts.SA;
    mcmcOptions = opts.MCMC;
    quenchOptions = opts.quench;
    
    % optional settings
    talkative = getFromStruct(opts, 'talkative', true);
    save_output = getFromStruct(opts, 'save_output', false);
    if ~isfield(opts, 'overwrite')
        opts.overwrite = false;
    end
    
    more_output = [];
    more_output.term_message = ''; % empty means normal


    % ==== initialization
    
    const_a = getFromStruct(saOptions, 'const_a', 0.5);
    T_init_range = getFromStruct(saOptions, 'T_init_range', []);
    
    [s_set, ~] = initial_state_generation(N);
    T = initial_temp_decision(costfun, s_set, const_a, T_init_range);
    if talkative
        disp(['Initial T=',num2str(T,'%1.1f')]);
    end
    if save_output
        % save initial information
        HS = costfun(s_set);
        flag = save_initial_info(opts, s_set, HS, T);
        if flag
            msg = 'failed to save initial_info';
            disp(msg);
            more_output.term_message = msg;
            return;
        end
    end


    % ==== simulated annealing

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
        more_output.samp_term_status = term_status_MCMC;

        % choose best sample (cost minimizer)
        [HS,ipick] = min(HS_list);
        s_set = s_list(:,ipick);

        % store results from current step
        T_track(iter) = T;
        HS_track(iter,:) = [HS,mean(HS_list,'omitnan'),std(HS_list,'omitnan')];

        if (save_output && rem(iter,20) == 0)
            % save sampling output from this T
            flag = save_samp_output_fixedT(opts, iter, HS_list, s_set, t_elap, T);
            if flag
                msg = ['failed to save sampling output at T=', num2str(T)];
                disp(msg);
                more_output.term_message = msg;
                return;
            end
        end

        if talkative
            disp(['Step #',num2str(iter),...
                ', T=',num2str(T,'%1.1f'),', H=',num2str(HS,'%1.1f'),...
                ' (',num2str(t_elap,'%.1f'),'s)']);
        end

        if quit_cond_SA
            msg = '* single-domain solution reached -- force quitting SA.';
            if talkative
                disp(msg);
            end
            more_output.term_message = msg;
            break;
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
    t_elap = toc;
    t_elap_total = t_elap_total + t_elap;
    
    if (term_status==1)
        if talkative
            disp(' - bad sign: descending too far. trial limit reached.');
        end
        %warning('maybe SA was not enough. check temperature schedule.');
    elseif (term_status==2)
        if talkative
            disp('Minimum reached successfully.');
        end
    end
    more_output.quench_term_status = term_status;
    
    if save_output
        flag = save_quench_output(opts, HS_list, s_set, t_elap_total);
        if flag
            msg = 'failed to save quench output';
            disp(msg);
            more_output.term_message = msg;
            return;
        end
    end
    
    more_output.HS_list = HS_list;

end

function flag = save_initial_info(opts, s_init, H_init, T_init)
    % save initial information
    save_folder_name = opts.out_path;
    % save_name = fullfile(save_folder_name, 'initial_info.mat');
    save_name = fullfile(save_folder_name, 'SA_init.mat');
    if ~opts.overwrite && exist(save_name, 'file')
        warning(['file already exists: ', save_name]);
        flag = 1;
        return
    end

    % match legacy namespace
    % eval(['save ',savename,' samplingtype r s_set HS_sum initial_temp comb_ratio '])
    initial_info.samplingtype = 'Single_mut';
    initial_info.r = opts.lambda;
    initial_info.s_set = s_init;
    initial_info.HS_sum = H_init;
    initial_info.initial_temp = T_init;
    initial_info.comb_ratio = 0.5; % no longer used
    save(save_name, '-struct', 'initial_info');
    flag = 0;

end

function flag = save_samp_output_fixedT(opts, iter, HS_list, s_set, t_elap, T)
    save_folder_name = opts.out_path;
    % save_name = fullfile(save_folder_name, ['results', num2str(iter), '.mat']);
    save_name = fullfile(save_folder_name, ['SA_it', num2str(iter, '%03d'), '.mat']);
    if ~opts.overwrite && exist(save_name, 'file')
        warning(['file already exists: ', save_name]);
        flag = 1;
        return
    end

    % match legacy namespace
    % eval(['save ',save_name,' HS_list s_set t_elap tmp'])
    output = [];
    output.HS_list = HS_list;
    output.s_set = s_set;
    output.t_elap = t_elap;
    output.tmp = T;
    save(save_name, '-struct', 'output');
    flag = 0;

end

function flag = save_quench_output(opts, HS_list, s_set, t_elap_total)
    % NOTE: also save every 1000 trials as /ZeroT# ?
    save_folder_name = opts.out_path;
    % save_name = fullfile(save_folder_name, 'ZeroTfinal.mat');
    save_name = fullfile(save_folder_name, 'SA_final.mat');
    if ~opts.overwrite && exist(save_name, 'file')
        warning(['file already exists: ', save_name]);
        flag = 1;
        return
    end

    % match legacy namespace
    % eval(['save ',save_name,' s_set time_cost HS_list']);
    output = [];
    output.s_set = s_set;
    output.HS_list = HS_list;
    output.time_cost = t_elap_total;
    save(save_name, '-struct', 'output');
    flag = 0;

end
