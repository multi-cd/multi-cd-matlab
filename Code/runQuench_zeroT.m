function [s_set,HS,HS_list,term_status] = runQuench_zeroT(costfun,s_set_start,quenchOptions)
    % runQuench_zeroT: Runs zero-temperature "quenching" to quickly reach the 
    % nearest local mininum of the energy landscape being sampled
    % 
    % INPUT:
    %   - costfun: a function handle for the cost function, with one argument 
    %              that is the domain state vector -- e.g., costfun(s)
    %   - s_set_start: initial state (N-vector)
    %   - quenchOptions: a struct variable that contains options
    % 
    % OUTPUT:
    %   - s_set: the final state (N-vector)
    %   - HS: the final value of the cost function (the "energy")
    %   - HS_list: the trajectory of cost function values during quenching 
    %             [M+1 1] vector, where M is the total number of descents made
    %             first element stores the initial H value
    %   - term_status: reports by what termination condition the program exits
    % ------------------------------------------------------------------------

    % Copyright 2018 Min Hyeok Kim & Ji Hyun Bak


    % unpack input options
    cnt_acp_limit = quenchOptions.cnt_acp_limit;
    
    % track descending H values
    HS_start = costfun(s_set_start);
    HS_list = [HS_start; NaN(cnt_acp_limit, 1)]; % prepend initial value
    
    % initialize
    s_set = s_set_start;
    cnt_acp = 0; % initialize acceptance counter
    term_status = 0; % termination status (0 means running)
    
    while (term_status == 0)
    
        % find a descent along the H landscape
        [s_set, HS, found_descent] = next_descent(costfun, s_set);
        
        % Stopping condition 2: no more descents after "enough" trials
        % NOTE: this is the "normal" termination, should have been condition 1...
        if ~found_descent
            term_status = 2;
            HS_list = HS_list(1:1+cnt_acp); % only return values up to this point
            break;
        end
        
        % report accepted H value
        cnt_acp = cnt_acp + 1;
        HS_list(1+cnt_acp) = HS; % with 1-row offset for the initial value
        
        % Stopping condition 1: too many descents (too far from local min)
        if cnt_acp >= cnt_acp_limit
            term_status = 1; % a bad sign (premature end of program)
            % break;
        end
        
    end
    
end

function [s_set, HS, found_descent] = next_descent(costfun, s_set)
    % sample the neighborhood and take the first descending move.
    
    found_descent = false;
    trial_num_limit = set_num_attempts(s_set);
    
    for try_num = 1:trial_num_limit    
        % propose a move
        s_set_cdd = single_mutation(s_set);
        s_set_cdd = renumber_clusters(s_set_cdd);
        
        % accept only when strictly descending (a "zero-temperature" behavior)
        HSpos = costfun(s_set_cdd);
        delta_HS = HSpos - HS;
        if delta_HS < 0
            s_set = s_set_cdd;
            HS = HSpos;
            found_descent = true;
            return;
        end
    end
end

function trial_num_limit = set_num_attempts(s_set)
    % set a limit to the # non-descending trials before calling an end,
    % according to the current state (s_set)
    
    N = numel(s_set); % number of loci
    K = max(s_set); % number of clusters in the current state
    trial_num_limit = N * K; % number of possible moves from this state
end
