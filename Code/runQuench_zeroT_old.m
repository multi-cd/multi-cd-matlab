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

%% unpack input

% starting state
s_set = s_set_start;
HS = costfun(s_set);

% set a limit to the # non-descending trials before calling an end
trial_num_limit = numel(s_set)*max(s_set); % number of possible moves

% unpack input options
cnt_acp_limit = quenchOptions.cnt_acp_limit;

% termination status
term_status = 0; % 0 means running


%% quickly descend along the energy landscape

cnt_acp = 0; % initialize acceptance number 
try_num = 0; % initialize trial number 

HS_list = NaN(1+cnt_acp_limit,1); % track descending H values
HS_list(1+cnt_acp) = HS;

while 1
    
    try_num = try_num+1;
    s_set_cdd = single_mutation(s_set);
    s_set_cdd = renumber_clusters(s_set_cdd);
    
    HSpos = costfun(s_set_cdd);
    delta_HS = HSpos-HS;
    
    if delta_HS<0 % a "zero-temperature" behavior
        
        cnt_acp = cnt_acp+1;
        s_set = s_set_cdd;
        
        HS = HSpos;
        HS_list(1+cnt_acp) = HS;
        
        % Stopping condition 1: too many descents (too far from local min)
        if cnt_acp>=cnt_acp_limit 
            term_status = 1; % a bad sign (premature end of program)
            break;
        end
        
        % reset try_num and update trial limit 
        try_num=0; 
        trial_num_limit = numel(s_set)*max(s_set); % update # possible moves
        
    end
    
    % Stopping condition 2: no more descents after "enough" trials
    if try_num>trial_num_limit 
        term_status = 2;
        HS_list = HS_list(1:cnt_acp); % trim list
        break;
    end
    
end

end
