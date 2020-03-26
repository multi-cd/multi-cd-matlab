function [s_best, H_best] = call_MultiCD(Cmat, lambda, opts, num_repeat)
    % Iterations for Multi-CD
    
    % set handle for the cost function, given the correlation matrix Cmat
    costfun = @(svec) HS_calculation_all(Cmat,svec,lambda);
    
    N = size(Cmat,1); % total number of loci  
    
    % run simulated annealing starting from multiple initializations
    s_set_list = zeros(N, num_repeat);
    H_list = zeros(num_repeat, 1);
    for nr = 1:num_repeat
        disp(['Run ', num2str(nr)]);
        [s_set, HS, ~] = runSimulatedAnnealing(costfun, N, opts);    
        s_set_list(:, nr) = s_set;
        H_list(nr) = HS;
    end
    
    % report the best configuration
    [H_best, i_best] = min(H_list);
    s_best = s_set_list(:, i_best);
    
end
