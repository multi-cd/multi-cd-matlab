function [s_init, K_init] = initial_state_generation(N)
    % generate an initial state for simulated annealing.

    % randomly select the number of clusters
    % (to allow sampling of highly clustered or highly fragmented solutions)
    K_init = randperm(N,1);

    % randomly assign each loci to a cluster
    s_init = randsample(K_init,N,true);
    s_init = renumber_clusters(s_init);

    K_init = max(s_init); % fix in case some domains were never assigned
end
