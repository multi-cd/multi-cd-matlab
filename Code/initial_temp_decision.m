function T_init = initial_temp_decision(costfun,s_init,varargin)
% INITIAL_TEMP_DECISION: Determines a suitable initial temperature
% for simulation annealing, given an initial state and the cost function.
%
% Input:
%   - costfun: a function handle for the cost function, with one argument
%              that is the domain state vector -- e.g., costfun(s)
%   - s_init: initial state (N-vector)
%   - const_a: [OPTIONAL] acceptance probability for the "worst" move
%   - T_init_range: [OPTIONAL] a range [T_min, T_max] to constrain T_init
%
% Output:
%   - T_init: initial temperature (single number)
% ------------------------------------------------------------------------

% Copyright 2018 Min Hyeok Kim & Ji Hyun Bak

%% unpack input

% specify the desired acceptance probability for the "worst" move
if (nargin > 2)
    const_a = varargin{1}; % may be passed as input
else
    const_a = 0.5; % default value
end

if (nargin > 3)
    T_init_range = varargin{2};
else
    T_init_range = []; % unconstrained by default
end

% unpack dimensions
N = size(s_init,1); % number of loci
Kmax = max(s_init); % number of existing clusters

% NOTE: the number of all possible single mutations is N*Kmax.



%% scan through all possible single-mutation moves

HS = costfun(s_init); % initial value of the cost function

dH_list = -Inf(N,Kmax);

for j = 1:N

    s_j = s_init(j); % locus to mutate

    for s = 1:Kmax

        % ----- make a single mutation ------------
        s_mut = s; % move to the selected domain
        if(s==s_j)
            s_mut = Kmax+1; % if self, create a new domain
        end
        s_set_cdd = s_init;
        s_set_cdd(j) = s_mut; % apply the mutation
        % ------------------------------------------

        s_set_cdd = renumber_clusters(s_set_cdd); % renumber

        HSpos = costfun(s_set_cdd);
        deltaHS = HSpos-HS; % the "energy" difference
        dH_list(j,s) = deltaHS;

    end
end


% calculate the temperature that corresponds to the
% "worst" acceptance probability specified above

T_init = -max(dH_list(:))/log(const_a);

% NOTE: this line was in SimulatedAnnealing_accelerate
% optionally constrain within a range
if ~isempty(T_init_range) && numel(T_init_range) == 2
    T_min = T_init_range(1);
    T_max = T_init_range(2);
    T_init = min(T_max, max(T_min, T_init));
end

end
