function [corrMat,S,X] = gen_corrMat_multiScale(N,c_num_list,g_mean_list,g_std_list,cont_option_list,numDraw)
% generate a correlation matrix with multi-scale structures

% 2018 Min Hyeok Kim & Ji Hyun Bak


%% draw states at each level

s_set0 = (1:N)'; % site indices
s_prev = s_set0;

S = [];
X = [];
x_list_all = zeros(N,1);

for lv = 1:numel(c_num_list)
    
    % unpack input parameters
    c_num = c_num_list(lv); % number of clusters
    g = struct('mean',g_mean_list(lv),'std',g_std_list(lv)); % clustering strength
    forceContinuous = logical(cont_option_list(lv)); % force continuous domains or not
    
    % generate states
    s_set = buildDomains(s_prev,c_num,forceContinuous); % build domains
    x_list = addNoise(s_set,g,numDraw); % add noise
    
    % pack 
    S.(['Level',num2str(lv)]) = s_set;
    X.(['Level',num2str(lv)]) = x_list;
    s_prev = s_set;
    
    % stack up multi-scale structures
    x_list_all = x_list_all + x_list;
    
end

% get a correlation matrix
corrMat = corr(x_list_all');

end

% =======================================================================

function s_set = buildDomains(s_set_lower,c_num,forceContinuous)

if(forceContinuous)
    % build continuous domains
    s_set = buildDomains_cont(s_set_lower,c_num);
else
    % allow non-continuous domains (random assignment)
    s_set = buildDomains_rand(s_set_lower,c_num);
end

end

function s_set = buildDomains_cont(s_set_lower,c_num)
% build domains such that each domain has a continuous territory
% along the linear index

% --- unpack input
N = numel(s_set_lower);
c_num_lower = numel(unique(s_set_lower));

% --- set up domain boundaries
while 1
    bound = randperm(c_num_lower-1,c_num-1);
    bound = sort(bound);
    if bound~=1
        break
    end
end
bound=bound';
bound(:,2) = [bound(2:end)-1;c_num_lower];
bound = [[1,bound(1)-1];bound];

% --- set up domain state vector
s_set = zeros(N,1);
if(isempty(s_set_lower)) % finest-level domains
    for s=1:size(bound,1)
        s_set(bound(s,1):bound(s,2))=s;
    end
    
else
    for i=1:size(bound,1)
        ix = [];
        for s = bound(i,1):bound(i,2)
            ix = [ix;find(s_set_lower==s)];
        end
        s_set(ix) = i;
    end

end

end

function s_set = buildDomains_rand(s_set_lower,c_num)
% build domains by random assignment of lower-level units

% --- unpack input
N = numel(s_set_lower);
c_num_lower = numel(unique(s_set_lower));

% --- domain assignment
s_set = zeros(N,1);
for s=1:c_num_lower
    ix = (s_set_lower==s);
    s_set(ix) = randperm(c_num,1);
end

end

function [x_list,g_r] = addNoise(s_set,g,numDraw)

% --- unpack input
N = numel(s_set); % number of sites
c_num = max(s_set); % number of domains

% --- generate random numbers from normal distribution
eta_r_temp = normrnd(0,1,c_num,numDraw);
eps_r = normrnd(0,1,N,numDraw);

% --- set up domain-level fluctuations
g_list = (g.mean) + (g.std)*rand(c_num,1); % strength of each cluster
eta_r = zeros(N,numDraw);
g_r = zeros(N,numDraw);
for s=1:c_num
    ix = find(s_set==s);
    eta_r(ix,:) = repmat(eta_r_temp(s,:),[size(ix,1),1]);
    g_r(ix,:) = g_list(s);
end

% --- get the "observed" state

% noise model [cf: JD Noh (2000), PRE]
myfun = @(g,eta,eps) (sqrt(g./(g+1)).*eta+1./sqrt(1+g).*eps);

% add up
x_list = myfun(g_r,eta_r,eps_r);

end
