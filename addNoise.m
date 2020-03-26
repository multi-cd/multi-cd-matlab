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
