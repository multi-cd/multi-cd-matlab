% demo2_MultiCD.m
% 
% Demo script for Multi-CD at a specified parameter lambda, 
% which runs a simulated annealing to find 
% the chromatin domain (CD) solution at the corresponding scale.
% 
% - Input: correlation matrix (C) & Multi-CD parameter (lambda)
% - Output: CD solution (s), as a vector of domain indices
% ------------------------------------------------------------------------

% Copyright 2018 Min Hyeok Kim & Ji Hyun Bak


%% initialize 

clear;clc;close all;

setpaths; % add path to custom functions

% ====== set Multi-CD parameter lambda ===================================
lambda = 1; % larger lambda prefers simpler solutions (less # clusters)
% ========================================================================


%% load HiC correlation matrix 

useRealData = false; % if false, generate a synthetic correlation matrix

% specify input data location (if applicable)
Cmat_filename = 'Data/Cmat_test.mat';

if(useRealData && exist(Cmat_filename,'file')) % --- load from existing file
    
    tempvar = load(Cmat_filename);
    Cmat_full = tempvar.Cmat;
    
    % choose subset
    idx_sub = 150+(1:50); % select subset range
    Cmat = Cmat_full(idx_sub,idx_sub);
    
else % --- generate a synthetic model
    
    % parameter setting
    Ngen = 50; % system size
    config = 100; % initialization number
    c_num_list = [10 5 3]; % number of clusters at each scale
    g_mean_list = [3 2 1]; % mean clustering strength
    g_std_list = [3 2 1]; % stdev of clustering strength
    cont_option_list = [true true false]; % continuous domains or not
    
    % sample a copy of correlation matrix
    [Cmat,~,~] = gen_corrMat_multiScale(Ngen,c_num_list,g_mean_list,g_std_list,cont_option_list,config);
    
end


% show input correlation matrix 
clf;
subplot(2,2,1)
imagesc(Cmat)
colormap(gca,'jet')
axis square
title('data C')
colorbar eastoutside
drawnow;

% set handle for the cost function, given the correlation matrix Cmat
costfun = @(svec) HS_calculation_all(Cmat,svec,lambda); 


%% initialization

disp(['Multi-CD at lambda=',num2str(lambda)]);

disp(' ');
disp('=== Initialization ===');

N = size(Cmat,1); % total number of loci

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

disp(['initial temp. = ',num2str(T_init,'%1.1f')]);


%% simulated annealing

% ===== parameter setup

% parameters for temperature schedule in SA
c_cool = 0.8; % cooling factor
% -- termination conditions for SA: terminate at whatever is reached first
maxIterSA = 100; % max number of steps in SA
T_stop = 0.1; % target temperature for escaping the SA

% MCMC options
mcmcOptions = [];
mcmcOptions.minIterMCMC = 100;
mcmcOptions.maxIterMCMC = 5000;
mcmcOptions.numSampAtEq = 500;


% ==== simulated annealing 

disp(' ');
disp('=== Simulated Annealing ===');

% initialize
T = T_init;
s_set = s_init;
disp(['Initial T=',num2str(T,'%1.1f')]);

HS_track = NaN(maxIterSA,3); % [best mean std]
T_track = NaN(maxIterSA,3);

tic;
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

    disp(['Step #',num2str(r),...
        ', T=',num2str(T,'%1.1f'),', H=',num2str(HS,'%1.1f')]);
    
    % plot best solution
    Bmat = bsxfun(@eq,s_set,s_set');
    subplot(2,2,4)
    imagesc(Bmat);
    colormap(gca,[1 1 1;0 0 0])
    axis square
    title(['best at T=',num2str(T,3)])
    
    % plot ensemble
    Bmat_ens = mean(bsxfun(@eq,permute(s_list,[1 3 2]),permute(s_list,[3 1 2])),3);
    subplot(2,2,3)
    imagesc(Bmat_ens);
    colormap(gca,bsxfun(@times,(1:-0.01:0)',[1 1 1]))
    axis square
    title(['ensemble at T=',num2str(T,3)])
    
    % track cost function
    subplot(2,2,2)
    plot(T_track(1:r),HS_track(1:r,1),'k.-','markersize',10,'linewidth',1)
    set(gca,'xscale','log')
    xlim([T_stop T_init])
    xlabel('temperature T (log scale)')
    ylabel('cost function')
    title(['annealing #',num2str(r),', T=',num2str(T,3)])
    
    drawnow;
    
    if(quit_cond_SA)
        disp('* single-domain solution reached -- force quitting SA.');
        break
    end
    
    % update temperature
    T = T*c_cool; % cooling by a constant factor
    if(T<T_stop) 
        break
    end
    
end

disp('done.');


%% Zero temperature sampling (Gradient descent method) 

disp(' ');
disp('=== Zero Temperature Quenching ===');

% quench options
quenchOptions = [];
quenchOptions.cnt_acp_limit = 1e5;

% zero temperature quenching
[s_set,HS,HS_list,term_status] = runQuench_zeroT(costfun,s_set,quenchOptions);
if(term_status==1)
    disp(' - bad sign: descending too far. trial limit reached.');
    %warning('maybe SA was not enough. check temperature schedule.');
elseif(term_status==2)
    disp('Minimum reached successfully.');
end

telaps = toc;
disp(['Total elapsed time: ',num2str(telaps,'%1.1f'),'sec']);

% plot solution
Bmat = bsxfun(@eq,s_set(:),s_set(:)');
subplot(2,2,4)
imagesc(Bmat);
colormap(gca,[1 1 1;0 0 0])
axis square
title('final')
drawnow;
