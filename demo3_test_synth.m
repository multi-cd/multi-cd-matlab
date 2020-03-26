% demo3_test_synth.m
% 
% Demo script for Multi-CD on a synthetic correlation matrix.
% This is just a rearranged combination of the other two demo scripts,
% `demo2_MultiCD.m` and `Code/gen_model/demo_gen_corrMat.m`.
% ------------------------------------------------------------------------

% Copyright 2020 Ji Hyun Bak


%% initialize 

clear;clc;close all;

setpaths; % add path to custom functions



%% generate a model (synthetic) correlation matrix 

% parameters for correlation matrix
Ngen = 50; % system size
config = 100; % initialization number
c_num_list = [10 2]; % number of clusters at each scale
g_mean_list = [2 2]; % mean clustering strength
g_std_list = [1 1]; % stdev of clustering strength
cont_option_list = [true false]; % continuous domains or not

% sample a copy of correlation matrix
[Cmat, S_gen, ~] = gen_corrMat_multiScale(Ngen,c_num_list,g_mean_list,g_std_list,cont_option_list,config);


% ===== hidden domain states

% unpack states
s_set1 = S_gen.Level1;
s_set2 = S_gen.Level2;

% plot correlation matrices at individual levels

clf;
% set(gcf, 'position', [134 310 1118 389])
colormap(flip(gray))

subplot(2,5,[1:2 6:7])
imagesc(Cmat)
axis square
title('correlation matrix')
drawnow;

subplot(2,5,3)
biM = bsxfun(@eq,s_set1,s_set1');
imagesc(biM)
axis square
title('Level 1 (true)')

subplot(2,5,4)
biM = bsxfun(@eq,s_set2,s_set2');
imagesc(biM)
axis square
title('Level 2 (true)')



%% use Multi-CD to find domain solutions

% ===== parameter setup

opts = [];

% parameters for temperature schedule in SA
opts.SA.c_cool = 0.8; % cooling factor
% -- termination conditions for SA: terminate at whatever is reached first
opts.SA.maxIterSA = 1000; % max number of steps in SA
opts.SA.T_stop = 0.01; % target temperature for escaping the SA

% MCMC options
opts.MCMC.minIterMCMC = 100;
opts.MCMC.maxIterMCMC = 1000;
opts.MCMC.numSampAtEq = 500;

% quench options
opts.quench.cnt_acp_limit = 1e5;

% other options
opts.talkative = false;


% ===== Multi-CD loop

lambda_list = [0, 5, 10];
num_repeat = 10; % number of repeated runs at each lambda

for nl = 1:numel(lambda_list)
    
    lambda = lambda_list(nl);
    disp(['Multi-CD at lambda=',num2str(lambda)]);
    
    % find best domain solution using Multi-CD
    [s_set, H_value] = call_MultiCD(Cmat, lambda, opts, num_repeat);

    % plot solution
    Bmat = bsxfun(@eq,s_set(:),s_set(:)');
    subplot(2,5, 7 + nl)
    imagesc(Bmat);
    colormap(gca,[[1 1 1];[0 0 0]])
    axis square
    title(['\lambda = ', num2str(lambda)])
    drawnow;

end

set(findall(gcf, '-property', 'fontsize'), 'fontsize', 14)

% figname = 'fig_test_synth.eps';
% print(figname, '-depsc')
