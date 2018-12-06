% demo_gen_corrMat.m
% 
% Demo script for the auxiliary function gen_corrMat_multiScale, 
% which generates a correlation matrix with multi-scale domains 
% at specified distributions of clustering strengths and noise.
% ------------------------------------------------------------------------

% Copyright 2018 Min Hyeok Kim & Ji Hyun Bak

%% initialize

clear;close all;clc;

setpaths; % add path to custom functions

% parameter setting 
N = 50; % system size
config = 1000; % initialization number
c_num_list = [10 5 3]; % number of clusters
g_mean_list = [4 2 1]; % mean clustering strength
g_std_list = [2 2 1]; % stdev of clustering strength
cont_option_list = [true true false]; % continuous domains or not


%% generate correlation matrix with multi-scale domains

% sample a copy of correlation matrix
[corrA,S,X] = gen_corrMat_multiScale(N,c_num_list,g_mean_list,g_std_list,cont_option_list,config);


% ===== hidden domain states

% unpack states
s_set1 = S.Level1;
s_set2 = S.Level2;
s_set3 = S.Level3;

% plot correlation matrices at individual levels

clf;
colormap(flip(gray))

subplot(4,3,1)
biM = bsxfun(@eq,s_set1,s_set1');
imagesc(biM)
axis square
title('Level 1 (true)')

subplot(4,3,2)
biM = bsxfun(@eq,s_set2,s_set2');
imagesc(biM)
axis square
title('Level 2 (true)')

subplot(4,3,3)
biM = bsxfun(@eq,s_set3,s_set3');
imagesc(biM)
axis square
title('Level 3 (true)')

drawnow;


% ===== "observed" states with noise

% unpack states
x_list1 = X.Level1;
x_list2 = X.Level2;
x_list3 = X.Level3;

% plot correlation matrices at individual levels

subplot(4,3,3+1)
corr1=corr(x_list1');
imagesc(corr1)
axis square
title('Level 1 (obs)')

subplot(4,3,3+2)
corr2=corr(x_list2');
imagesc(corr2)
axis square
title('Level 2 (obs)')

subplot(4,3,3+3)
corr3=corr(x_list3');
imagesc(corr3)
axis square
title('Level 3 (obs)')

drawnow;


% ===== plot the combined correlation matrix

subplot(4,3,7:12)
imagesc(corrA)
axis square
title('correlation matrix')
drawnow;

