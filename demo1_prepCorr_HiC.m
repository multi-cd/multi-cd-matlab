% demo1_prepCorr_HiC.m
% 
% Demo script for the pre-processing of Hi-C data (M),
% which is first translated to the contact probability (P), 
% then to the auxiliary matrix of inverse distance scales (G) and 
% finally to the correlation matrix (C) to be used for Multi-CD.
% 
% - Input: Hi-C matrix M
% - Output: correlation matrix C
% ------------------------------------------------------------------------

% Copyright 2018 Min Hyeok Kim & Ji Hyun Bak


%% initialize 

clear;close all;clc;

setpaths; % add path to custom functions


%% load Hi-C data

disp('Pre-processing for Multi-CD');

% ========================================================================
% DATA INPUT:
% replace with the path to your own Hi-C data, and load as HiC_input.
% ------------------------------------------------------------------------
inputfilename = 'Data/HiC_test.mat';
if(exist(inputfilename,'file'))
    tempvar = load(inputfilename);
    HiC_input = tempvar.HiC_test;
else
    error('input data could not be loaded. check file path.');
end
% ========================================================================

disp(' ');
disp('input data loaded.');


%% prepare the Hi-C matrix M

disp(' ');
disp('=== Hi-C matrix M ===');

% ===== check input HiC data 

% remove all NaNs (replace with 0's)
idx_nan = sum(~isnan(HiC_input),2)==0; 
HiC_full = HiC_input(~idx_nan,~idx_nan);
HiC_full(isnan(HiC_full))=0;

% should be a non-negative matrix
if(any(HiC_full<0))
    error('input Hi-C should be non-negative.');
end

% should be either a symmetric matrix, or a triangular matrix
if(~isequal(HiC_full,HiC_full'))
    % try to symmetrize
    uppertri = triu(HiC_full,1); % above diagonal
    lowertri = tril(HiC_full,-1); % below diagonal
    if(~any(uppertri(:)>0))
        HiC_full = HiC_full + lowertri';
    elseif(~any(lowertri(:)>0))
        HiC_full = HiC_full + uppertri';
    else
        error('input Hi-C should be a symmetric matrix, or only contain upper/lower triangular elements.');
    end
end

disp('input Hi-C format OK.');
disp(['size of input matrix: ',num2str(size(HiC_full,1))]);

% ===== simple treatment

% remove empty rows/columns
idx_empty = all(HiC_full==0,2);
Mmat = HiC_full(~idx_empty,~idx_empty);

% detect size
N = size(Mmat,1); % size of Hi-C (number of genomic loci)

disp(['size after removing empty loci: N=',num2str(N)]);

% remove diagonal (self-contact)
Mmat(logical(eye(N))) = 0; 


% ===== plot M 

clf;
subplot(2,2,1)
imagesc(log(Mmat))
axis square
title('log (Hi-C matrix M)')
colormap(gca,'jet')
colorbar eastoutside
drawnow;


%% the contact probability matrix P

disp(' ');
disp('=== contact probability matrix P ===');

% === normalization

% set custom parameter
p1fix = 0.9; % nearest-neighbor contact probability

% normalize to fix the mean of nearest-neighbor contact probability
meanMdiag1 = mean(diag(Mmat,1)); % mean of +1 diagonal
Z = meanMdiag1/p1fix; 
Pmat = Mmat/Z; % such that mean of +1 diagonal = p1fix

disp('normalized to fix the mean of +1st diagonal to p1');
disp(['custom value: p1=',num2str(p1fix)]);


% === more treatment

% fill diagonal with NaN's (the diagonal will be fixed below with C=1)
Pmat(logical(eye(N))) = NaN;

% suppress occasional "outliers" with P>1
idx_highP = (Pmat(:)>1); % each pair is found twice (symmetric)
Pmat(idx_highP) = 1; % suppress down to 1

% check that the outliers are only a small fraction
numOutliers = sum(idx_highP)/2; % count # outliers (P>1)
numTotalPairs = N*(N-1)/2;
disp(['removed outliers (',num2str(numOutliers), '/',num2str(numTotalPairs),...
    '=',num2str(numOutliers/numTotalPairs,'%1.1e'),')']);


% === plot P

subplot(2,2,2)
imagesc(Pmat)
axis square
title('contact prob. matrix P')
colormap(gca,'jet')
colorbar eastoutside
drawnow;



%% the matrix G (Gamma) of inverse square distance scales

disp(' ');
disp('=== inverse square distance scales G ===');

% Each element of P is given as an integral over the inter-loci distance r
% from 0 to r_cutoff, with the elements of G (gamma) in the integrand.
% Here we change the variable of integration such that:
% 
% t = gamma * r^2
% 
% which is a dimensionless parameter.

% set custom parameter
r_cutoff = 1; % distance cutoff for inter-loci contact


% === prepare an inversion table

% the integrand
myfun = @(t) sqrt(2/pi)*(t.^2).*exp(-(1/2)*t.^2);

% set up t grid
dt = 0.05;
tgrid = (0:dt:5)'; % grid endpoints

% calculate the definite integrals from 0 to t
tmid = tgrid(1:end-1)+dt/2; % midpoints
ymid = myfun(tmid); % evaluate at midpoints
Yinteg = [0;cumsum(ymid*dt)]; % integrate using rectangle method


% === invert to get G from P

% prepare a non-redundant list of contact probability values
[P_list,~,ic] = unique(Pmat(:));

% invert P to the auxiliary variable t using a simple table search
t_list = table_search_YtoT(P_list,Yinteg,tgrid);

% transform to G, in the unit of inverse square distance
G_list = (t_list./r_cutoff).^2;

disp('inverted by table search. rearranging in matrix form...');

% rearrange back in matrix form (2/14/2019 update, with speedup)
Gmat_vectorized = G_list(ic); % map back using stored index correspondence
Gmat_vectorized(isnan(G_list)) = 0; % fill diagonals with 0's
Gmat = reshape(Gmat_vectorized,[N,N]); % rearrange to a matrix

disp('done.');


% === plot G

subplot(2,2,3)
imagesc(Gmat)
axis square
title('inv. distance scales G')
colormap(gca,'jet')
colorbar eastoutside
drawnow;


%% the correlation matrix C

disp(' ');
disp('=== correlation matrix C ===');

% === normalize to get the correlation matrix

% determine the normalizer
idx_lowertri = logical(tril(ones(N))-eye(N)); % lower triangular part, minus diagonal
gfix = median(Gmat(idx_lowertri)); % we used median (not the unique choice)
% note: this gfix corresponds to 1/(2*sigma_c) in the manuscript

% fix occasional cases where median is zero (when Hi-C is sparse)
if(gfix==0) 
    gfix = min(Gmat(Gmat>0)); % use the smallest positive value instead
end

% get the correlation matrix C
Cmat = 1 - gfix./Gmat; % by normalizing G
% Cmat = Cmat + Cmat'; % (bug fix 2/19/2019: already symmetric)
Cmat(logical(eye(N))) = 1; % fix diagonal (self-correlation)

% apply lower cutoff at -1
idx_lowC = (Cmat(:)<-1); % long tail due to finite # measurements in Hi-C
Cmat(idx_lowC)=-1; % cutoff to -1; this is still a strong anticorrelation

disp('correlation matrix ready.');

% === plot C

subplot(2,2,4)
imagesc(Cmat)
axis square
title('correlation matrix C')
colormap(gca,'jet')
colorbar eastoutside
drawnow;


% === optionally, save the correlation matrix 
% === to use as Multi-CD input

% outfilename = 'Data/Cmat_test.mat';
% save(outfilename,'Cmat');

