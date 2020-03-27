% demo0_corrMat_mixed.m
    
%% groups

n_box = 20;
box = ones(n_box, 1);

s_set11 = [1*box; 2*box; 3*box; 4*box; 5*box; 5*box];
s_set12 = [1*box; 1*box; 2*box; 3*box; 3*box; 3*box];

s_set21 = [1*box; 1*box; 2*box; 3*box; 4*box; 4*box];
s_set22 = [1*box; 2*box; 3*box; 4*box; 4*box; 4*box];

scenario1 = {s_set11, s_set12};
scenario2 = {s_set21, s_set22};

scenario0 = scenario1; % equivalent

scenarios = {scenario1, scenario2, scenario0};

corr_mats = cell(numel(scenarios), 1);

N = numel(s_set11);

% unpack input parameters
% c_num = c_num_list(lv); % number of clusters
% g = struct('mean',g_mean_list(lv),'std',g_std_list(lv)); % clustering strength
% forceContinuous = logical(cont_option_list(lv)); % force continuous domains or not

for np = 1:numel(scenarios)
    
    scenario = scenarios{np};
    x_list_all = zeros(N,1);

    for ns = 1:numel(scenario)
        
        s_set = scenario{ns};
        
        g = struct('mean', 3, 'std', 0);
        numDraw = 500;
        
        % generate states
        % s_set = buildDomains(s_prev,c_num,forceContinuous); % build domains
        x_list = addNoise(s_set,g,numDraw); % add noise
        
        % pack
        % S.(['Level',num2str(lv)]) = s_set;
        % X.(['Level',num2str(lv)]) = x_list;
        % s_prev = s_set;
        
        % stack up multi-scale structures
        x_list_all = x_list_all + x_list;
        
    end
    
    % get a correlation matrix
    corrMat = corr(x_list_all');
    corr_mats{np} = corrMat;
    
end



%% plot

clf;

color0 = [0.7 0.7 0.1; 0.3 0.3 0.05];
color1 = [0.1 0.2 0.7];
color2 = [0.7 0.2 0.1];

subplot(3,3,2)
biM = bsxfun(@eq,s_set11,s_set11') + bsxfun(@eq,s_set12,s_set12');
imagesc(biM)
colormap(gca, [1 1 1; color0])
axis square
title('population 0')
set(gca, 'xtick', [], 'ytick', [])

subplot(3,3,4)
biM = bsxfun(@eq,s_set11,s_set11');
imagesc(biM)
colormap(gca, [1 1 1; color1])
axis square
title('population 1-1')
set(gca, 'xtick', [], 'ytick', [])

subplot(3,3,5)
biM = bsxfun(@eq,s_set12,s_set12');
imagesc(biM)
colormap(gca, [1 1 1; color1])
axis square
title('population 1-2')
set(gca, 'xtick', [], 'ytick', [])

subplot(3,3,7)
biM = bsxfun(@eq,s_set21,s_set21');
imagesc(biM)
colormap(gca, [1 1 1; color2])
axis square
title('population 2-1')
set(gca, 'xtick', [], 'ytick', [])

subplot(3,3,8)
biM = bsxfun(@eq,s_set22,s_set22');
imagesc(biM)
colormap(gca, [1 1 1; color2])
axis square
title('population 2-2')
set(gca, 'xtick', [], 'ytick', [])


subplot(3,3,3)
corrMat = corr(x_list_all'); % re-sample
imagesc(corrMat)
axis square
colormap(gca, flip(gray))
title('correlation matrix 0')
set(gca, 'xtick', [], 'ytick', [])

subplot(3,3,6)
imagesc(corr_mats{1})
axis square
colormap(gca, flip(gray))
title('correlation matrix 1')
set(gca, 'xtick', [], 'ytick', [])

subplot(3,3,9)
imagesc(corr_mats{2})
axis square
colormap(gca, flip(gray))
title('correlation matrix 2')
set(gca, 'xtick', [], 'ytick', [])


set(findall(gcf, '-property', 'fontsize'), 'fontsize', 14)

figname = 'fig_mixed_population';
print(figname, '-depsc')
