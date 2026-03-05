clc;clear;
addpath 'user/MultipleTestingToolbox'

path ='/project/';
LEAP_ROI = readtable([path 'atlas_group.xlsx'],'Sheet','LEAP_243ROI');
reorderT = sortrows(LEAP_ROI,{'network_order','Group_Order_Buch','Number'},{'ascend','ascend','ascend'});

load('/user/T_sub.mat');
asd_idx = find(T.t1_diagnosis==2);
T_asd = T(asd_idx,:);
load('/user/sub1000_split.mat');

load('colorbar(bbr).mat');
respath = '/user/RCCA_results_243ROI/';
load([respath 'mean_fc_loading.mat']);

ROI_num = size(fc_loading_all{1},1);

n = 256;%samples when calclating correlation
df = n-2;
mask1 = flipud(tril(ones(ROI_num),-1));
thr = 0.05;
% thr = 0.001;

for cv = 1:size(fc_loading_all,2)
    cur_r = fc_loading_all{:,cv};
    cpval_new1 = zeros(ROI_num,ROI_num);
    cpval_sig = zeros(ROI_num,ROI_num);

    for i = 1:size(cur_r,1)
        for j = 1:size(cur_r,1)
                r = cur_r(i,j);
                t = r * sqrt(df/(1-r^2));
                p = 2 * (1-tcdf(abs(t),df));
                pval(i,j) = p;
        end
    end
    [cpval1,~,~] = fdr_BH(pval(mask1==1),thr);
    cpval_new1(mask1==1) = cpval1;
    cpval = cpval_new1 + fliplr(flipud(cpval_new1'))+fliplr(eye(ROI_num));
    cpval_sig(cpval<thr)=1;

    fc_loading_sig{cv} = fc_loading_all{cv}.*cpval_sig;
end

%% merge_plot
cv_num = 1;
reverse_fc = flipud(fc_loading_sig{cv_num});
mapping = reorderT.Group_Order_Buch;
region_num = length(unique(reorderT.Group_Order_Buch));
newfc = zeros(region_num,region_num);

for i = 1:region_num
    for j = 1:region_num
        rows = find(mapping == i);
        cols = find(mapping == j);
        smallfc = reverse_fc(rows, cols);
        newfc(i, j) = mean(smallfc(smallfc~=0));
        newfc(j, i) = newfc(i, j); 
    end
end
fc_min = flipud(newfc);
fc_min(isnan(fc_min)) = 0;

%% plot
figure
set(gcf, 'unit', 'centimeters', 'position', [2 2 18 18]); 
ax3 = axes('Position', [0.16, 0.23, 0.7, 0.7]); 
%% fig2
imagesc(fc_min);%RRB
clim([-0.4 0.4]);
colormap(ax3,map);

% h = colorbar;
% h.FontSize = 14;
title('CV1 - RRB related fc','FontSize',18)
set(gca, 'XTick', [], 'YTick', []);

hold on
%% network assignment
netnum = 39;
nLim = 9;
nDMN = 8;
nFPTC = 3;
nSN = 4;
nCOTC = 3;
nSMN = 2;
nAud = 1;
nVis = 5;
nCBM = 1;
nBS = 3;

xline([nLim+0.5 nLim+nDMN+0.5 nLim+nDMN+nFPTC+nSN+nCOTC+0.5],'-','LineWidth',1.2,'Color','w')
yline([nSMN+nAud+nVis+nCBM+nBS+0.5 nSMN+nAud+nVis+nCBM+nBS+nFPTC+nSN+nCOTC+0.5 nSMN+nAud+nVis+nCBM+nBS+nFPTC+nSN+nCOTC+nDMN+0.5],'-','LineWidth',1.2,'Color','w')
load('colorbar_brain.mat')

ax1 = axes('Position', [0.12, 0.23, 0.03, 0.7]);
imagesc([1:netnum]'); 
colors = [repmat(brain_map(10,:),nLim,1);repmat(brain_map(9,:),nDMN,1);repmat(brain_map(8,:),nFPTC,1);...
    repmat(brain_map(7,:),nSN,1);repmat(brain_map(6,:),nCOTC,1);repmat(brain_map(5,:),nSMN,1);...
    repmat(brain_map(4,:),nAud,1);repmat(brain_map(3,:),nVis,1);repmat(brain_map(2,:),nCBM,1);repmat(brain_map(1,:),nBS,1)];
colormap(ax1, colors);

region_name = unique(reorderT.Group_Buch,'stable');
net_name = unique(reorderT.network_Buch,'stable');
set(gca, 'XTick',[], 'YTick',[1:netnum],'YTicklabel',region_name,'TickLength',[0,0],'FontSize',9);
axis tight;
set(gca, 'YDir', 'normal'); 

ax2 = axes('Position', [0.16, 0.19, 0.7, 0.03]);
imagesc([1:netnum]); 
colormap(ax2, colors);
set(gca, 'XTick', [nLim/2+0.5 nLim+nDMN/2+0.5 nLim+nDMN+nFPTC/2+0.5 nLim+nDMN+nFPTC+nSN/2+0.5 ...
    nLim+nDMN+nFPTC+nSN+nCOTC/2+0.5 nLim+nDMN+nFPTC+nSN+nCOTC+nSMN/2+0.5 nLim+nDMN+nFPTC+nSN+nCOTC+nSMN+nAud/2+0.5 ...
    nLim+nDMN+nFPTC+nSN+nCOTC+nSMN+nAud+nVis/2+0.5 nLim+nDMN+nFPTC+nSN+nCOTC+nSMN+nAud+nVis+nCBM/2+0.5 ...
    nLim+nDMN+nFPTC+nSN+nCOTC+nSMN+nAud+nVis+nCBM+nBS/2+0.5],'XTicklabel',net_name, 'XTickLabelRotation', 45,...
    'YTick', [],'TickLength',[0,0],'FontSize',10);
axis tight;
