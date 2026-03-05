clc;clear;
addpath 'user/MultipleTestingToolbox'

path ='/project/';

LEAP_ROI = readtable([path 'atlas_group.xlsx'],'Sheet','LEAP_243ROI');
reorderT = sortrows(LEAP_ROI,{'network_order','Group_Order_Buch','Number'},{'ascend','ascend','ascend'});

datapath = '/user/FC_asd_data/';
load('/user/T_sub.mat');
asd_idx = find(T.t1_diagnosis==2);
T_asd = T(asd_idx,:);
load('/user/sub1000_split.mat');

respath = '/user/RCCA_results_243ROI/';
load([respath 'CV_brain_permuted_1000replicates.mat']);
ROI_num = size(LEAP_ROI,1);
mask = flipud(tril(ones(ROI_num),-1));
for s = 1:size(T_asd,1)
    load([datapath strcat(num2str(T_asd.subjects(s)),'_corr.mat')]);
    vector_fc = fc(mask==1);
    allfc(s,:) = vector_fc';
end
for i = 1:length(Sbrain_permuted)
    trsub_fc = zscore(allfc(train_idx1(i,:),:));
    [r_val{i},p_val{i}] = corr(Sbrain_permuted{i},trsub_fc);
end
all_r = cat(3,r_val{:});

mean_loading = mean(all_r,3);
std_loading = std(all_r,0,3);

%% reshape 
fc_loading = zeros(ROI_num,ROI_num);
fc_loading_std = zeros(ROI_num,ROI_num);

for p = 1:size(mean_loading,1)
    fc_loading(mask==1) = mean_loading(p,:);
    fc_loading_all{p} = fc_loading + fliplr(flipud(fc_loading'));

    fc_loading_std(mask==1) = std_loading(p,:);
    fc_loading_std_all{p} = fc_loading_std + fliplr(flipud(fc_loading_std'));
end

save([respath 'mean_fc_std.mat'],'fc_loading_std_all');
save([respath 'mean_fc_loading.mat'],'fc_loading_all');

%% plot Fig. s2
load([respath 'mean_fc_loading.mat']);
ROI_num = size(fc_loading_all{1},1);

n = 256;%samples when calclating correlation
df = n-2;
mask1 = flipud(tril(ones(ROI_num),-1));
thr = 0.05;

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
    [cpval1,~,~] = fdr_BH(pval(mask1==1), thr);
    cpval_new1(mask1==1) = cpval1;
    cpval = cpval_new1 + fliplr(flipud(cpval_new1'))+fliplr(eye(ROI_num));
    cpval_sig(cpval<thr)=1;

    fc_loading_sig{cv} = fc_loading_all{cv}.*cpval_sig;
end

cv_num = 1;
fc_loading1 = fc_loading_sig{cv_num}; 

figure
set(gcf, 'unit', 'centimeters', 'position', [2 2 18 18]); 
ax3 = axes('Position', [0.13, 0.1, 0.7, 0.7]); 
imagesc(fc_loading1);%RRB
clim([-0.4 0.4]);
colormap(ax3,map);
% colorbar

title('CV1 - RRB related fc','FontSize',18)
set(gca, 'XTick', [], 'YTick', []);
hold on
%% network assignment
nLim = 26;
nDMN = 64;
nFPTC = 21;
nSN = 32;
nCOTC = 19;
nSMN = 24;
nAud = 11;
nVis = 35;
nCBM = 8;
nBS = 3;

xline([nLim+0.5 nLim+nDMN+0.5 nLim+nDMN+nFPTC+nSN+nCOTC+0.5],'-','LineWidth',1.2,'Color','k')
yline([nSMN+nAud+nVis+nCBM+nBS+0.5 nSMN+nAud+nVis+nCBM+nBS+nFPTC+nSN+nCOTC+0.5 nSMN+nAud+nVis+nCBM+nBS+nFPTC+nSN+nCOTC+nDMN+0.5],'-','LineWidth',1.2,'Color','k')

load('colorbar_brain.mat');

ax1 = axes('Position', [0.1, 0.1, 0.03, 0.7]);
imagesc([1:ROI_num]'); 
colors = [repmat(brain_map(10,:),nLim,1);repmat(brain_map(9,:),nDMN,1);repmat(brain_map(8,:),nFPTC,1);...
    repmat(brain_map(7,:),nSN,1);repmat(brain_map(6,:),nCOTC,1);repmat(brain_map(5,:),nSMN,1);...
    repmat(brain_map(4,:),nAud,1);repmat(brain_map(3,:),nVis,1);repmat(brain_map(2,:),nCBM,1);repmat(brain_map(1,:),nBS,1)];
colormap(ax1, colors);
set(gca, 'XTick', [], 'YTick', []);
axis tight;
set(gca, 'YDir', 'normal'); 

ax2 = axes('Position', [0.13, 0.07, 0.7, 0.03]);
imagesc([1:ROI_num]); 
colormap(ax2, colors);
set(gca, 'XTick', [], 'YTick', []);
axis tight;

