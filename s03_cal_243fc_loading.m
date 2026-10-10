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

thr = 0.05;
for i = 1:length(Sbrain_permuted)
    trsub_fc = zscore(allfc(train_idx1(i,:),:));
    [r_val{i},p_val{i}] = corr(Sbrain_permuted{i},trsub_fc);

    [c_pval1{i}, ~, ~] = fdr_BH(p_val{i},thr);   
    c_pval{i} = reshape(c_pval1{i},3,[]);

    cp_mask = zeros(size(c_pval{i}));
    cp_mask(c_pval{i}<thr) = 1;
    cp_mask_com{i} = cp_mask;
end
cp_mask_sum = sum(cat(3,cp_mask_com{:}),3);
save([respath 'cp_mask_sum0.05.mat'],'cp_mask_sum');

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

