clc;clear;
addpath 'user/MultipleTestingToolbox'

path ='/project/';
LEAP_ROI = readtable([path 'atlas_group.xlsx'],'Sheet','LEAP_243ROI');
reorderT = sortrows(LEAP_ROI,{'network_order','Group_Order_Buch','Number'},{'ascend','ascend','ascend'});

load('/user/T_sub.mat');
asd_idx = find(T.t1_diagnosis==2);
T_asd = T(asd_idx,:);
load('/user/sub1000_split.mat');

respath = '/user/RCCA_results_243ROI/';
load([respath 'mean_fc_loading.mat']);

ROI_num = size(fc_loading_all{1},1);

n = 256;%samples when calclating correlation
df = n-2;
mask1 = flipud(tril(ones(ROI_num),-1));
% thr = 0.05;
thr = 0.001;

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

cv_num = 1;
fc_loading1 = fc_loading_sig{cv_num};

reverse_fc = flipud(fc_loading1);
mapping = reorderT.Group_Order_Buch;
region_num = length(unique(reorderT.Group_Order_Buch));
newfc = zeros(region_num,region_num);

for i = 1:region_num
    for j = i:region_num
        rows = find(mapping == i);
        cols = find(mapping == j);
        smallfc = reverse_fc(rows, cols);
        newfc(i, j) = mean(smallfc(smallfc~=0));
        newfc(j, i) = newfc(i, j); 
    end
end

mask2 = ones(region_num);
mask2(eye(region_num) == 1) = 0;
fc_dis = newfc .* mask2;
fc_dis(isnan(fc_dis)) = 0;

%% save for Rplot
region39 = readtable([path '39region_label.xlsx'],'Sheet','243ROI');

zeroRows = all(fc_dis == 0, 2);
zero_idx = find(zeroRows);
fc_dis(zero_idx,:) = [];
fc_dis(:,zero_idx) = [];
region39(zero_idx,:) = [];

save([respath 'zero_idx_rrb_0.01.mat'],"zero_idx");
save([respath 'rrb_circle_0.01.mat'],"fc_dis");
writetable(region39,[respath 'rrb_ccT_0.01.xlsx']);
